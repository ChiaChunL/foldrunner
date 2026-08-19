"""Client for the MMseqs2 alignment web services.

The ColabFold and Protenix endpoints run the same server, so one client covers
both and either can stand in when the other is busy. Requests carry a batch of
sequences; the reply packs one alignment block per query into a single file,
separated by null bytes and numbered from 101.

Two behaviours of these services need handling rather than merely reporting.
A throttled request answers with an HTML page under a normal status code, so the
archive check has to happen before unpacking or the failure surfaces as a
corrupt-file error far from its cause. And the whole service is a shared
academic resource sized for a few thousand alignments a day across all users,
so a panel of any size has to back off rather than retry immediately.
"""

from __future__ import annotations

import io
import random
import tarfile
import time
from dataclasses import dataclass, field

from foldrunner import __version__
from foldrunner.msa.a3m import species_id
from foldrunner.msa.backends.base import MsaResult, RateLimited, SearchError
from foldrunner.msa.cache import Recipe, looks_like_archive

PROTENIX_HOST = "https://protenix-server.com/api/msa"
COLABFOLD_HOST = "https://api.colabfold.com"

# Query records are numbered from 101 by the service; the reply uses the same
# numbers to label each block.
FIRST_QUERY_ID = 101

# The two services return different archives, verified against both live
# endpoints rather than assumed to match.
#
# ColabFold names its members after the databases searched. Protenix numbers
# them by query and adds uniref_tax.m8, a table joining each hit accession to a
# taxon — which is the only place a per-sequence alignment can get its taxonomy
# from, and therefore the only route to pairing without a second search.
PAIRED_MEMBER = "pair.a3m"
COLABFOLD_MEMBERS = ("uniref.a3m", "bfd.mgnify30.metaeuk30.smag30.a3m")
TAXONOMY_MEMBER = "uniref_tax.m8"

# Only ColabFold serves ticket/pair; the Protenix host answers it with 404, so
# the two are not interchangeable for pairing however alike they look otherwise.
PAIRING_HOSTS = (COLABFOLD_HOST,)

RETRY_STATUSES = frozenset({"UNKNOWN", "RATELIMIT"})
WAIT_STATUSES = frozenset({"UNKNOWN", "RUNNING", "PENDING"})


def split_blocks(text: str) -> dict[int, str]:
    """Separate a multi-query reply into one alignment per query number.

    A null byte marks the end of a block; the header that follows names the
    query it belongs to.
    """
    blocks: dict[int, list[str]] = {}
    current: int | None = None
    expect_header = True
    for line in text.splitlines(keepends=True):
        if "\x00" in line:
            line = line.replace("\x00", "")
            expect_header = True
        if line.startswith(">") and expect_header:
            try:
                current = int(line[1:].strip().split()[0])
            except (ValueError, IndexError) as error:
                raise SearchError(f"unexpected block header {line!r}") from error
            expect_header = False
            blocks.setdefault(current, [])
            continue
        if current is not None:
            blocks[current].append(line)
    return {number: "".join(lines) for number, lines in blocks.items()}


def pairable(a3m: str) -> str:
    """Keep the query and every hit whose organism is known.

    Pairing groups by organism, so a hit without one can never take part. The
    metagenomic databases that give an unpaired alignment its depth carry no
    taxonomy at all, which is why the two halves differ in size.
    """
    lines = a3m.splitlines(keepends=True)
    out: list[str] = []
    keep = False
    for index, line in enumerate(lines):
        if line.startswith(">"):
            keep = index == 0 or species_id(line) is not None
        if keep:
            out.append(line)
    return "".join(out)


def parse_taxonomy(text: str) -> dict[str, str]:
    """Read the accession-to-taxon table the Protenix service ships.

    Columns are query id, accession, taxon id, species name, lineage. Only the
    middle two matter here.
    """
    mapping: dict[str, str] = {}
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3 and fields[1]:
            mapping.setdefault(fields[1], fields[2])
    return mapping


def annotate_taxonomy(a3m: str, taxonomy: dict[str, str]) -> str:
    """Fold taxon ids into a3m headers, as the pairing code expects to find them.

    Hits come back as bare accessions, and pairing needs the organism. Rewriting
    them into ``UniRef100_<accession>_<taxid>/`` puts the taxonomy where the
    parser already looks, so a per-sequence search becomes pairable locally
    without asking the service to pair anything.
    """
    out = []
    for line in a3m.splitlines(keepends=True):
        if not line.startswith(">"):
            out.append(line)
            continue
        rest = line[1:]
        accession = rest.split("\t")[0].split()[0] if rest.strip() else ""
        taxon = taxonomy.get(accession)
        if taxon and not accession.startswith("UniRef"):
            out.append(f">UniRef100_{accession}_{taxon}/" + rest[len(accession):])
        else:
            out.append(line)
    return "".join(out)


def _rewrite_query_header(body: str, sequence: str) -> str:
    """Restore a readable query header on a block.

    The service strips the submitted header and the block starts straight at the
    query row, which downstream parsers expect to find labelled.
    """
    if body.startswith(">"):
        return body
    return f">query\n{body}" if body else f">query\n{sequence}\n"


@dataclass
class WebBackend:
    """Fetches alignments from one of the MMseqs2 services."""

    host: str = PROTENIX_HOST
    fallback: str | None = COLABFOLD_HOST
    contact: str = ""
    use_env: bool = True
    pairing_strategy: str = "greedy"
    max_attempts: int = 5
    poll_seconds: float = 5.0
    timeout: float = 60.0
    session: object | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        if not self.contact:
            # The services ask every client to identify itself and warn when one
            # does not; the warning is documented as becoming an error.
            raise SearchError(
                "set contact to an email address: the alignment services require "
                "callers to identify themselves"
            )

    @property
    def user_agent(self) -> str:
        return f"foldrunner/{__version__} {self.contact}"

    @property
    def recipe(self) -> Recipe:
        databases = ("uniref30", "colabfold_envdb") if self.use_env else ("uniref30",)
        backend = "protenix-api" if self.host == PROTENIX_HOST else "colabfold-api"
        return Recipe(
            backend=backend,
            databases=databases,
            tool_version="mmseqs2-service",
            params=(("pairing_strategy", self.pairing_strategy),),
        )

    def supports_pairing(self) -> bool:
        return self.host in PAIRING_HOSTS

    def _members_for(self, host: str, paired: bool) -> tuple[str, ...]:
        if paired:
            return (PAIRED_MEMBER,)
        # Protenix numbers its output per query; ColabFold names it per database.
        return COLABFOLD_MEMBERS + tuple(f"{index}.a3m" for index in range(64))

    def _mode(self, paired: bool) -> str:
        if not paired:
            return "env" if self.use_env else ""
        mode = "pairgreedy" if self.pairing_strategy == "greedy" else "paircomplete"
        return f"{mode}-env" if self.use_env else mode

    def _client(self):
        if self.session is not None:
            return self.session
        import requests

        return requests

    def _jitter(self, base: float) -> float:
        """Spread retries so a batch of clients does not resynchronise on the service."""
        return base + random.uniform(0, base * 0.6)

    def _sleep(self, attempt: int) -> None:
        time.sleep(self._jitter(min(60.0, self.poll_seconds * (2**attempt))))

    def _submit(self, host: str, sequences: list[str], mode: str) -> dict:
        query = "".join(
            f">{FIRST_QUERY_ID + index}\n{sequence}\n" for index, sequence in enumerate(sequences)
        )
        endpoint = "ticket/pair" if mode.startswith("pair") else "ticket/msa"
        response = self._client().post(
            f"{host}/{endpoint}",
            data={"q": query, "mode": mode},
            timeout=self.timeout,
            headers={"User-Agent": self.user_agent},
        )
        try:
            return response.json()
        except ValueError as error:
            raise RateLimited(
                f"{host} answered with a non-JSON body, which is how it reports "
                "throttling; wait a few minutes and retry"
            ) from error

    def _poll(self, host: str, ticket: str) -> dict:
        state = {"status": "PENDING"}
        waited = 0.0
        while state["status"] in WAIT_STATUSES:
            time.sleep(self._jitter(self.poll_seconds))
            waited += self.poll_seconds
            response = self._client().get(
                f"{host}/ticket/{ticket}",
                timeout=self.timeout,
                headers={"User-Agent": self.user_agent},
            )
            try:
                state = response.json()
            except ValueError as error:
                raise RateLimited(f"{host} stopped answering with JSON while polling") from error
        return state

    def _download(self, host: str, ticket: str) -> dict[str, str]:
        response = self._client().get(
            f"{host}/result/download/{ticket}",
            timeout=self.timeout * 2,
            headers={"User-Agent": self.user_agent},
        )
        blob = response.content
        if not looks_like_archive(blob):
            raise RateLimited(
                f"{host} returned a {len(blob)} byte body that is not an archive. "
                "This is how the service reports throttling, not a damaged download; "
                "the same request usually succeeds after a few minutes."
            )
        members: dict[str, str] = {}
        with tarfile.open(fileobj=io.BytesIO(blob)) as archive:
            for member in archive.getmembers():
                handle = archive.extractfile(member)
                if handle is not None:
                    members[member.name] = handle.read().decode()
        return members

    def _run(self, host: str, sequences: list[str], mode: str) -> dict[str, str]:
        for attempt in range(self.max_attempts):
            state = self._submit(host, sequences, mode)
            status = state.get("status", "UNKNOWN")
            if status in RETRY_STATUSES:
                self._sleep(attempt)
                continue
            if status == "MAINTENANCE":
                raise SearchError(f"{host} is under maintenance")
            if status == "ERROR":
                raise SearchError(f"{host} rejected the batch; check the sequences are protein")
            final = self._poll(host, state["id"])
            if final.get("status") == "COMPLETE":
                return self._download(host, state["id"])
            if final.get("status") == "ERROR":
                raise SearchError(f"{host} failed while running the search")
            self._sleep(attempt)
        raise RateLimited(f"{host} did not accept the batch after {self.max_attempts} attempts")

    def fetch(self, sequences: list[str], *, paired: bool = False) -> list[MsaResult]:
        """Search a batch, falling back to the second host if the first refuses.

        Asking for the pairing half does not call the pairing endpoint. That
        endpoint pairs server-side across whatever was submitted together, so it
        answers a single sequence with nothing at all and answers a complex with
        rows that carry no taxonomy — usable only for the exact complex it was
        asked about, which puts the cost back at one search per pair.

        A search that carries taxonomy can be paired locally instead, for any
        combination, so the pairing half is derived from the same hits: those
        whose organism is known.
        """
        if paired:
            return [
                MsaResult(sequence=r.sequence, paired=pairable(r.unpaired or ""))
                for r in self.fetch(sequences, paired=False)
            ]

        unique: list[str] = []
        for sequence in sequences:
            if sequence not in unique:
                unique.append(sequence)

        mode = self._mode(paired)
        hosts = [self.host] + ([self.fallback] if self.fallback else [])
        errors: list[str] = []
        members: dict[str, str] | None = None
        for host in hosts:
            try:
                members = self._run(host, unique, mode)
                break
            except RateLimited as error:
                errors.append(f"{host}: {error}")
        if members is None:
            raise RateLimited("; ".join(errors))

        taxonomy = parse_taxonomy(members.get(TAXONOMY_MEMBER, ""))
        merged: dict[int, list[str]] = {}
        for name in self._members_for(self.host, paired):
            if name not in members:
                continue
            text = members[name]
            if taxonomy:
                text = annotate_taxonomy(text, taxonomy)
            blocks = split_blocks(text)
            if len(blocks) == 1 and name.endswith(".a3m") and name[0].isdigit():
                # Protenix writes one numbered file per query rather than one
                # file holding every query's block.
                blocks = {FIRST_QUERY_ID + int(name.split(".")[0]): next(iter(blocks.values()))}
            for number, body in blocks.items():
                merged.setdefault(number, []).append(body)

        by_sequence: dict[str, MsaResult] = {}
        for index, sequence in enumerate(unique):
            blocks = merged.get(FIRST_QUERY_ID + index, [])
            if not blocks:
                raise SearchError(f"no alignment returned for query {index}")
            text = _rewrite_query_header(_join_blocks(blocks), sequence)
            result = MsaResult(sequence=sequence)
            if paired:
                result.paired = text
            else:
                result.unpaired = text
            by_sequence[sequence] = result
        return [by_sequence[sequence] for sequence in sequences]


def _join_blocks(blocks: list[str]) -> str:
    """Concatenate per-database blocks, keeping only the first copy of the query."""
    if len(blocks) == 1:
        return blocks[0]
    joined = [blocks[0]]
    for block in blocks[1:]:
        lines = block.splitlines(keepends=True)
        # Each database repeats the query as its first record.
        joined.append("".join(lines[2:]) if len(lines) > 2 else "")
    return "".join(joined)
