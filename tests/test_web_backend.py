"""The web client is tested against a stand-in service, never the real one."""

import io
import tarfile

import pytest

from foldrunner.msa.backends import RateLimited, SearchError, WebBackend
from foldrunner.msa.backends.web import COLABFOLD_HOST, PROTENIX_HOST, split_blocks

SEQ_A = "MKVLAAGIVG"
SEQ_B = "PQRSTVWYAC"

# Two query blocks in one reply, separated by a null byte and numbered from 101.
UNIREF_REPLY = (
    ">101\n"
    f"{SEQ_A}\n"
    ">UniRef100_A_9606/\t1\t2\t3\n"
    f"{SEQ_A}\n\x00"
    ">102\n"
    f"{SEQ_B}\n"
    ">UniRef100_B_10090/\t1\t2\t3\n"
    f"{SEQ_B}\n"
)


def make_archive(members: dict[str, str]) -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        for name, text in members.items():
            data = text.encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


class Reply:
    def __init__(self, payload=None, content=b""):
        self._payload = payload
        self.content = content

    def json(self):
        if self._payload is None:
            raise ValueError("not json")
        return self._payload


class FakeService:
    """Answers the ticket protocol with a canned archive."""

    def __init__(self, members=None, submit_status="COMPLETE", download=None):
        self.members = members or {"uniref.a3m": UNIREF_REPLY}
        self.submit_status = submit_status
        self.download = download
        self.posts = []
        self.headers = []

    def post(self, url, data, timeout, headers):
        self.posts.append((url, data))
        self.headers.append(headers)
        return Reply({"id": "t1", "status": self.submit_status})

    def get(self, url, timeout, headers):
        if "/result/download/" in url:
            if self.download is not None:
                return Reply(content=self.download)
            return Reply(content=make_archive(self.members))
        return Reply({"id": "t1", "status": "COMPLETE"})


@pytest.fixture
def backend():
    return WebBackend(contact="someone@example.org", poll_seconds=0, session=FakeService())


def test_contact_is_required():
    with pytest.raises(SearchError, match="identify themselves"):
        WebBackend()


def test_blocks_are_split_on_the_null_separator():
    blocks = split_blocks(UNIREF_REPLY)
    assert sorted(blocks) == [101, 102]
    assert SEQ_A in blocks[101] and SEQ_B in blocks[102]


def test_each_query_gets_its_own_alignment(backend):
    results = backend.fetch([SEQ_A, SEQ_B])
    assert [r.sequence for r in results] == [SEQ_A, SEQ_B]
    assert SEQ_A in results[0].unpaired and SEQ_B not in results[0].unpaired


def test_repeated_sequences_are_submitted_once(backend):
    backend.fetch([SEQ_A, SEQ_B, SEQ_A])
    submitted = backend.session.posts[0][1]["q"]
    assert submitted.count(SEQ_A) == 1


def test_results_follow_the_order_asked_for(backend):
    results = backend.fetch([SEQ_B, SEQ_A, SEQ_B])
    assert [r.sequence for r in results] == [SEQ_B, SEQ_A, SEQ_B]


def test_query_row_is_labelled_for_downstream_parsers(backend):
    from foldrunner.msa.a3m import parse_a3m

    result = backend.fetch([SEQ_A])[0]
    assert parse_a3m(result.unpaired).query == SEQ_A


def test_pairing_returns_the_pairable_half_of_one_search():
    """Verified against both live services: the pairing endpoint pairs across
    whatever was submitted together, so it answers a single sequence with
    nothing and answers a complex with rows carrying no taxonomy. Deriving the
    pairing half from a taxonomy-carrying search keeps one search per sequence
    instead of one per pair."""
    service = FakeService(members={"uniref.a3m": UNIREF_REPLY})
    backend = WebBackend(contact="a@b.c", poll_seconds=0, session=service)
    results = backend.fetch([SEQ_A, SEQ_B], paired=True)
    assert all("ticket/msa" in url for url, _ in service.posts)
    assert results[0].paired is not None and results[0].unpaired is None


def test_unpaired_uses_the_msa_endpoint(backend):
    backend.fetch([SEQ_A])
    assert "ticket/msa" in backend.session.posts[0][0]
    assert backend.session.posts[0][1]["mode"] == "env"


def test_client_identifies_itself(backend):
    backend.fetch([SEQ_A])
    assert backend.session.headers[0]["User-Agent"].startswith("foldrunner/")
    assert "someone@example.org" in backend.session.headers[0]["User-Agent"]


def test_throttled_html_body_is_reported_as_throttling_not_corruption():
    service = FakeService(download=b"<html>Too many requests</html>")
    backend = WebBackend(
        contact="a@b.c", poll_seconds=0, fallback=None, max_attempts=1, session=service
    )
    with pytest.raises(RateLimited, match="not an archive"):
        backend.fetch([SEQ_A])


def test_maintenance_is_not_retried():
    service = FakeService(submit_status="MAINTENANCE")
    backend = WebBackend(contact="a@b.c", poll_seconds=0, session=service)
    with pytest.raises(SearchError, match="maintenance"):
        backend.fetch([SEQ_A])


def test_the_second_host_takes_over_when_the_first_throttles():
    class Failover(FakeService):
        def get(self, url, timeout, headers):
            if "/result/download/" in url and PROTENIX_HOST in url:
                return Reply(content=b"<html>rate limited</html>")
            return super().get(url, timeout, headers)

    service = Failover()
    backend = WebBackend(
        host=PROTENIX_HOST,
        fallback=COLABFOLD_HOST,
        contact="a@b.c",
        poll_seconds=0,
        max_attempts=1,
        session=service,
    )
    assert backend.fetch([SEQ_A])[0].unpaired is not None
    assert any(COLABFOLD_HOST in url for url, _ in service.posts)


def test_recipe_records_which_service_and_databases_were_used():
    protenix = WebBackend(contact="a@b.c", host=PROTENIX_HOST).recipe
    colabfold = WebBackend(contact="a@b.c", host=COLABFOLD_HOST).recipe
    assert protenix.backend == "protenix-api"
    assert colabfold.backend == "colabfold-api"
    assert protenix.recipe_id != colabfold.recipe_id
    assert "colabfold_envdb" in protenix.databases


def test_turning_off_the_environmental_database_changes_the_recipe():
    with_env = WebBackend(contact="a@b.c").recipe
    without = WebBackend(contact="a@b.c", use_env=False).recipe
    assert with_env.recipe_id != without.recipe_id


# --- what the live services actually return ---------------------------------

TAXONOMY = (
    "101\tA0A6C0JTL4\t1070528\tSome species\t-_cellular organisms\n"
    "101\tE9L6B9\t10455\tOther\t-_x\n"
)
PROTENIX_A3M = (
    ">101\n"
    f"{SEQ_A}\n"
    ">A0A6C0JTL4\t107\t0.894\t1.606E-24\n"
    f"{SEQ_A}\n"
    ">MGYP001411914921\t111\t0.855\t6.097E-26\n"
    f"{SEQ_A}\n"
)


def test_a_plain_tar_is_accepted_as_well_as_a_gzip_stream():
    """The services send Content-Encoding: gzip, so an HTTP client that honours
    it hands back an already-decompressed tar."""
    from foldrunner.msa.cache import looks_like_archive

    plain = make_archive({"uniref.a3m": UNIREF_REPLY})
    import gzip as gziplib

    assert looks_like_archive(plain)                       # gzip stream
    assert looks_like_archive(gziplib.decompress(plain))   # bare tar
    assert not looks_like_archive(b"<html>rate limited</html>")


def test_protenix_style_numbered_members_are_read():
    service = FakeService(members={"0.a3m": PROTENIX_A3M, "uniref_tax.m8": TAXONOMY})
    backend = WebBackend(contact="a@b.c", poll_seconds=0, host=PROTENIX_HOST, session=service)
    result = backend.fetch([SEQ_A])[0]
    assert result.unpaired is not None
    assert SEQ_A in result.unpaired


def test_taxonomy_is_folded_into_the_headers():
    """Hits come back as bare accessions; pairing needs the organism on them."""
    service = FakeService(members={"0.a3m": PROTENIX_A3M, "uniref_tax.m8": TAXONOMY})
    backend = WebBackend(contact="a@b.c", poll_seconds=0, host=PROTENIX_HOST, session=service)
    from foldrunner.msa.a3m import parse_a3m

    alignment = parse_a3m(backend.fetch([SEQ_A])[0].unpaired)
    assert "1070528" in alignment.species
    # The metagenomic hit has no taxonomy and keeps none.
    assert None in alignment.species[1:]


def test_the_pairing_half_is_derived_not_searched_again():
    """The pairing endpoint answers a single sequence with nothing, and answers a
    complex with rows usable only for that complex."""
    service = FakeService(members={"0.a3m": PROTENIX_A3M, "uniref_tax.m8": TAXONOMY})
    backend = WebBackend(contact="a@b.c", poll_seconds=0, host=PROTENIX_HOST, session=service)
    result = backend.fetch([SEQ_A], paired=True)[0]
    assert result.paired is not None
    assert all("ticket/pair" not in url for url, _ in service.posts)


def test_the_pairing_half_keeps_only_hits_with_a_known_organism():
    from foldrunner.msa.a3m import parse_a3m
    from foldrunner.msa.backends.web import pairable

    kept = parse_a3m(pairable(PROTENIX_A3M.replace(">A0A6C0JTL4", ">UniRef100_A0A6C0JTL4_9606/")))
    assert kept.depth == 2                    # query plus the one with a taxon
    assert kept.species[1] == "9606"


def test_only_one_service_offers_pairing():
    """They look interchangeable but the Protenix host answers ticket/pair with 404."""
    assert WebBackend(contact="a@b.c", host=COLABFOLD_HOST).supports_pairing()
    assert not WebBackend(contact="a@b.c", host=PROTENIX_HOST).supports_pairing()
