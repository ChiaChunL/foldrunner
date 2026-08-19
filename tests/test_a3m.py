import pytest

from foldrunner.msa.a3m import (
    pair_alignments,
    parse_a3m,
    render_block_diagonal,
    render_paired,
    species_id,
    ungapped_length,
)

# Header shapes taken from real service output: the pairing endpoint returns
# UniRef hits carrying the taxon, the unpaired endpoint returns metagenomic hits
# that carry none.
PAIRED_A = (
    ">query\nMKVL\n"
    ">UniRef100_UPI0008F9D4CD_9606/\t102\t0.351\t2.228E-22\nMKVI\n"
    ">UniRef100_K0A1Z2_10090/\t100\t0.355\t7.911E-22\nMKVV\n"
    ">UniRef100_A0A7M5VAS7_7227/\t99\t0.303\t1.491E-21\nMKAA\n"
)
PAIRED_B = (
    ">query\nPQRS\n"
    ">UniRef100_A0A5K1UL34_9606/\t325\t0.985\t6.248E-94\nPQRT\n"
    ">UniRef100_UPI0013F3A27B_10090/\t302\t0.659\t8.586E-86\nPQRR\n"
    ">UniRef100_ZZZ_9823/\t120\t0.5\t1E-30\nPQAA\n"
)
UNPAIRED_A = ">query\nMKVL\n>ERR1719370_563322\t87\t0.311\t3.413E-17\nMKCC\n"


def test_taxon_is_read_from_uniref_headers():
    assert species_id(">UniRef100_UPI0008F9D4CD_7038/\t102") == "7038"
    assert species_id("UniRef100_A0A7M5VAS7_252671/") == "252671"


def test_uniprot_style_headers_are_also_understood():
    assert species_id("tr|A0A123|SOMENAME_HUMAN extra") == "HUMAN"


def test_metagenomic_headers_carry_no_taxon():
    assert species_id(">ERR1719370_563322\t87\t0.311") is None


def test_parse_strips_the_null_separator_used_between_service_blocks():
    alignment = parse_a3m(">query\x00\nMKVL\n>UniRef100_A_9606/\nMKVI\n")
    assert alignment.depth == 2
    assert alignment.query == "MKVL"


def test_insertions_do_not_count_toward_alignment_width():
    assert ungapped_length("MKvvVL") == 4


def test_pairing_keeps_shared_taxa_and_drops_single_chain_ones():
    a, b = parse_a3m(PAIRED_A), parse_a3m(PAIRED_B)
    pairing = pair_alignments([a, b])
    # Row 0 is always the pair of queries.
    assert pairing[0] == [0, 0]
    taxa = {a.species[entry[0]] for entry in pairing[1:] if entry[0] is not None}
    assert taxa == {"9606", "10090"}  # 7227 only in A, 9823 only in B


def test_require_all_chains_is_stricter_than_greedy():
    a = parse_a3m(PAIRED_A)
    b = parse_a3m(">query\nPQRS\n>UniRef100_X_9606/\nPQRT\n")
    c = parse_a3m(">query\nWYFG\n>UniRef100_Y_10090/\nWYFH\n")
    greedy = pair_alignments([a, b, c])
    strict = pair_alignments([a, b, c], require_all_chains=True)
    assert len(greedy) > len(strict)
    assert len(strict) == 1  # no taxon appears in all three


def test_absent_chains_are_gapped_not_dropped():
    a = parse_a3m(PAIRED_A)
    b = parse_a3m(">query\nPQRS\n>UniRef100_X_9606/\nPQRT\n>UniRef100_Y_10090/\nPQRR\n")
    c = parse_a3m(">query\nWYFG\n>UniRef100_Z_9606/\nWYFH\n")
    rows = render_paired([a, b, c], pair_alignments([a, b, c]))
    assert rows[0] == "MKVLPQRSWYFG"
    # 10090 is missing from chain C, which is gapped rather than excluded.
    assert any(row.endswith("----") for row in rows[1:])


def test_unpaired_block_is_diagonal():
    a, b = parse_a3m(UNPAIRED_A), parse_a3m(">query\nPQRS\n>ERR9_1\nPQAA\n")
    rows = render_block_diagonal([a, b])
    assert rows == ["MKCC----", "----PQAA"]


def test_pairing_needs_at_least_two_chains():
    with pytest.raises(ValueError, match="at least two chains"):
        pair_alignments([parse_a3m(PAIRED_A)])


def test_max_rows_caps_the_paired_block():
    a, b = parse_a3m(PAIRED_A), parse_a3m(PAIRED_B)
    assert len(pair_alignments([a, b], max_rows=2)) == 2
