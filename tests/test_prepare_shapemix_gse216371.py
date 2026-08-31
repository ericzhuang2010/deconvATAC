import xml.etree.ElementTree as ET

import pytest

from scripts.prepare_shapemix_gse216371 import (
    cell_text,
    excel_column_index,
    normalize_stage,
    parse_ccre_name,
)


@pytest.mark.parametrize(
    ("reference", "expected"),
    [("A1", 0), ("Z9", 25), ("AA10", 26), ("AB10", 27)],
)
def test_excel_column_index(reference, expected):
    assert excel_column_index(reference) == expected


@pytest.mark.parametrize("value", ["E13.5", "e13.5", "13.5", "13.50"])
def test_normalize_stage_accepts_author_e13_5_forms(value):
    assert normalize_stage(value) == "E13.5"


def test_cell_text_resolves_shared_and_inline_strings():
    shared = ET.fromstring('<c r="A1" t="s"><v>1</v></c>')
    inline = ET.fromstring('<c r="B1" t="inlineStr"><is><t>label</t></is></c>')
    assert cell_text(shared, ["zero", "one"]) == "one"
    assert cell_text(inline, []) == "label"


def test_excel_column_index_rejects_missing_letters():
    with pytest.raises(ValueError, match="no column letters"):
        excel_column_index("123")


def test_parse_ccre_name_preserves_contig_and_builds_canonical_id():
    assert parse_ccre_name("chrUn_GL456239_100_600") == (
        "chrUn_GL456239",
        100,
        600,
        "chrUn_GL456239:100-600",
    )


@pytest.mark.parametrize("value", ["chr1_bad_500", "1_100_600", "chr1_600_100"])
def test_parse_ccre_name_rejects_invalid_values(value):
    with pytest.raises(ValueError, match="Invalid author cCRE"):
        parse_ccre_name(value)
