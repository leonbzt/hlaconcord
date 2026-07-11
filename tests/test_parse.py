import pytest

from hlaharm.nomenclature import ParseError, Resolution, parse_allele


def test_two_field_with_and_without_prefix():
    a = parse_allele("A*02:01")
    assert a.gene == "A"
    assert a.fields == ("02", "01")
    assert a.expression is None and a.group is None
    assert a.had_prefix is False
    assert a.resolution is Resolution.FIELD_2
    assert a.name() == "A*02:01"
    assert a.name(prefix=True) == "HLA-A*02:01"

    b = parse_allele("HLA-A*02:01")
    assert b.had_prefix is True
    assert b.name() == "A*02:01"


def test_full_resolution_with_expression_suffix():
    a = parse_allele("HLA-A*02:01:01:02L")
    assert a.fields == ("02", "01", "01", "02")
    assert a.expression == "L"
    assert a.group is None
    assert a.resolution is Resolution.FIELD_4
    assert a.name() == "A*02:01:01:02L"


def test_null_allele_is_flagged_and_preserved():
    a = parse_allele("A*02:11N")
    assert a.is_null is True
    assert a.expression == "N"
    assert a.name() == "A*02:11N"
    # reducing drops the suffix, but null-ness was readable on the original
    assert a.truncated(1).name() == "A*02"


@pytest.mark.parametrize(
    "text,fields,group",
    [
        ("A*02:01:01G", ("02", "01", "01"), "G"),
        ("HLA-B*07:02:01G", ("07", "02", "01"), "G"),
        ("A*02:01P", ("02", "01"), "P"),
    ],
)
def test_group_designations(text, fields, group):
    a = parse_allele(text)
    assert a.fields == fields
    assert a.group == group
    assert a.expression is None
    expected = Resolution.G_GROUP if group == "G" else Resolution.P_GROUP
    assert a.resolution is expected


def test_g_locus_is_not_a_group_suffix():
    # HLA-G is a locus; the leading G before '*' must parse as the gene.
    a = parse_allele("G*01:01:01:01")
    assert a.gene == "G"
    assert a.group is None
    assert a.fields == ("01", "01", "01", "01")


def test_single_field_low_resolution():
    a = parse_allele("A*02")
    assert a.fields == ("02",)
    assert a.resolution is Resolution.FIELD_1
    assert a.legacy is False


def test_legacy_colonless_form():
    a = parse_allele("A*0201")
    assert a.fields == ("02", "01")
    assert a.legacy is True
    assert a.name() == "A*02:01"


def test_legacy_with_expression_suffix():
    a = parse_allele("A*0211N")
    assert a.fields == ("02", "11")
    assert a.expression == "N"
    assert a.legacy is True


def test_truncation_bounds():
    a = parse_allele("A*02:01:01:01")
    assert a.truncated(2).name() == "A*02:01"
    with pytest.raises(ValueError):
        a.truncated(5)
    with pytest.raises(ValueError):
        a.truncated(0)


@pytest.mark.parametrize(
    "bad",
    [
        "",
        "   ",
        "not-an-allele",
        "A*",
        "HLA-",
        "A*02:01:01:01:01",  # too many fields
        "A*02:01X",  # unknown suffix
        "A*02010",  # odd-length legacy run needs the DB to disambiguate
    ],
)
def test_malformed_raises(bad):
    with pytest.raises(ParseError):
        parse_allele(bad)
