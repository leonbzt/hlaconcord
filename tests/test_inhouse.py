from hlaharm.nomenclature import ReductionBasis, ValidationStatus

# -- reference loading --------------------------------------------------------

def test_reference_loads(reference):
    assert reference.version == "3.55.0"
    assert "A*02:01:01:01" in reference.alleles
    assert reference.g_group_of["A*02:01:01:02L"] == "A*02:01:01G"
    assert "A*02:01:01G" in reference.g_groups
    assert reference.accession_of("A*02:01:01:01") == "HLA00005"


# -- validation ---------------------------------------------------------------

def test_validate_exact(nom):
    assert nom.validate(nom.parse("A*02:01:01:01")) is ValidationStatus.EXACT
    assert nom.validate(nom.parse("A*02:11N")) is ValidationStatus.EXACT


def test_validate_valid_reduction(nom):
    # A real 2-/3-field prefix of an existing allele, not a full allele itself.
    assert nom.validate(nom.parse("A*02:01")) is ValidationStatus.VALID_REDUCTION
    assert nom.validate(nom.parse("A*02:01:01")) is ValidationStatus.VALID_REDUCTION


def test_validate_group_designations(nom):
    assert nom.validate(nom.parse("A*02:01:01G")) is ValidationStatus.G_GROUP
    assert nom.validate(nom.parse("A*02:01P")) is ValidationStatus.P_GROUP
    assert nom.validate(nom.parse("A*99:99:99G")) is ValidationStatus.UNKNOWN


def test_validate_renamed_via_history(nom):
    # Present in allele history but not valid in the target release -> flagged, not coerced.
    assert nom.validate(nom.parse("A*66:01N")) is ValidationStatus.RENAMED


def test_validate_unknown(nom):
    assert nom.validate(nom.parse("A*99:99")) is ValidationStatus.UNKNOWN


# -- reduction ----------------------------------------------------------------

def test_reduce_two_field(nom):
    assert nom.reduce(nom.parse("A*02:01:01:01"), ReductionBasis.TWO_FIELD) == "A*02:01"
    assert nom.reduce(nom.parse("A*02:01"), ReductionBasis.TWO_FIELD) == "A*02:01"


def test_reduce_g_group_exact_hit(nom):
    assert nom.reduce(nom.parse("A*02:01:01:01"), ReductionBasis.G) == "A*02:01:01G"
    # a null/low-expression member still maps to its G-group
    assert nom.reduce(nom.parse("A*02:01:01:02L"), ReductionBasis.G) == "A*02:01:01G"


def test_reduce_g_group_identity(nom):
    assert nom.reduce(nom.parse("A*02:01:01G"), ReductionBasis.G) == "A*02:01:01G"


def test_reduce_g_group_via_unambiguous_prefix(nom):
    # A 2-field call maps to a G-group only when all fuller alleles agree on one.
    assert nom.reduce(nom.parse("B*07:02"), ReductionBasis.G) == "B*07:02:01G"


def test_reduce_g_group_singleton_returns_allele_unchanged(nom):
    # C*07:02:01:01 has no multi-member G-group -> it is its own ARD representative.
    # Returning it unchanged (not 2-field) avoids over-merging distinct ARD sequences.
    assert nom.reduce(nom.parse("C*07:02:01:01"), ReductionBasis.G) == "C*07:02:01:01"


def test_reduce_never_drops_null_suffix(nom):
    # A*02:11N has no G-group; every basis must keep it null, never collapse to A*02:11.
    null = nom.parse("A*02:11N")
    assert nom.reduce(null, ReductionBasis.G) == "A*02:11N"
    assert nom.reduce(null, ReductionBasis.LGX) == "A*02:11N"
    assert nom.reduce(null, ReductionBasis.TWO_FIELD) == "A*02:11N"


def test_reduce_lgx(nom):
    assert nom.reduce(nom.parse("A*02:01:01:01"), ReductionBasis.LGX) == "A*02:01"
    assert nom.reduce(nom.parse("A*02:01:01G"), ReductionBasis.LGX) == "A*02:01"


def test_reduce_p_group(nom):
    assert nom.reduce(nom.parse("A*02:01:01:01"), ReductionBasis.P) == "A*02:01P"


# -- accession ----------------------------------------------------------------

def test_accession_lookup(nom):
    assert nom.accession(nom.parse("A*02:01:01:01")) == "HLA00005"
    assert nom.accession(nom.parse("A*99:99")) is None
