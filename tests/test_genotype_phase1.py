import numpy as np
import pytest

from statgen.genotype import load_genotype, load_genotype_cache, save_genotype_cache
from statgen.reference import load_reference
from tests.conftest import (
    FIXTURES_DIR,
    GENOTYPE_CHR1_BIM,
    GENOTYPE_CHR2_BIM,
    GENOTYPE_SAMPLES,
    write_plink_bed,
    write_plink_bim,
    write_plink_fam,
)


REF_SHARDED = FIXTURES_DIR / "reference/sharded/@.bim"
REF_NONSHARDED = FIXTURES_DIR / "reference/nonsharded/all.bim"
G_SHARDED = FIXTURES_DIR / "genotype/sharded/@"
G_NONSHARDED = FIXTURES_DIR / "genotype/nonsharded/all"


def test_sharded_genotype_metadata_loads():
    """@-sharded bfile loads correctly; chrX ploidy reflects X.ploidy sidecar (male=1, female=2)."""
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)  # X.ploidy exists — no warning expected

    assert panel.source_layout == "sharded"
    assert [s.label for s in panel.shards] == ["1", "X"]
    assert panel.num_snp == 8
    assert panel.num_sample == 4
    assert panel.is_present.tolist() == [True] * 8
    assert panel.source_row0.tolist() == [0, 1, 2, 3, 4, 0, 1, 2]
    # chr1: default diploid; chrX: hemizygous males (1), diploid females (2)
    assert panel.ploidy_male.tolist() == [2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0]
    assert panel.ploidy_female.tolist() == [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
    assert panel.fid.tolist() == ["FAM1", "FAM1", "FAM2", "FAM2"]
    assert panel.sex.tolist() == [1, 2, 1, 2]
    assert panel.is_male.tolist() == [True, False, True, False]
    assert panel.is_female.tolist() == [False, True, False, True]
    assert panel.is_subject_present("X").tolist() == [True] * 4


def test_nonsharded_genotype_loads_against_multi_shard_reference(tmp_path):
    """Non-sharded bfile produces multiple GenotypeShard objects all pointing to one .bed."""
    prefix = tmp_path / "all"
    prefix.with_suffix(".bim").write_text((FIXTURES_DIR / "reference/nonsharded/all.bim").read_text())
    prefix.with_suffix(".fam").write_text((FIXTURES_DIR / "genotype/sharded/1.fam").read_text())
    write_plink_bed(prefix.with_suffix(".bed"), num_snp=8, num_sample=4)

    ref = load_reference(REF_NONSHARDED)
    with pytest.warns(RuntimeWarning, match="chrX genotype source has no .ploidy"):
        panel = load_genotype(prefix, ref)

    assert panel.source_layout == "non_sharded"
    assert [s.bed_path for s in panel.shards] == [prefix.with_suffix(".bed")] * 2
    assert [s.source_num_snp for s in panel.shards] == [8, 8]
    assert panel.source_row0.tolist() == list(range(8))


def test_nonsharded_committed_fixture_loads():
    """Committed non-sharded bfile (with .ploidy) loads against a non-sharded reference."""
    ref = load_reference(REF_NONSHARDED)
    panel = load_genotype(G_NONSHARDED, ref)

    assert panel.source_layout == "non_sharded"
    assert len(panel.shards) == 2
    assert panel.num_snp == 8
    assert panel.is_present.tolist() == [True] * 8
    # chr1 rows (0-4): diploid (2, 2); chrX rows (5-7): male=1, female=2 from all.ploidy
    assert panel.ploidy_male.tolist() == [2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0]
    assert panel.ploidy_female.tolist() == [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]


def test_chrx_subset_fam_maps_to_panel_axis(tmp_path):
    """chrX FAM subset in non-autosomal order maps correctly to the panel sample axis."""
    for label in ("1", "X"):
        for suffix in (".bim", ".bed", ".fam"):
            src = FIXTURES_DIR / f"genotype/sharded/{label}{suffix}"
            dst = tmp_path / f"{label}{suffix}"
            dst.write_bytes(src.read_bytes())

    (tmp_path / "X.fam").write_text(
        "FAM2\tIND3\t0\t0\t1\t-9\n"
        "FAM1\tIND1\t0\t0\t1\t-9\n"
    )

    ref = load_reference(REF_SHARDED)
    panel = load_genotype(str(tmp_path / "@"), ref)

    assert panel.is_subject_present("X").tolist() == [True, False, True, False]
    assert panel.shards[1].source_subject_row0.tolist() == [1, -1, 0, -1]


def test_genotype_cache_roundtrip_and_subset(tmp_path):
    """Cache round-trip preserves all accessors; shard subset preserves panel sample axis."""
    ref = load_reference(REF_SHARDED)
    panel = load_genotype(G_SHARDED, ref)
    cache = tmp_path / "genotype.npz"
    save_genotype_cache(panel, cache)

    loaded = load_genotype_cache(cache, shards=["X"])
    assert loaded.source_layout == "sharded"
    assert [s.label for s in loaded.shards] == ["X"]
    assert loaded.num_snp == 3
    assert loaded.num_sample == panel.num_sample
    assert loaded.source_row0.tolist() == [0, 1, 2]
    assert np.array_equal(loaded.fid, panel.fid)
    assert loaded.ploidy_male.tolist() == [1.0, 1.0, 1.0]
    assert loaded.ploidy_female.tolist() == [2.0, 2.0, 2.0]


def test_missing_source_shard_fails(tmp_path):
    """@-sharded load fails clearly when a required source shard file is absent."""
    for suffix in (".bim", ".bed", ".fam"):
        (tmp_path / f"1{suffix}").write_bytes(
            (FIXTURES_DIR / f"genotype/sharded/1{suffix}").read_bytes()
        )
    # X shard files intentionally not copied

    ref = load_reference(REF_SHARDED)
    with pytest.raises(FileNotFoundError, match="Missing genotype source.*shard.*X"):
        load_genotype(str(tmp_path / "@"), ref)


def test_reference_subset_with_source_subset_loads(tmp_path):
    """Reference subset plus matching @-source subset loads without requiring other shards."""
    for suffix in (".bim", ".bed", ".fam"):
        (tmp_path / f"1{suffix}").write_bytes(
            (FIXTURES_DIR / f"genotype/sharded/1{suffix}").read_bytes()
        )

    ref = load_reference(REF_SHARDED, shards=["1"])
    panel = load_genotype(str(tmp_path / "@"), ref)

    assert panel.source_layout == "sharded"
    assert [s.label for s in panel.shards] == ["1"]
    assert panel.num_snp == 5


def test_extra_files_on_disk_ignored(tmp_path):
    """Extra bfile triplets on disk that are not in the reference are silently ignored."""
    for label in ("1", "X"):
        for suffix in (".bim", ".bed", ".fam"):
            src = FIXTURES_DIR / f"genotype/sharded/{label}{suffix}"
            (tmp_path / f"{label}{suffix}").write_bytes(src.read_bytes())

    # Extra chr2 bfile that the reference (chr1 + X) does not reference
    write_plink_bim(tmp_path / "2.bim", GENOTYPE_CHR2_BIM)
    write_plink_fam(tmp_path / "2.fam", GENOTYPE_SAMPLES)
    write_plink_bed(tmp_path / "2.bed", num_snp=2, num_sample=4)

    ref = load_reference(REF_SHARDED)  # only chr1 + X
    panel = load_genotype(str(tmp_path / "@"), ref)
    assert [s.label for s in panel.shards] == ["1", "X"]


def test_duplicate_fid_iid_fails(tmp_path):
    """FAM file with duplicate (fid, iid) pair fails at load time."""
    write_plink_bim(tmp_path / "dup.bim", GENOTYPE_CHR1_BIM)
    (tmp_path / "dup.fam").write_text(
        "FAM1\tIND1\t0\t0\t1\t-9\n"
        "FAM1\tIND1\t0\t0\t2\t-9\n"  # same (fid, iid) as row 1
    )
    write_plink_bed(tmp_path / "dup.bed", num_snp=5, num_sample=2)

    ref = load_reference(REF_SHARDED, shards=["1"])
    with pytest.raises(ValueError, match="duplicate FAM subject"):
        load_genotype(str(tmp_path / "dup"), ref)


def test_invalid_sex_fails(tmp_path):
    """FAM file with sex value outside {0, 1, 2} fails at load time."""
    write_plink_bim(tmp_path / "badsex.bim", GENOTYPE_CHR1_BIM)
    (tmp_path / "badsex.fam").write_text("FAM1\tIND1\t0\t0\t3\t-9\n")  # sex=3 is invalid
    write_plink_bed(tmp_path / "badsex.bed", num_snp=5, num_sample=1)

    ref = load_reference(REF_SHARDED, shards=["1"])
    with pytest.raises(ValueError, match="FAM sex must be one of"):
        load_genotype(str(tmp_path / "badsex"), ref)


def test_autosomal_fam_mismatch_fails_on_exposed_columns(tmp_path):
    """Sharded autosomal shards with differing sex in FAM fail at load time."""
    # Build a two-autosome reference (chr1 + chr2)
    write_plink_bim(tmp_path / "ref1.bim", GENOTYPE_CHR1_BIM)
    write_plink_bim(tmp_path / "ref2.bim", GENOTYPE_CHR2_BIM)
    ref = load_reference(str(tmp_path / "ref@.bim"))

    # chr1 genotype with standard FAM
    write_plink_bim(tmp_path / "1.bim", GENOTYPE_CHR1_BIM)
    write_plink_fam(tmp_path / "1.fam", GENOTYPE_SAMPLES)
    write_plink_bed(tmp_path / "1.bed", num_snp=5, num_sample=4)

    # chr2 genotype with first sample sex changed (1 → 2)
    samples_diff_sex = list(GENOTYPE_SAMPLES)
    samples_diff_sex[0] = ("FAM1", "IND1", 0, 0, 2, -9)
    write_plink_bim(tmp_path / "2.bim", GENOTYPE_CHR2_BIM)
    write_plink_fam(tmp_path / "2.fam", samples_diff_sex)
    write_plink_bed(tmp_path / "2.bed", num_snp=2, num_sample=4)

    with pytest.raises(ValueError, match="autosomal FAM mismatch"):
        load_genotype(str(tmp_path / "@"), ref)


def test_autosomal_fam_phenotype_mismatch_does_not_fail(tmp_path):
    """Sharded autosomal FAM phenotype differences are ignored; load succeeds."""
    write_plink_bim(tmp_path / "ref1.bim", GENOTYPE_CHR1_BIM)
    write_plink_bim(tmp_path / "ref2.bim", GENOTYPE_CHR2_BIM)
    ref = load_reference(str(tmp_path / "ref@.bim"))

    write_plink_bim(tmp_path / "1.bim", GENOTYPE_CHR1_BIM)
    write_plink_fam(tmp_path / "1.fam", GENOTYPE_SAMPLES)  # pheno = -9
    write_plink_bed(tmp_path / "1.bed", num_snp=5, num_sample=4)

    samples_diff_pheno = [(*s[:5], 0) for s in GENOTYPE_SAMPLES]  # pheno = 0 (different but not exposed)
    write_plink_bim(tmp_path / "2.bim", GENOTYPE_CHR2_BIM)
    write_plink_fam(tmp_path / "2.fam", samples_diff_pheno)
    write_plink_bed(tmp_path / "2.bed", num_snp=2, num_sample=4)

    panel = load_genotype(str(tmp_path / "@"), ref)
    assert panel.num_sample == 4


def test_chrx_absent_subject_fails(tmp_path):
    """chrX FAM containing a subject absent from the autosomal axis fails."""
    for label in ("1", "X"):
        for suffix in (".bim", ".bed", ".fam"):
            src = FIXTURES_DIR / f"genotype/sharded/{label}{suffix}"
            (tmp_path / f"{label}{suffix}").write_bytes(src.read_bytes())

    (tmp_path / "X.fam").write_text(
        "FAM1\tIND1\t0\t0\t1\t-9\n"
        "FAM9\tINDX\t0\t0\t1\t-9\n"  # subject absent from chr1 FAM
    )

    ref = load_reference(REF_SHARDED)
    with pytest.raises(ValueError, match="chrX FAM subject.*absent from the autosomal sample axis"):
        load_genotype(str(tmp_path / "@"), ref)


def test_chrx_metadata_mismatch_fails(tmp_path):
    """chrX FAM sex differing from the autosomal value for a matched subject fails."""
    for label in ("1", "X"):
        for suffix in (".bim", ".bed", ".fam"):
            src = FIXTURES_DIR / f"genotype/sharded/{label}{suffix}"
            (tmp_path / f"{label}{suffix}").write_bytes(src.read_bytes())

    # IND1 is male (sex=1) in chr1 FAM; claim female (sex=2) in chrX FAM
    (tmp_path / "X.fam").write_text("FAM1\tIND1\t0\t0\t2\t-9\n")

    ref = load_reference(REF_SHARDED)
    with pytest.raises(ValueError, match="chrX FAM metadata mismatch"):
        load_genotype(str(tmp_path / "@"), ref)


def test_chrx_only_load_defines_sample_axis():
    """When only chrX is loaded, the chrX FAM defines the panel-level sample axis."""
    ref = load_reference(REF_SHARDED).select_shards(["X"])
    panel = load_genotype(G_SHARDED, ref)

    assert panel.source_layout == "sharded"
    assert [s.label for s in panel.shards] == ["X"]
    assert panel.num_sample == 4
    assert panel.is_subject_present("X").tolist() == [True] * 4


def test_chrx_without_ploidy_warns_and_defaults_to_diploid(tmp_path):
    """chrX load without .ploidy emits the documented warning; matched rows default to (2, 2)."""
    for label in ("1", "X"):
        for suffix in (".bim", ".bed", ".fam"):
            src = FIXTURES_DIR / f"genotype/sharded/{label}{suffix}"
            (tmp_path / f"{label}{suffix}").write_bytes(src.read_bytes())
    # X.ploidy intentionally not copied

    ref = load_reference(REF_SHARDED)
    with pytest.warns(RuntimeWarning, match="chrX genotype source has no .ploidy"):
        panel = load_genotype(str(tmp_path / "@"), ref)

    chrx = panel.shards[1]
    assert chrx.ploidy_male.tolist() == [2.0, 2.0, 2.0]
    assert chrx.ploidy_female.tolist() == [2.0, 2.0, 2.0]


def test_is_present_false_for_absent_variants(tmp_path):
    """Reference SNPs not present in the genotype source have is_present=False and NaN ploidy."""
    # chr1-only reference with one extra SNP (rs1006) absent from the genotype fixture
    extra_bim = (
        "1\trs1001\t0\t100\tA\tG\n"
        "1\trs1002\t0\t200\tC\tT\n"
        "1\trs1003\t0\t300\tA\tC\n"
        "1\trs1004\t0\t400\tG\tA\n"
        "1\trs1005\t0\t500\tT\tC\n"
        "1\trs1006\t0\t600\tA\tC\n"
    )
    (tmp_path / "extra.bim").write_text(extra_bim)
    ref = load_reference(tmp_path / "extra.bim")  # non-sharded, chr1 only, 6 SNPs

    panel = load_genotype(G_SHARDED, ref)  # source has only 5 chr1 SNPs

    assert panel.num_snp == 6
    assert panel.is_present.tolist() == [True, True, True, True, True, False]
    assert panel.source_row0.tolist() == [0, 1, 2, 3, 4, -1]
    assert panel.ploidy_male[:5].tolist() == [2.0, 2.0, 2.0, 2.0, 2.0]
    assert np.isnan(panel.ploidy_male[5])
    assert panel.ploidy_female[:5].tolist() == [2.0, 2.0, 2.0, 2.0, 2.0]
    assert np.isnan(panel.ploidy_female[5])
