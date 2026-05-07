import pytest

from statgen.annotations import AnnotationPanel
from statgen.genotype import GenotypePanel
from statgen.ld import LDPanel
from statgen.reference import ReferencePanel
from statgen.sumstats import Sumstats
from tests.conftest import matlab_data_lines, run_octave, skipif_no_octave


def test_python_panel_constructors_reject_zero_shards():
    cases = [
        (ReferencePanel, ([],)),
        (LDPanel, ([],)),
        (AnnotationPanel, ([], ["baseline"])),
        (Sumstats, ([],)),
    ]
    for constructor, args in cases:
        with pytest.raises(ValueError, match="requires at least one"):
            constructor(*args)

    with pytest.raises(ValueError, match="requires at least one"):
        GenotypePanel(
            [],
            fid=[],
            iid=[],
            father_id=[],
            mother_id=[],
            sex=[],
            source_layout="sharded",
        )


@pytest.mark.octave
@skipif_no_octave
def test_octave_panel_constructors_reject_zero_shards():
    script = (
        "ok = zeros(1, 5); "
        "try; statgen.ReferencePanel({}); catch; ok(1) = 1; end; "
        "try; statgen.LDPanel({}); catch; ok(2) = 1; end; "
        "try; statgen.AnnotationPanel({}, {'baseline'}); catch; ok(3) = 1; end; "
        "try; statgen.Sumstats({}); catch; ok(4) = 1; end; "
        "try; statgen.GenotypePanel({}, {}, {}, {}, {}, [], 'sharded'); catch; ok(5) = 1; end; "
        "fprintf('%d ', ok); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 1 1 1 1"]
