import pytest

from statgen.annotations import AnnotationPanel
from statgen.genotype import GenotypePanel
from statgen.ld import LDPanel
from statgen.reference import ReferencePanel
from statgen.sumstats import Sumstats
from tests.conftest import FIXTURES_DIR, matlab_data_lines, run_octave, skipif_no_octave


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


@pytest.mark.octave
@skipif_no_octave
def test_octave_panel_and_shard_display_methods_are_safe():
    fixture_dir = str(FIXTURES_DIR)
    script = (
        f"fixture_dir = '{fixture_dir}'; "
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "sumstats = statgen.load_sumstats([fixture_dir '/sumstats/traits.tsv.gz'], ref); "
        "annotations = statgen.load_annotations({[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}, ref); "
        "genotype = statgen.load_genotype([fixture_dir '/genotype/sharded/@'], ref); "
        "ld = statgen.load_ld([fixture_dir '/ld/matlab']); "
        "objects = {ref, ref.shards{1}, sumstats, sumstats.shards{1}, annotations, annotations.shards{1}, genotype, genotype.shards{1}, ld, ld.shards{1}}; "
        "needles = {'ReferencePanel', 'ReferenceShard', 'Sumstats', 'SumstatsShard', 'AnnotationPanel', 'AnnotationShard', 'GenotypePanel', 'GenotypeShard', 'LDPanel', 'LDShard'}; "
        "ok = zeros(1, numel(objects)); "
        "for i = 1:numel(objects); txt = evalc('disp(objects{i})'); ok(i) = ~isempty(strfind(txt, needles{i})); end; "
        "fprintf('%d ', ok); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert matlab_data_lines(result.stdout) == ["1 1 1 1 1 1 1 1 1 1"]
