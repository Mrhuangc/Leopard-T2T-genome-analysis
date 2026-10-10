from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]


def _scheduler_rows() -> list[tuple[str, int, int, str, int, int, int]]:
    script = (ROOT / "bin" / "run_species_pipeline.sh").read_text(encoding="utf-8")
    matches = re.findall(
        r'^\s+"([a-z][a-z0-9_]+)\|(\d+)\|(\d+)\|'
        r'(bundle|distributed)\|(\d+)\|(\d+)\|(\d+)"$',
        script,
        flags=re.MULTILINE,
    )
    return [
        (name, int(first), int(last), mode, int(cpus), int(memory), int(hours))
        for name, first, last, mode, cpus, memory, hours in matches
    ]


def test_optimized_scheduler_covers_every_stage_once() -> None:
    rows = _scheduler_rows()
    assert rows
    covered = [stage for _, first, last, *_ in rows for stage in range(first, last + 1)]
    assert covered == list(range(1, 20))


def test_scheduler_does_not_shadow_bash_groups_and_has_final_gate() -> None:
    script = (ROOT / "bin" / "run_species_pipeline.sh").read_text(encoding="utf-8")
    assert "TASK_GROUPS=(" in script
    assert "\nGROUPS=(" not in script
    assert "FINAL STAGE-GATE AUDIT" in script
    assert '((SELECTED_STAGES == EXPECTED_STAGES))' in script


def test_only_exonerate_scatter_stages_are_distributed() -> None:
    rows = _scheduler_rows()
    distributed = {(first, last) for _, first, last, mode, *_ in rows if mode == "distributed"}
    assert distributed == {(10, 10), (12, 12)}


def test_bundles_have_bounded_positive_resources() -> None:
    for _, _, _, mode, cpus, memory, hours in _scheduler_rows():
        if mode == "bundle":
            assert 1 <= cpus <= 16
            assert 4 <= memory <= 64
            assert 1 <= hours <= 36
        else:
            assert (cpus, memory, hours) == (0, 0, 0)


def test_batch_profile_is_strictly_serial_and_hardlink_is_default() -> None:
    config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    assert "publish_mode  = 'link'" in config
    batch = config.split("    batch {", maxsplit=1)[1].split("    slurm {", maxsplit=1)[0]
    assert "executor = 'local'" in batch
    assert "maxForks = 1" in batch
    assert "executor.queueSize = 1" in batch
