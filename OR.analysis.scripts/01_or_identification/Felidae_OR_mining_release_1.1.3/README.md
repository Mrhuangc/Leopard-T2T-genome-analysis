# Felidae OR Mining 1.1.3

Versioned OR-identification workflow with a Python package, Nextflow stages, Conda environments and a resumable launcher.

Input parameters and commands for this repository are in the [OR analysis guide](../../README.md). Detailed identification usage is in [USAGE_CN.md](docs/USAGE_CN.md), and output directories are described in [OUTPUTS_CN.md](docs/OUTPUTS_CN.md).

```bash
bash bootstrap_control_env.sh
conda activate felidae-or-control
cp params.template.yml params.Ppa.yml
# Set external input paths and output/cluster settings.
bash bin/run_species_pipeline.sh --params params.Ppa.yml --preflight-only
bash bin/run_species_pipeline.sh --params params.Ppa.yml --profile slurm,conda
```

The `tests/data/initial_search/` FASTAs are small synthetic software-test fixtures.

