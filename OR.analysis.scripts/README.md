# Felid olfactory-receptor analyses

Scripts for OR identification, subfamily assignment, copy-number comparison, regional synteny and manuscript figures. The module contains source code, environment files, input configuration and usage instructions. Genome FASTAs, HiFi reads, protein databases and study input/result snapshots are supplied separately.

[中文说明](README_CN.md)

## Assemblies and layout

| Analysis unit | Species | Assembly label |
| --- | --- | --- |
| Fca | *Felis catus* | Fca |
| Lge | *Leopardus geoffroyi* | Lge |
| Pbe | *Prionailurus bengalensis* | Pbe |
| Ple | *Panthera leo* | Ple |
| Ppa | *Panthera pardus* | Ppa_T2T / PPA |
| Ppa_or | *Panthera pardus orientalis* | PpoHapG |
| Pti | Sumatran tiger | Pti |

Keep these case-sensitive analysis-unit codes in OR tables and FASTA identifiers. Ppa and Ppa_or are separate assembly units. PPA inputs use the final assembly coordinates, including the finalized ChrY orientation.

| Directory | Contents |
| --- | --- |
| `01_or_identification/` | OR-mining workflow 1.1.3, Python package, Nextflow modules and tests |
| `02_or_subfamily_assignment/` | DIAMOND candidates, family trees and tree-assisted assignments |
| `03_main_figure/` | Figure 5A–D and OR6C regional-alignment scripts |
| `04_supplementary_figures/` | Figures S1–S6 and chromosome-region alignment scripts |
| `05_supporting_analysis/` | Functional OR copy-number matrices and heatmap |
| `config/` | Seven assembly parameter templates, species metadata and table schemas |

[`scripts_manifest.tsv`](scripts_manifest.tsv) lists module inputs and outputs. [`config/input_schemas.tsv`](config/input_schemas.tsv) records input-table filenames and column names; it contains no data rows. Coordinates in OR interval manifests and `start_1`/`end_1` columns are 1-based inclusive; PAF coordinates are 0-based half-open.

## Software

Use Linux and Bash. The identification workflow supplies a control Conda environment and per-process environments. The versions used in the manuscript are BLAST+ 2.16.0, EMBOSS 6.6.0, MAFFT 7.505, IQ-TREE 2.3.6 for identification, CD-HIT 4.8.1, Exonerate 2.4.0, DIAMOND 2.2.5, trimAI 1.5, IQ-TREE 3.1.3 for family trees, and minimap2 2.24 for regional alignment. InterProScan and SAMtools are also required by the relevant stages. SLURM submission wrappers use the partition and resources in their configurations.

Downstream Python scripts use Biopython; S6 assembly uses Pillow and pypdf:

```bash
python -m pip install biopython Pillow pypdf
```

R plotting packages:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "scales",
                   "patchwork", "jsonlite", "svglite", "ragg"))
```

## 1. Identify OR loci

From the repository root:

```bash
export OR_ROOT="$(pwd)/OR.analysis.scripts"
export OR_DATA_ROOT=/absolute/path/OR_data
export OR_RESULT_ROOT=/absolute/path/OR_results
MINING="$OR_ROOT/01_or_identification/Felidae_OR_mining_release_1.1.3"
cd "$MINING"
bash bootstrap_control_env.sh
conda activate felidae-or-control
cp "$OR_ROOT/config/species/params.Ppa.example.yml" params.Ppa.yml
```

Edit the copied YAML: genome FASTA, OR/GPCR protein databases, reference alignment, TYPE I/II reference ID files, clade-labelled OR database, InterProScan installation, output directory, Conda cache and SLURM settings. Build the matching FASTA index with `samtools faidx`. The output directory's last component must match the analysis-unit code.

```bash
bash bin/run_species_pipeline.sh --params params.Ppa.yml --preflight-only
bash bin/run_species_pipeline.sh --params params.Ppa.yml \
  --profile slurm,conda --from-stage 1 --to-stage 19
```

`initial_search_mode: generate` searches the reference OR database against the assembly before stage 1. Set `reuse` and `input_merged` to use an existing 16-column TBLASTN table. The launcher resumes completed stages. Repeat with the other assembly templates. See the [identification usage guide](01_or_identification/Felidae_OR_mining_release_1.1.3/docs/USAGE_CN.md) and [output guide](01_or_identification/Felidae_OR_mining_release_1.1.3/docs/OUTPUTS_CN.md).

The workflow extends candidate intervals by 1 kb, extracts ORFs, applies similarity and phylogenetic filters, recovers assembly-specific copies and multi-exon models, and writes complete OR catalogs and physical-cluster tables. The default physical-cluster gap is 1 Mb; the PPA distribution plot highlights clusters with at least 50 loci. Sequence-based functional, pseudogene and truncated classifications describe sequence evidence.

## 2. Assign OR subfamilies

Prepare these external inputs:

```text
$OR_DATA_ROOT/
  joint/functional_OR_unique_proteins.prot.fa
  joint/final_joint_OR_representative_assignments.tsv
  joint/functional_OR_identical_sequence_map.tsv
  reference/human_OR_reference.prot.fa
  reference/human_OR_reference.metadata.tsv
```

Use the 417 human OR reference proteins and their HGNC/UniProt family/subfamily metadata. See the schema inventory for table fields. Copy configurations beside their examples so their relative module paths resolve correctly:

```bash
CANDIDATE="$OR_ROOT/02_or_subfamily_assignment/OR_subfamily_candidate_stage_v1.0.2_20260901"
TREES="$OR_ROOT/02_or_subfamily_assignment/OR_family_trees_v1.0.0_20260902"
ASSIGN="$OR_ROOT/02_or_subfamily_assignment/OR_subfamily_tree_assignment_v1.0.0_20260902"
cp "$CANDIDATE/config.example.sh" "$CANDIDATE/config.sh"
cp "$TREES/config.example.sh" "$TREES/config.sh"
cp "$ASSIGN/config.example.sh" "$ASSIGN/config.sh"
bash "$CANDIDATE/bin/run_candidate_stage.sh" "$CANDIDATE/config.sh"
```

The candidate stage writes `candidate_families/family_manifest.tsv`, combined human/felid family FASTAs and `candidate_assignments.tsv`. The family-tree configuration points to these outputs. Run all families locally:

```bash
for list in "$TREES"/lists/*.txt; do
  while IFS= read -r family; do
    [ -n "$family" ] || continue
    bash "$TREES/bin/run_family_tree.sh" "$TREES/config.sh" "$family"
  done < "$list"
done
bash "$TREES/bin/finalize.sh" "$TREES/config.sh"
bash "$ASSIGN/bin/run_assignment.sh" "$ASSIGN/config.sh"
```

For SLURM, use `submit_formal_family_trees.sh` and `submit_assignment.sh` instead of the corresponding local runs. A pre-existing augmented family set can be selected by changing `BRIDGE_ROOT` and `MANIFEST`; its existing completion-marker paths can also be specified.

Family trees use MAFFT, trimAI `-automated1`, IQ-TREE ModelFinder, 1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates. High-confidence defaults include SH-aLRT ≥ 80, UFBoot ≥ 95 and a nearest-alternative distance ratio ≥ 1.10. Final labels also require concordance with similarity evidence; the medium-confidence distance ratio is 1.02. Primary named-subfamily summaries retain high/medium assignments and exclude unresolved loci.

## 3. OR6C regional comparison

```bash
MICRO="$OR_ROOT/03_main_figure/fig5c_or6c_micro_synteny/analysis/OR6C_synteny_audit_v1_20260902"
cp "$MICRO/config.example.sh" "$MICRO/config.sh"
# Edit genome and output paths, then load the configuration.
source "$MICRO/config.sh"
bash "$MICRO/bin/run_OR6C_synteny_audit.sh"
bash "$MICRO/bin/run_OR6C_detailed_alignment.sh"
python "$MICRO/bin/parse_or6c_detailed_paf.py" \
  "$OUT_ROOT/alignments/Ppa_vs_Ppa_or.OR6C.asm5.detailed.paf" \
  "$OUT_ROOT/tables/OR6C.structural_events.min50bp.tsv" \
  "$OUT_ROOT/qc/OR6C.structural_events.summary.json" \
  --minimum-mapq 30 --minimum-event-bp 50
```

The extraction script specifies the study's Ppa ChrB4 and Ppa_or homologous intervals and adds 500-kb flanks. Event parsing retains forward-collinear blocks and reports insertions, deletions and unaligned gaps. `map_or_loci_between_assemblies.py --help` lists the detailed PAF, extracted intervals, subfamily-context table and two OR catalogs needed for locus projection.

The S5 plotting script takes external read-depth and junction-support tables. For these inputs, the manuscript uses Winnowmap `map-pb` with a repetitive 15-mer list, primary alignments with MAPQ ≥ 20 and at least 2 kb of aligned reference sequence on each side, a 4,001-bp junction window and two 10-kb flanks. Keep all alignment classes in regional BAMs used for IGV inspection.

## 4. Chromosome-region comparison

```bash
MACRO="$OR_ROOT/04_supplementary_figures/figs6_chrd1_macrosynteny/analysis_package/OR_ChrD1_macro_synteny_v1_20260903"
cp "$MACRO/config.example.sh" "$MACRO/config.sh"
# Edit config.sh and the genome paths in input/regions.tsv.
source "$MACRO/config.sh"
python "$MACRO/bin/validate_package.py"
bash "$MACRO/bin/run_macro_synteny.sh"
```

`ARCHITECTURE_DIR` contains the seven externally generated `<unit>.OR_operational_clusters.tsv` and `<unit>.OR_locus_architecture.tsv` tables. `input/regions.tsv` contains the study's two-largest-cluster intervals plus 1-Mb flanks; update coordinates and cluster IDs when using another assembly version. `REGIONS_TSV` selects an alternative interval manifest. For SLURM, use `bin/submit_slurm.sh` after loading the configuration.

The summary retains primary blocks ≥ 10 kb, identity ≥ 90% and MAPQ ≥ 20. Coverage is the union of retained intervals; weighted identity is total matching bases divided by total alignment length.

## 5. Figures and copy-number matrices

See [plotting commands and input layout](FIGURES.md). Input-table column names are listed in [`config/input_schemas.tsv`](config/input_schemas.tsv); JSON files and depth files are generated or supplied with the corresponding external analysis outputs.

