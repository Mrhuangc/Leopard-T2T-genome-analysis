# Leopard T2T genome analysis

Analysis scripts for the *Panthera pardus* telomere-to-telomere (T2T) genome project, including genome annotation, assembly comparison, centromere characterization, satellite evolution, chromosome Y analysis, and olfactory-receptor (OR) gene analyses across felids.

[中文说明](README_CN.md)

## Species

| Code | Species |
| --- | --- |
| PPA | *Panthera pardus* — leopard |
| FCA | *Felis catus* — domestic cat |
| PTI | *Panthera tigris* — tiger |
| CFA | Domestic dog, used for synteny comparison |

## Repository structure

```text
*.sh                         Analysis entry scripts
scripts/                     Genome annotation, assembly and centromere helpers
CDR.scripts/                 Methylation, CDR selection, enrichment and plots
SAT.evolution.scripts/       Satellite phase normalization, variants and subtype analysis
ChrY.analysis.scripts/
  annotation.scripts/        Coding-gene and pseudogene annotation
  PAR.scripts/               PAR detection and validation
  partition.scripts/         Y sequence classes and gene-context analysis
  repeat.units.scripts/      Repeat units and pairwise identity
OR.analysis.scripts/         OR identification, subfamilies, synteny and manuscript figures
envs/requirements.txt        Python dependencies
```

## Analysis modules

| Entry script | Analysis | Main outputs |
| --- | --- | --- |
| [repeat.and.gene.annotation.sh](repeat.and.gene.annotation.sh) | EarlGrey repeats; ab initio, homology and RNA-seq gene evidence | Repeat annotations and gene predictions |
| [segment.duplication.annotation.sh](segment.duplication.annotation.sh) | BISER segmental duplication detection, filtering and statistics | Filtered SD pairs, genomic span and chromosome summaries |
| [PUR.analysis.sh](PUR.analysis.sh) | Previously unresolved regions from older assembly/T2T comparison | PUR BED and retained alignments |
| [assembly.err.analysis.sh](assembly.err.analysis.sh) | SyRI assembly comparison and selected rearrangement evidence | Structural variant calls and locus plots |
| [centromere.region.detect.sh](centromere.region.detect.sh) | Centromere prediction, evidence integration and assembly evaluation | Centromere intervals, QV/GCI and evidence figures |
| [SAT_array_pipeline_v2.sh](SAT_array_pipeline_v2.sh) | Satellite monomer and array divergence | Copy/array metrics, phase reports and consensus sequences |
| [HOR.detect.sh](HOR.detect.sh) | HiCAT higher-order repeat analysis | Monomer and HOR annotations |
| [CDR.detect.sh](CDR.detect.sh) | Modified-base processing and centromere dip region analysis | Methylation tracks, Primary/Secondary CDRs and plots |
| [run.sat.evolution.sh](run.sat.evolution.sh) | Cross-species SAT113/SAT483 variants and subtypes | Phase-normalized variants, subtype diagnostics and figures |
| [ChrY.PAR.detect.sh](ChrY.PAR.detect.sh) | X/Y alignment and PAR analysis | PAR candidates and reviewed-interval validation |
| [ChrY.partition.region.sh](ChrY.partition.region.sh) | Y sequence classes and X-derived gene context | Window classes, reviewed partitions and gene-order tables |
| [ChrY.annotation.sh](ChrY.annotation.sh) | Coding-gene curation and pseudogene annotation | GFF3s, validation and pseudogene tables |
| [ChrY.repeat.units.analysis.sh](ChrY.repeat.units.analysis.sh) | Ampliconic and repeat-rich units | Oriented unit sequences, pairwise identity and figures |
| [OR analysis guide](OR.analysis.scripts/README.md) | Felid OR loci, subfamilies, copy number, regional synteny and manuscript figures | OR catalogs, tree-assisted assignments, alignment metrics and plots |

The entry scripts organize the commands used for each analysis. Select the relevant stages and set paths to the corresponding input files and software installations.

## Requirements

Use Linux with Bash and the external tools required by the selected module. Selected workflows use SLURM, Snakemake, or CUDA GPUs.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r envs/requirements.txt
chmod +x scripts/*.pl
export PATH="$(pwd)/scripts:$PATH"
```

The supplied SAT-specific Conda environment can be installed separately:

```bash
conda env create -f SAT.evolution.scripts/environment.yml
conda activate sat_variant_origin
```

R plotting packages:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "patchwork", "scales"))
```

| Analysis | External tools |
| --- | --- |
| Genome annotation | EarlGrey, ANNEVO, Tiberius, Augustus, GeMoMa, fastp, PASA, EVM; HISAT2, Trinity, StringTie, TransDecoder, CD-HIT and PASA aligners for RNA-seq |
| Assembly comparison and SD | minimap2, Winnowmap, meryl, SAMtools, SyRI, BISER, BEDTools |
| Centromeres | TRASH2, Centromics, NGenomeSyn, ModDotPlot, Merqury, GCI, NucFreq, StainedGlass |
| Satellite arrays and HOR | MAFFT, HiCAT, SeqKit, ParaFly |
| Satellite variants | MAFFT, VSEARCH, IQ-TREE for the optional tree stage |
| Methylation and CDRs | pbmm2, modkit, tabix, Dorado, CDR-Finder, Snakemake and Conda |
| ChrY | minimap2, miniprot, mappy, RepeatMasker annotations and a C compiler (`cc`) |

The PTI raw-signal workflow uses Dorado 0.8.2 and the 4 kHz simplex/modified-base models specified in its script. The coding-gene runner compiles the bundled `protein_align_batch.c`.

## Input files and coordinates

Prepare genome FASTA/FAI, reads, repeat and gene annotations, and the interval/evidence files used by the selected module. Examples under `annotation.scripts/config/` and `partition.scripts/` define the ChrY configuration layout.

BED and PAF coordinates are 0-based and half-open; GFF3 coordinates are 1-based and inclusive. Some study interval tables include a header and additional class/evidence columns. Match chromosome identifiers and final assembly coordinates across files.

SAT-array input is a headerless TRASH2 CSV with 11 columns:

```text
"seqID","arrayID","start","end","strand","score","eval","width","class","score_template","sequence"
```

The array entry uses `sequence-orientation as_is`, corresponding to monomer sequences already supplied in their annotation orientation. SAT manifests contain per-species copy tables, sequence sources and phase reports; relative paths resolve against the manifest directory. Relative paths in coding/pseudogene JSON configurations resolve against the JSON directory. Paths in partition YAML configurations refer to the analysis working directory.

## Example commands

### SAT-array analysis

From the repository root:

```bash
INPUT=/absolute/path/PPA.SAT113.repeats_with_seq.csv \
OUTDIR=results/PPA.SAT113.array_pipeline \
PREFIX=PPA.SAT113 THREADS=20 RUN_PLOTS=1 \
bash SAT_array_pipeline_v2.sh
```

The pipeline produces results in `07_results/` and figures in `08_plots/`. Its `--resume` option reuses existing per-array alignments.

### Cross-species satellite analysis

Create the working manifests from the examples and edit the data paths:

```bash
cp SAT.evolution.scripts/sat113_manifest.example.tsv SAT.evolution.scripts/sat113_manifest.tsv
cp SAT.evolution.scripts/sat483_manifest.example.tsv SAT.evolution.scripts/sat483_manifest.tsv
bash run.sat.evolution.sh
```

The SAT113 primary run excludes zero core-edge sites and uses 8 model starts / 2 CV starts. The sensitivity run excludes five core-edge sites and uses 24 model starts / 4 CV starts. Both use up to 18 subtypes and five CV folds. The SAT483 run uses up to 12 subtypes, 8 model starts / 2 CV starts, five CV folds and five excluded core-edge sites.

### ChrY annotation

Copy and edit the JSON templates in place:

```bash
cp ChrY.analysis.scripts/annotation.scripts/config/FCA.example.json ChrY.analysis.scripts/annotation.scripts/config/FCA.json
cp ChrY.analysis.scripts/annotation.scripts/config/PPA.example.json ChrY.analysis.scripts/annotation.scripts/config/PPA.json
cp ChrY.analysis.scripts/annotation.scripts/config/config.standalone.example.json ChrY.analysis.scripts/annotation.scripts/config/config.json
bash ChrY.annotation.sh
```

Coding curation can use reference proteins independently of the other target species' new annotation. An existing cross-species GFF/FASTA reference can also be supplied using the optional `reference_species`, `reference_gff` and `reference_genome` fields together.

Reviewed ChrY partitions can be generated before coding annotation using:

```bash
python3 ChrY.analysis.scripts/partition.scripts/reclassify_xdeg_ampliconic.py --partitions-only
```

These reviewed intervals are written under `ChrY.analysis.scripts/results/04_comparative/` and used by the annotation templates. Candidate-window analysis uses the PPA partition YAML, with `PARTITION_CONFIG` available to select a working configuration. Gene-context analysis uses the corresponding coding annotations and comparative evidence tables.

Centromere and PAR boundaries include study curation. Assembly rearrangements are evaluated with read-support tracks and IGV inspection. Keep these steps aligned with the study methods.

## Olfactory-receptor analyses

The OR module covers seven assemblies from six felid species. `Ppa` denotes Ppa_T2T (PPA), and `Ppa_or` denotes PpoHapG. See the [OR guide](OR.analysis.scripts/README.md) for software environments, input configuration and commands, or the [Chinese guide](OR.analysis.scripts/README_CN.md). Genome assemblies, reads, protein databases and study result tables are supplied separately.
