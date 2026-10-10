# OR 鉴定流程使用方法

## 安装与输入

在本目录的上一级工作流根目录执行：

```bash
bash bootstrap_control_env.sh
conda activate felidae-or-control
felidae-or version
cp params.template.yml params.Ppa.yml
```

填写参数文件中的基因组 FASTA、80% identity 去冗余 OR 蛋白库、GPCR/chemoreceptor 蛋白库、参考比对、TYPE I/II ID 文件、带类群标签的 OR 参考蛋白库和 InterProScan 路径。生物学输入使用绝对路径；基因组索引通过 `samtools faidx` 生成。`outdir` 最后一级目录须与 `species` 相同。

仓库级 `config/species/` 提供七套组装的模板。物种代码由属名首字母和种加词前两位组成；亚种加下划线及亚种前两位。组装标签单独写入 `assembly_id`。

## 初始搜索

新搜索使用 `initial_search_mode: generate`。已有初始 TBLASTN 表时使用 `reuse`，并填写 `input_merged`。复用表无表头、以制表符分隔，列顺序为：

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen
```

## 执行

```bash
bash bin/run_species_pipeline.sh --params params.Ppa.yml --preflight-only
bash bin/run_species_pipeline.sh --params params.Ppa.yml \
  --profile slurm,conda --from-stage 1 --to-stage 19
```

配置共享 Conda 缓存及 SLURM 的 partition、account 和 QoS。启动器自动启用恢复；`generate` 模式会在第一阶段前运行初始搜索。只运行初始搜索可用 `--initial-only`；只准备 Conda 环境可用 `--prepare-conda-only`。

```bash
bash bin/run_species_pipeline.sh --list-stages
bash bin/run_species_pipeline.sh --params params.Ppa.yml \
  --profile slurm,conda --from-stage 10 --to-stage 19
```

保留 `work/` 和 `.nextflow/` 以恢复已完成任务。单步工作流位于 `tests/nextflow/`；启动器、控制脚本、Python 模块和对应环境共同构成运行流程。

