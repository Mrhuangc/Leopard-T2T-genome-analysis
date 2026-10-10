# OR 鉴定输出

每套组装使用独立的输出目录。TSV 用于位点及分类结果，FASTA 用于序列，GFF3 用于注释；summary JSON 记录参数和计数。

| 目录 | 主要内容 |
| --- | --- |
| `00_initial_search/` | 查询分片及合并的初始 TBLASTN 表 |
| `03_OR_Finder_results/` | 候选区域、ORF、GPCR 筛选和系统发育过滤 |
| `04_species_OR_expansion/` | 组装自身功能 OR 的二次搜索 |
| `05_multiexon_OR/` | 多外显子候选、Exonerate 模型和分类 |
| `06_OR_catalog/pre7tm/` | 7TM 筛选前的合并目录 |
| `07_7tm_validation/` | InterProScan 结构域结果 |
| `08_final_OR_catalog/` | 结构域筛选后的目录 |
| `09_NC_final_homology/` | 最终同源性复核 |
| `10_single_exon_nonfunctional/` | 单外显子非功能序列恢复 |
| `11_complete_OR_catalog/` | 完整 OR 目录、GFF3 和序列 |
| `12_OR_clade_classification/` | OR 类群和 Class I/II 结果 |
| `13_single_species_OR_analysis/` | 位点结构、物理簇和染色体计数 |

`functional`、`pseudogene`、`truncated` 和 `ambiguous` 为序列证据分类。物理簇按相邻位点间隔不超过 1 Mb 的链式规则确定。

`start_1/end_1` 及 GFF3 为 1-based inclusive；BED 为 0-based half-open。转换为 BED 时起点减 1，终点不变。

