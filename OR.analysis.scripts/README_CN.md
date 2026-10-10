# 猫科嗅觉受体（OR）分析

本模块包含 OR 鉴定、亚家族分类、拷贝数比较、区域共线性及论文绘图脚本，并提供软件环境和输入配置。基因组、HiFi reads、参考蛋白序列及分析数据快照单独准备。

分析包含六种猫科动物、七套组装。OR 脚本使用区分大小写的代码：`Fca`（家猫）、`Lge`（乔氏猫）、`Pbe`（豹猫）、`Ple`（狮）、`Ppa`（Ppa_T2T/PPA）、`Ppa_or`（PpoHapG）和 `Pti`（苏门答腊虎）。PPA 使用最终组装坐标及已调整方向的 ChrY。

## 输入配置

- `config/species/params.<unit>.example.yml`：七套组装的 OR 鉴定参数模板。复制至鉴定工作流目录，填写基因组、参考数据库、InterProScan、输出和集群路径。
- `config/species_metadata.tsv`：分析单元和物种信息。
- `config/input_schemas.tsv`：外部输入表的文件名和列名，不包含数据行。
- 亚家族三个子模块的 `config.example.sh`：候选输入、家族树和最终赋值路径。
- OR6C 与染色体共线性模块的 `config.example.sh`：基因组、区域及结果路径；染色体区域配置见 `input/regions.tsv`。

## 运行顺序

从仓库根目录设置路径：

```bash
export OR_ROOT="$(pwd)/OR.analysis.scripts"
export OR_DATA_ROOT=/absolute/path/OR_data
export OR_RESULT_ROOT=/absolute/path/OR_results
```

1. 在 `01_or_identification/Felidae_OR_mining_release_1.1.3/` 创建并激活控制环境，逐组装运行鉴定流程，得到 OR 目录、类群和物理簇。
2. 在 `02_or_subfamily_assignment/` 依次运行候选筛选、家族树、亚家族赋值。示例配置已经衔接三个模块的结果路径。
3. 对 OR6C 区域执行序列提取、详细比对、事件解析和位点投影；对两个最大 OR 簇执行染色体区域比较。
4. 用外部分析结果运行主图、补充图及拷贝数脚本。

完整命令、输入布局和依赖版本见 [英文使用说明](README.md)；绘图命令见 [FIGURES.md](FIGURES.md)。鉴定流程的简要说明和输出目录见 [USAGE_CN.md](01_or_identification/Felidae_OR_mining_release_1.1.3/docs/USAGE_CN.md) 和 [OUTPUTS_CN.md](01_or_identification/Felidae_OR_mining_release_1.1.3/docs/OUTPUTS_CN.md)。

