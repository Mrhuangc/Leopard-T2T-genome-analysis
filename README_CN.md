# 金钱豹 T2T 基因组项目分析代码

本仓库包含金钱豹 T2T 基因组研究的分析脚本，涵盖重复序列与基因注释、组装比较、着丝粒、卫星序列演化、 ChrY 分析及嗅觉受体分析。完整模块、依赖、输入格式和运行示例见 [英文 README](README.md)。

物种代码：PPA 为金钱豹，FCA 为家猫，PTI 为虎，CFA 为家犬。



## 使用方法

在 Linux 环境中安装相应模块的工具，设置数据与软件路径，再运行相关入口。Python 包列表位于 `envs/requirements.txt`；卫星变异分析的 Conda 环境位于 `SAT.evolution.scripts/environment.yml`。

配置文件路径的基准：SAT manifest 和 ChrY JSON 按配置文件所在目录解析相对路径；partition YAML 中的数据路径按分析工作目录使用。

SAT113 主分析使用 0 个 core-edge 位点排除、8 次模型初始化和 2 次 CV 初始化；敏感性分析使用 5 个位点排除、24 次模型初始化和 4 次 CV 初始化。SAT483 的入口及 manifest 示例已包含在仓库中。

ChrY 人工审阅分区可先通过 `reclassify_xdeg_ampliconic.py --partitions-only` 生成，再用于蛋白编码基因和假基因注释。配置示例在 `ChrY.analysis.scripts/annotation.scripts/config/`。

## 嗅觉受体（OR）分析

OR 模块包含六种猫科动物、七套基因组的 OR 鉴定、亚家族分类、拷贝数、共线性分析及论文绘图脚本。`Ppa` 对应 Ppa_T2T（PPA），`Ppa_or` 对应 PpoHapG。输入配置和运行方法见 [OR 中文说明](OR.analysis.scripts/README_CN.md)。
