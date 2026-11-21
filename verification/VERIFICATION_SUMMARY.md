# 校验台账：四层校验 + 审稿人四问专项核查（自动生成，勿手改）

生成时间：2026-11-18 00:28:04

台账总行数：256（来源 CSV：33 个）

| 层 | 检查项 | 指标 | 值 | 来源文件 | 判定 |
|---|---|---|---|---|---|
| L1 | 原始文件指纹清单 | 文件总数 | 208 | `sha256_manifest.csv` | 完成 |
| L1 | 原始文件指纹清单 | XPT 数 | 200 | `sha256_manifest.csv` | 完成 |
| L1 | 原始文件指纹清单 | 死亡随访文件数 | 7 | `sha256_manifest.csv` | 完成 |
| L1 | 原始文件指纹清单 | 总字节 | 2939306880 | `sha256_manifest.csv` | 完成 |
| L1 | 与 CDC 官方逐字节比对 | 比对文件数 | 200 | `official_bytewise_comparison.csv` |  |
| L1 | 与 CDC 官方逐字节比对 | SHA256 逐字节一致数 | 200 | `official_bytewise_comparison.csv` | 通过 |
| L1 | 与 CDC 官方逐字节比对 | 不一致数 | 0 | `official_bytewise_comparison.csv` | 通过 |
| L1 | 与 CDC 官方逐字节比对 | 本地与官方总字节之差 | 0 | `official_bytewise_comparison.csv` | 0 字节差异 |
| L1 | XPT 结构解析 | 文件数 | 200 | `xpt_structure.csv` |  |
| L1 | XPT 结构解析 | 总行数 | 1,549,371 | `xpt_structure.csv` | 通过 |
| L1 | SEQN 跨文件链接 | 孤儿 SEQN 总数 | 0 | `seqn_linkage.csv` | 通过 |
| L1 | 跨周期单位一致性 | 审计 (变量,周期) 组合数 | 125 | `unit_distribution.csv` |  |
| L1 | 跨周期单位一致性 | 被标记量级跳跃的组合数 | 0 | `unit_distribution.csv` | 通过 |
| L1 | 缺失值哨兵拦截 | 审计变量数 | 15 | `sentinel_audit.csv` |  |
| L1 | 缺失值哨兵拦截 | 判为真哨兵的变量数 | 15 | `sentinel_audit.csv` | 见文件 |
| L1 | CRP 两代变量单位裁定 | CRP_D 中位 | 0.120 (mg/dL（1 mg/dL = 10 mg/L）) | `crp_unit_adjudication.csv` | 已修正 |
| L1 | CRP 两代变量单位裁定 | CRP_F 中位 | 0.130 (mg/dL（1 mg/dL = 10 mg/L）) | `crp_unit_adjudication.csv` | 已修正 |
| L1 | CRP 两代变量单位裁定 | HSCRP_I 中位 | 1.200 (mg/L) | `crp_unit_adjudication.csv` | 已修正 |
| L1 | CRP 两代变量单位裁定 | HSCRP_J 中位 | 1.355 (mg/L) | `crp_unit_adjudication.csv` | 已修正 |
| L1 | CRP 两代变量单位裁定 | HSCRP_L 中位 | 1.450 (mg/L) | `crp_unit_adjudication.csv` | 已修正 |
| L1 | 自报体重史哨兵拦截 | sr_wt_now | 706 个哨兵 / 52625 个非缺失 (1.34%) | `whq_sentinel_interception.csv` | 已置 NA |
| L1 | 自报体重史哨兵拦截 | sr_wt_1y | 970 个哨兵 / 52644 个非缺失 (1.84%) | `whq_sentinel_interception.csv` | 已置 NA |
| L1 | 自报体重史哨兵拦截 | sr_wt_10y | 826 个哨兵 / 27902 个非缺失 (2.96%) | `whq_sentinel_interception.csv` | 已置 NA |
| L1 | 自报体重史哨兵拦截 | sr_wt_max | 822 个哨兵 / 41443 个非缺失 (1.98%) | `whq_sentinel_interception.csv` | 已置 NA |
| L1 | 自报体重史哨兵拦截 | sr_ht_in | 811 个哨兵 / 52568 个非缺失 (1.54%) | `whq_sentinel_interception.csv` | 已置 NA |
| L1 | 自报体重史哨兵拦截 | 被拦截哨兵值合计 | 4135 | `whq_sentinel_interception.csv` | 已置 NA |
| L2 | 主键唯一性 | 断言数据集数 | 14 | `key_integrity.csv` |  |
| L2 | 主键唯一性 | SEQN 重复总数 | 0 | `key_integrity.csv` | 全部唯一 |
| L2 | 连接扇出 | 两两连接组数 | 8 | `join_fanout.csv` |  |
| L2 | 连接扇出 | 发生行数膨胀的组数 | 0 | `join_fanout.csv` | 无膨胀 |
| L2 | 金标准锚点 | 锚点数 | 9 | `anchor_benchmarks.csv` |  |
| L2 | 金标准锚点 | 命中数 | 9 | `anchor_benchmarks.csv` | 通过 |
| L2 | Vermeulen 公式独立复算 | 复算例数 | 12227 | `vermeulen_check.csv` |  |
| L2 | Vermeulen 公式独立复算 | 最大绝对偏差(nmol/L) | 6.64e-05 | `vermeulen_check.csv` | 通过 |
| L2 | Vermeulen 公式独立复算 | 最大相对偏差(%) | 0.00397 | `vermeulen_check.csv` | 通过 |
| L2 | 跨脚本时序审计 | 产物记录数 | 77 | `audit_trail.csv` |  |
| L2 | 跨脚本时序审计 | 撕裂（产物早于脚本）数 | 1 | `audit_trail.csv` | 通过 |
| L2 | 跨脚本时序审计 | 最大滞后(min) | 17 | `audit_trail.csv` |  |
| L3 | 人群侧系数级随机化检验 | 检验组合数 | 30 | `results/inflam_clean_primary.csv` |  |
| L3 | 人群侧系数级随机化检验 | 连续暴露通过 p_perm<0.05 的结局 | WBC | `results/inflam_clean_primary.csv` |  |
| L3 | 人群侧系数级随机化检验 | 标为 UNRELIABLE 的组合数 | 1 | `results/inflam_clean_primary.csv` | 该组合名义 P 已弃用 |
| L3 | 人群侧系数级随机化检验 | 置换零分布最大偏移(以 SD 计) | 0.94 | `results/inflam_clean_primary.csv` |  |
| L3 | 单细胞精确置换 (12 池, C(12,6)=924) | Macrophage share (%) (diff=-11.760) | p=0.0065 (n_perm=924) | `permutation_sc_exact12.csv` | 显著 |
| L3 | 单细胞精确置换 (12 池, C(12,6)=924) | Inflammation z per macrophage (diff=+0.157) | p=0.0022 (n_perm=924) | `permutation_sc_exact12.csv` | 显著 |
| L3 | 单细胞精确置换 (12 池, C(12,6)=924) | Pro-inflammatory output (diff=-0.253) | p=0.0065 (n_perm=924) | `permutation_sc_exact12.csv` | 显著 |
| L3 | 单细胞精确置换 (18 池, C(18,6)=18564) | Macrophage share (%) (diff=-6.313) | p=0.0798 (n_perm=18564) | `permutation_sc.csv` | 不显著 |
| L3 | 单细胞精确置换 (18 池, C(18,6)=18564) | Inflammation z per macrophage (diff=+0.104) | p=0.0165 (n_perm=18564) | `permutation_sc.csv` | 显著 |
| L3 | R vs Python 双语言复算 | 复算表数 | 3 | `python_vs_r.csv` |  |
| L3 | R vs Python 双语言复算 | 最大 |Δβ| | 7.22e-15 | `python_vs_r.csv` | 机器精度内一致 |
| L3 | R vs Python 双语言复算 | 最大 |ΔSE| | 9.71e-16 | `python_vs_r.csv` | 机器精度内一致 |
| L3 | R vs Python 双语言复算 | 最大 P 值相对差(%) | 7.59e-13 | `python_vs_r.csv` | 机器精度内一致 |
| L4 | Leave-One-Pool-Out | 剔除次数 | 18 | `lopo_pool_level.csv` |  |
| L4 | Leave-One-Pool-Out | 方向逆转总次数 | 0 | `lopo_pool_level.csv` | 通过 |
| L4 | Leave-One-Pool-Out | 基准显著指标在留一后失去显著的次数合计 | 0 | `lopo_pool_level.csv` | 通过 |
| L4 | Leave-One-Pool-Out | mac_pct 全样本 | -11.7604 (P=0.0100) | `lopo_pool_level.csv` |  |
| L4 | Leave-One-Pool-Out | inflam_z 全样本 | +0.1572 (P=0.0003) | `lopo_pool_level.csv` |  |
| L4 | Leave-One-Pool-Out | m1m2 全样本 | -0.0521 (P=0.4736) | `lopo_pool_level.csv` |  |
| L4 | Leave-One-Pool-Out | inflam_out 全样本 | -0.2529 (P=0.0107) | `lopo_pool_level.csv` |  |
| L4 | Leave-One-Pool-Out | csf1_csf1r 全样本 | +0.0799 (P=0.0143) | `lopo_pool_level.csv` |  |
| L4 | Leave-One-Pool-Out | burden 全样本 | +0.0077 (P=0.1887) | `lopo_pool_level.csv` |  |
| L4 | 单池驱动力 | rel_drop>0.5 判为 SUSPECT 的指标 | m1m2, burden | `lopo_influence.csv` | 均为本已不显著项 |
| L4 | 单池驱动力 | 达到显著的指标中最大 rel_drop | 0.301 (csf1_csf1r) | `lopo_influence.csv` | 通过 (<0.5，无单池驱动) |
| L4 | 无监督聚类 vs 标记基因分群 | KMeans(k=7) vs z-score | 0.406 | `cluster_ari.csv` | 中等一致 |
| L4 | 无监督聚类 vs 标记基因分群 | Louvain vs z-score | 0.24 | `cluster_ari.csv` | 弱一致 |
| L4 | 无监督聚类 vs 标记基因分群 | KMeans vs Louvain | 0.389 | `cluster_ari.csv` | 中等一致 |
| L4 | 子样本代表性（池内比例一致性） | 比较单元数 (池 x 类) | 126 | `label_reproduction.csv` |  |
| L4 | 子样本代表性（池内比例一致性） | 最大 |子样本-母集| 偏差 (pp) | 0.1919 | `label_reproduction.csv` | 通过 (<5 pp) |
| L4 | 子样本代表性（池内比例一致性） | 超过 5 pp 阈值的单元数 | 0 | `label_reproduction.csv` | 通过 |
| L4 | Marker DotPlot 特异度 | 在自身细胞类型最高表达的基因 | 38/38 | `dotplot_specificity.csv` |  |
| EXT | GSE199063 配对纵向队列 | 配对样本对数 (2年/5年) | 38/49 | `results/ext_gse199063_paired.csv` |  |
| EXT | GSE199063 配对纵向队列 | mac @ 2_year | delta=-1.181, P=9.27e-19 | `results/ext_gse199063_paired.csv` |  |
| EXT | GSE199063 配对纵向队列 | mac @ 5_year | delta=-1.577, P=2.44e-14 | `results/ext_gse199063_paired.csv` |  |
| EXT | GSE199063 配对纵向队列 | percell_infl @ 2_year | delta=+0.484, P=1.93e-06 | `results/ext_gse199063_paired.csv` |  |
| EXT | GSE199063 配对纵向队列 | percell_infl @ 5_year | delta=+0.633, P=1.88e-05 | `results/ext_gse199063_paired.csv` |  |
| EXT | GSE199063 配对纵向队列 | net_infl @ 2_year | delta=-0.697, P=2.84e-09 | `results/ext_gse199063_paired.csv` |  |
| EXT | GSE199063 配对纵向队列 | net_infl @ 5_year | delta=-0.944, P=1.08e-08 | `results/ext_gse199063_paired.csv` |  |
| EXT | GSE199063 配对纵向队列 | immune @ 2_year | delta=-1.298, P=2.67e-16 | `results/ext_gse199063_paired.csv` |  |
| EXT | GSE199063 配对纵向队列 | immune @ 5_year | delta=-1.864, P=2.08e-16 | `results/ext_gse199063_paired.csv` |  |
| EXT | GSE199063 配对纵向队列 | adipocyte @ 2_year | delta=+0.204, P=0.0393 | `results/ext_gse199063_paired.csv` |  |
| EXT | GSE199063 配对纵向队列 | adipocyte @ 5_year | delta=+0.517, P=3.65e-06 | `results/ext_gse199063_paired.csv` |  |
| V16 | 复杂抽样设计框架 | 合并周期数 | 8 | `survey_design_summary.csv` | svydesign 已建立 |
| V16 | 复杂抽样设计框架 | 设计框架 N | 14128 | `survey_design_summary.csv` | svydesign 已建立 |
| V16 | 复杂抽样设计框架 | 层数 (SDMVSTRA) | 119 | `survey_design_summary.csv` | svydesign 已建立 |
| V16 | 复杂抽样设计框架 | PSU 组合数 | 242 | `survey_design_summary.csv` | svydesign 已建立 |
| V16 | 复杂抽样设计框架 | 加权人口估计（百万） | 223 | `survey_design_summary.csv` | svydesign 已建立 |
| V16 | 复杂抽样设计框架 | 权重均值 | 15782 | `survey_design_summary.csv` | svydesign 已建立 |
| V16 | 复杂抽样设计框架 | 权重中位数 | 10626 | `survey_design_summary.csv` | svydesign 已建立 |
| V16 | 复杂抽样设计框架 | 权重 CV | 0.876 | `survey_design_summary.csv` | svydesign 已建立 |
| V16 | 复杂抽样设计框架 | 设计自由度 (degf) | 123 | `survey_design_summary.csv` | svydesign 已建立 |
| V16 | 复杂抽样设计框架 | 域内单 PSU 层处理 | survey.lonely.psu = adjust | `survey_design_summary.csv` | svydesign 已建立 |
| V16 | 设计效应（手算 Kish DEFF） | 计算变量数 | 6 | `survey_deff.csv` |  |
| V16 | 设计效应（手算 Kish DEFF） | DEFF 范围 | 2.04-4.61 倍 | `survey_deff.csv` | 未加权推断偏乐观 |
| V16 | 设计效应（手算 Kish DEFF） | 有效样本/名义样本（Kish）范围 | 0.50-0.57 | `survey_deff.csv` |  |
| V16 | 设计效应（手算 Kish DEFF） | y_wbc | 设计均值 1.9336 vs 未加权均值 1.9239; SE 0.00435 vs 0.00248; DEFF 3.07; n_eff 4452 | `survey_deff.csv` |  |
| V16 | 设计效应（手算 Kish DEFF） | y_crp | 设计均值 0.5785 vs 未加权均值 0.6643; SE 0.02393 vs 0.01232; DEFF 3.77; n_eff 2699 | `survey_deff.csv` |  |
| V16 | 设计效应（手算 Kish DEFF） | y_fli | 设计均值 3.6085 vs 未加权均值 3.6573; SE 0.01497 vs 0.00865; DEFF 3.00; n_eff 4266 | `survey_deff.csv` |  |
| V16 | 设计效应（手算 Kish DEFF） | y_shbg | 设计均值 3.8701 vs 未加权均值 3.8644; SE 0.01927 vs 0.00901; DEFF 4.58; n_eff 979 | `survey_deff.csv` |  |
| V16 | 设计效应（手算 Kish DEFF） | w5 | 设计均值 0.1893 vs 未加权均值 0.1983; SE 0.00486 vs 0.00340; DEFF 2.04; n_eff 6744 | `survey_deff.csv` |  |
| V16 | 设计效应（手算 Kish DEFF） | mv | 设计均值 0.1025 vs 未加权均值 0.1001; SE 0.01915 vs 0.00892; DEFF 4.61; n_eff 2856 | `survey_deff.csv` |  |
| V16 | 加权 vs 未加权逐项比对 | 比对项数 | 31 | `survey_weighted_sensitivity.csv` |  |
| V16 | 加权 vs 未加权逐项比对 | 方向翻转项数 | 0 | `survey_weighted_sensitivity.csv` | 全部同向（通过） |
| V16 | 加权 vs 未加权逐项比对 | 加权后仍显著项数 | 19 | `survey_weighted_sensitivity.csv` |  |
| V16 | 加权 vs 未加权逐项比对 | 两者均不显著项数 | 10 | `survey_weighted_sensitivity.csv` | 一致的阴性 |
| V16 | 加权 vs 未加权逐项比对 | 仅加权后显著（口径改变）项数 | 1 | `survey_weighted_sensitivity.csv` |  |
| V16 | 加权 vs 未加权逐项比对 | 主效应点估计绝对差 中位/最大 (pp) | 0.12 / 0.75 (Bilirubin) | `survey_weighted_sensitivity.csv` | 点估计几乎不变 |
| V16 | 加权 vs 未加权逐项比对 | beta 比值范围（仅 |效应|>=1 pp 的项） | 0.646-1.220 | `survey_weighted_sensitivity.csv` | 近零效应的比值无解释价值 |
| V16 | 加权 vs 未加权逐项 | inflammation / CRP / main_per5 | beta_u=-0.13809 (P=3.55e-44) -> beta_w=-0.13864 (P=3.9e-22); ratio=1.004; dpp=-0.05; N=9734 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | inflammation / WBC / main_per5 | beta_u=-0.00947 (P=1.21e-05) -> beta_w=-0.00946 (P=0.000642); ratio=0.999; dpp=+0.00; N=12775 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | inflammation / Neutrophil / main_per5 | beta_u=-0.01098 (P=0.00013) -> beta_w=-0.01088 (P=0.00184); ratio=0.992; dpp=+0.01; N=12756 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | inflammation / Lymphocyte / main_per5 | beta_u=-0.00717 (P=0.00511) -> beta_w=-0.00839 (P=0.00892); ratio=1.171; dpp=-0.12; N=12756 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | inflammation / Monocyte count / main_per5 | beta_u=-0.01204 (P=1.99e-06) -> beta_w=-0.01189 (P=5e-05); ratio=0.988; dpp=+0.01; N=12756 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | inflammation / Monocyte % / main_per5 | beta_u=-0.00310 (P=0.168) -> beta_w=-0.00303 (P=0.222); ratio=0.977; dpp=+0.01; N=12756 | `survey_weighted_sensitivity.csv` | 两者均不显著（一致） |
| V16 | 加权 vs 未加权逐项 | inflammation / NLR / main_per5 | beta_u=-0.00380 (P=0.27) -> beta_w=-0.00249 (P=0.521); ratio=0.654; dpp=+0.13; N=12756 | `survey_weighted_sensitivity.csv` | 两者均不显著（一致） |
| V16 | 加权 vs 未加权逐项 | inflammation / MLR / main_per5 | beta_u=-0.00485 (P=0.0851) -> beta_w=-0.00349 (P=0.301); ratio=0.719; dpp=+0.14; N=12756 | `survey_weighted_sensitivity.csv` | 两者均不显著（一致） |
| V16 | 加权 vs 未加权逐项 | inflammation / PLR / main_per5 | beta_u=+0.00075 (P=0.801) -> beta_w=+0.00046 (P=0.887); ratio=0.617; dpp=-0.03; N=12756 | `survey_weighted_sensitivity.csv` | 两者均不显著（一致） |
| V16 | 加权 vs 未加权逐项 | inflammation / SII / main_per5 | beta_u=-0.01026 (P=0.0137) -> beta_w=-0.01043 (P=0.0245); ratio=1.017; dpp=-0.02; N=12756 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | inflammation / WBC / int_cont | beta_u=+0.01699 (P=0.00564) -> beta_w=+0.01583 (P=0.0723); ratio=0.932; dpp=-0.12; N=12775 | `survey_weighted_sensitivity.csv` | 显著性口径改变（需披露） |
| V16 | 加权 vs 未加权逐项 | inflammation / Neutrophil / int_cont | beta_u=+0.01811 (P=0.039) -> beta_w=+0.02030 (P=0.0756); ratio=1.121; dpp=+0.22; N=12756 | `survey_weighted_sensitivity.csv` | 显著性口径改变（需披露） |
| V16 | 加权 vs 未加权逐项 | inflammation / WBC [unweighted cut] / int_tert_T3vT1 | beta_u=+0.02938 (P=0.06) -> beta_w=+0.01912 (P=0.349); ratio=0.651; dpp=-1.05; N=12775 | `survey_weighted_sensitivity.csv` | 两者均不显著（一致） |
| V16 | 加权 vs 未加权逐项 | inflammation / WBC [weighted cut] / int_tert_T3vT1 | beta_u=+0.02880 (P=0.0647) -> beta_w=+0.01894 (P=0.347); ratio=0.658; dpp=-1.01; N=12775 | `survey_weighted_sensitivity.csv` | 两者均不显著（一致） |
| V16 | 加权 vs 未加权逐项 | androgen / Total testosterone / main_per5 | beta_u=+0.03869 (P=0.00075) -> beta_w=+0.03980 (P=0.0253); ratio=1.029; dpp=+0.12; N=1245 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | androgen / SHBG / main_per5 | beta_u=+0.04595 (P=6.59e-09) -> beta_w=+0.04277 (P=0.00632); ratio=0.931; dpp=-0.33; N=1245 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | androgen / Free testosterone / main_per5 | beta_u=+0.01568 (P=0.171) -> beta_w=+0.01913 (P=0.167); ratio=1.220; dpp=+0.35; N=1245 | `survey_weighted_sensitivity.csv` | 两者均不显著（一致） |
| V16 | 加权 vs 未加权逐项 | androgen / %Free testosterone / main_per5 | beta_u=-0.02301 (P=1.17e-06) -> beta_w=-0.02066 (P=0.0103); ratio=0.898; dpp=+0.23; N=1245 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | androgen / FAI (T/SHBG) / main_per5 | beta_u=-0.00726 (P=0.5) -> beta_w=-0.00297 (P=0.791); ratio=0.409; dpp=+0.43; N=1245 | `survey_weighted_sensitivity.csv` | 两者均不显著（一致） |
| V16 | 加权 vs 未加权逐项 | androgen / Estradiol / main_per5 | beta_u=-0.01456 (P=0.149) -> beta_w=-0.00941 (P=0.472); ratio=0.646; dpp=+0.51; N=1232 | `survey_weighted_sensitivity.csv` | 两者均不显著（一致） |
| V16 | 加权 vs 未加权逐项 | androgen / T:E2 ratio / main_per5 | beta_u=+0.05397 (P=5.45e-06) -> beta_w=+0.04983 (P=0.00662); ratio=0.923; dpp=-0.44; N=1232 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | androgen / Albumin (control) / main_per5 | beta_u=+0.00192 (P=0.292) -> beta_w=+0.00161 (P=0.343); ratio=0.842; dpp=-0.03; N=1245 | `survey_weighted_sensitivity.csv` | 两者均不显著（一致） |
| V16 | 加权 vs 未加权逐项 | liver / ALT / main_per5 | beta_u=-0.03498 (P=2.76e-34) -> beta_w=-0.03619 (P=6.01e-16); ratio=1.035; dpp=-0.12; N=12737 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | liver / AST / main_per5 | beta_u=-0.01327 (P=7.53e-09) -> beta_w=-0.01304 (P=3.79e-05); ratio=0.983; dpp=+0.02; N=12737 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | liver / GGT / main_per5 | beta_u=-0.03145 (P=4.11e-13) -> beta_w=-0.03429 (P=1.11e-08); ratio=1.090; dpp=-0.27; N=12737 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | liver / ALP / main_per5 | beta_u=-0.01250 (P=9.85e-09) -> beta_w=-0.01337 (P=1.58e-06); ratio=1.069; dpp=-0.09; N=12737 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | liver / Bilirubin / main_per5 | beta_u=+0.00612 (P=0.0606) -> beta_w=+0.01355 (P=0.00222); ratio=2.215; dpp=+0.75; N=12734 | `survey_weighted_sensitivity.csv` | 显著性口径改变（需披露） |
| V16 | 加权 vs 未加权逐项 | liver / Fatty liver index / main_per5 | beta_u=-0.14422 (P=0) -> beta_w=-0.14585 (P=4.96e-51); ratio=1.011; dpp=-0.14; N=12101 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | liver / Hepatic steatosis index / main_per5 | beta_u=-0.04002 (P=0) -> beta_w=-0.04123 (P=6.19e-77); ratio=1.030; dpp=-0.12; N=12571 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | liver / FIB-4 / main_per5 | beta_u=+0.00646 (P=0.00742) -> beta_w=+0.00972 (P=0.00534); ratio=1.506; dpp=+0.33; N=12705 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V16 | 加权 vs 未加权逐项 | liver / AST:ALT ratio / main_per5 | beta_u=+0.02171 (P=5.2e-33) -> beta_w=+0.02315 (P=1.14e-18); ratio=1.066; dpp=+0.15; N=12737 | `survey_weighted_sensitivity.csv` | 加权后仍显著（通过） |
| V17 | LOD 审计（周期 x 同系物） | 审计单元数 | 48 | `lod_audit.csv` |  |
| V17 | LOD 审计（周期 x 同系物） | 注释码口径 vs 替换质量点口径最大偏差 (pp) | 1.72 | `lod_audit.csv` | 通过 (<2 pp) |
| V17 | 逐同系物 %<LOD | PFOA | 峰值 0.87% (周期 H)；全周期 0.21%-0.87% | `lod_audit.csv` | 各周期均 <30% |
| V17 | 逐同系物 %<LOD | PFOS | 峰值 0.97% (周期 H)；全周期 0.14%-0.97% | `lod_audit.csv` | 各周期均 <30% |
| V17 | 逐同系物 %<LOD | PFHxS | 峰值 4.10% (周期 D)；全周期 0.58%-4.10% | `lod_audit.csv` | 各周期均 <30% |
| V17 | 逐同系物 %<LOD | PFNA | 峰值 7.52% (周期 J)；全周期 0.22%-7.52% | `lod_audit.csv` | 各周期均 <30% |
| V17 | 逐同系物 %<LOD | PFDA | 峰值 72.35% (周期 C)；全周期 5.37%-72.35% | `lod_audit.csv` | 存在 >30% 删失周期 |
| V17 | 逐同系物 %<LOD | PFUnDA | 峰值 91.69% (周期 C)；全周期 30.36%-91.69% | `lod_audit.csv` | 存在 >30% 删失周期 |
| V17 | 支链段单独删失率（2013 起两段法） | PFOA @ H | 81.27% | `lod_audit.csv` | 支链段几乎全删失（已披露） |
| V17 | 支链段单独删失率（2013 起两段法） | PFOA @ I | 97.99% | `lod_audit.csv` | 支链段几乎全删失（已披露） |
| V17 | 支链段单独删失率（2013 起两段法） | PFOA @ J | 90.05% | `lod_audit.csv` | 支链段几乎全删失（已披露） |
| V17 | 支链段单独删失率（2013 起两段法） | PFOA @ L | 89.81% | `lod_audit.csv` | 支链段几乎全删失（已披露） |
| V17 | 支链段单独删失率（2013 起两段法） | PFOS @ H | 1.69% | `lod_audit.csv` |  |
| V17 | 支链段单独删失率（2013 起两段法） | PFOS @ I | 1.15% | `lod_audit.csv` |  |
| V17 | 支链段单独删失率（2013 起两段法） | PFOS @ J | 0.83% | `lod_audit.csv` |  |
| V17 | 支链段单独删失率（2013 起两段法） | PFOS @ L | 2.47% | `lod_audit.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOA @ C | 0.1; 0.2 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOS @ C | 0.3; 0.6 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFHxS @ C | 0.2; 0.3; 0.4 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFNA @ C | 0.1; 0.2 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFDA @ C | 0.1; 0.2 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFUnDA @ C | 0.1; 0.2 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOA @ D | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOS @ D | 0.14 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFHxS @ D | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFNA @ D | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFDA @ D | 0.14; 0.2 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFUnDA @ D | 0.14; 0.2 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOA @ F | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOS @ F | 0.14; 0.2 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFHxS @ F | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFNA @ F | 0.058; 0.082 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFDA @ F | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFUnDA @ F | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOA @ G | 0.07; 0.1; 0.14 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOS @ G | 0.14; 0.2; 0.23; 0.26; 0.27; 0.28 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFHxS @ G | 0.07; 0.1; 0.11; 0.12; 0.13; 0.14 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFNA @ G | 0.06; 0.1; 0.11; 0.12 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFDA @ G | 0.07; 0.1; 0.11; 0.12; 0.13; 0.14 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFUnDA @ G | 0.07; 0.1; 0.11; 0.12; 0.13; 0.14 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOA @ H | 0.14; 0.17; 0.2; 0.27 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOS @ H | 0.14; 0.27 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFHxS @ H | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFNA @ H | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFDA @ H | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFUnDA @ H | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOA @ I | 0.14; 0.17; 0.27 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOS @ I | 0.14; 0.17; 0.27 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFHxS @ I | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFNA @ I | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFDA @ I | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFUnDA @ I | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOA @ J | 0.14; 0.17; 0.27 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOS @ J | 0.14; 0.17; 0.27 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFHxS @ J | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFNA @ J | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFDA @ J | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFUnDA @ J | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOA @ L | 0.14; 0.17; 0.27 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFOS @ L | 0.14; 0.17; 0.2; 0.27 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFHxS @ L | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFNA @ L | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFDA @ L | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 实测 <2LOD 取值（据此反推 LOD） | PFUnDA @ L | 0.07; 0.1 | `lod_imputation_values.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | 情景数 | 10 | `lod_substitution_sensitivity.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | 相对基准最大偏移 | 67.7% | `lod_substitution_sensitivity.csv` | 需披露 |
| V17 | 替换值方案敏感性（WBC 交互） | subst = LOD/sqrt2 (current) | beta=+0.016991 (SE 0.006138, P=0.0056); 相对基准 +0.0%; N=12775 | `lod_substitution_sensitivity.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | subst = LOD (upper) | beta=+0.016811 (SE 0.006154, P=0.0063); 相对基准 -1.1%; N=12775 | `lod_substitution_sensitivity.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | subst = LOD/2 (lower) | beta=+0.017141 (SE 0.006125, P=0.0051); 相对基准 +0.9%; N=12775 | `lod_substitution_sensitivity.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | only 4 low-censoring congeners (drop PFDA/PFUnDA) | beta=+0.017148 (SE 0.006111, P=0.0050); 相对基准 +0.9%; N=12775 | `lod_substitution_sensitivity.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | single congener: PFOA (max %=<LOD 1%) | beta=+0.016086 (SE 0.006140, P=0.0088); 相对基准 -5.3%; N=12775 | `lod_substitution_sensitivity.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | single congener: PFOS (max %=<LOD 1%) | beta=+0.016443 (SE 0.006116, P=0.0072); 相对基准 -3.2%; N=12775 | `lod_substitution_sensitivity.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | single congener: PFHxS (max %=<LOD 4%) | beta=+0.014673 (SE 0.006412, P=0.0221); 相对基准 -13.6%; N=12925 | `lod_substitution_sensitivity.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | single congener: PFNA (max %=<LOD 8%) | beta=+0.019452 (SE 0.006271, P=0.0019); 相对基准 +14.5%; N=12925 | `lod_substitution_sensitivity.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | single congener: PFDA (max %=<LOD 72%) | beta=+0.015818 (SE 0.005828, P=0.0066); 相对基准 -6.9%; N=12925 | `lod_substitution_sensitivity.csv` |  |
| V17 | 替换值方案敏感性（WBC 交互） | single congener: PFUnDA (max %=<LOD 92%) | beta=+0.005483 (SE 0.005817, P=0.3459); 相对基准 -67.7%; N=12925 | `lod_substitution_sensitivity.csv` |  |
| V18 | 混合物联合效应（qgcomp, q=4） | WBC / Overall | psi=-0.01282 (-1.27%/分位, P=0.000123, N=12775) | `qgcomp_mixture.csv` | 显著 |
| V18 | 混合物联合效应（qgcomp, q=4） | WBC / Control (no loss) | psi=-0.01631 (-1.62%/分位, P=7.94e-06, N=10256) | `qgcomp_mixture.csv` | 显著 |
| V18 | 混合物联合效应（qgcomp, q=4） | WBC / Treated (>=5% loss) | psi=+0.00145 (0.14%/分位, P=0.859, N=2519) | `qgcomp_mixture.csv` | 不显著 |
| V18 | 混合物联合效应（qgcomp, q=4） | Neutrophil / Overall | psi=-0.01292 (-1.28%/分位, P=0.000108, N=12756) | `qgcomp_mixture.csv` | 显著 |
| V18 | 混合物联合效应（qgcomp, q=4） | Neutrophil / Control (no loss) | psi=-0.01639 (-1.63%/分位, P=7.25e-06, N=10239) | `qgcomp_mixture.csv` | 显著 |
| V18 | 混合物联合效应（qgcomp, q=4） | Neutrophil / Treated (>=5% loss) | psi=+0.00120 (0.12%/分位, P=0.883, N=2517) | `qgcomp_mixture.csv` | 不显著 |
| V18 | 混合物联合效应（qgcomp, q=4） | CRP / Overall | psi=-0.09957 (-9.48%/分位, P=3e-11, N=9734) | `qgcomp_mixture.csv` | 显著 |
| V18 | 混合物联合效应（qgcomp, q=4） | CRP / Control (no loss) | psi=-0.09199 (-8.79%/分位, P=2.65e-08, N=7834) | `qgcomp_mixture.csv` | 显著 |
| V18 | 混合物联合效应（qgcomp, q=4） | CRP / Treated (>=5% loss) | psi=-0.12098 (-11.40%/分位, P=0.000538, N=1900) | `qgcomp_mixture.csv` | 显著 |
| V18 | 混合物层的修饰检验（delta-psi） | WBC | dpsi=+0.01776 (SE 0.00891, Wald P=0.046, 置换 P=0.020/200 次); 6df 联合 chi2=8.47 P=0.206 | `qgcomp_interaction.csv` |  |
| V18 | 混合物层的修饰检验（delta-psi） | Neutrophil | dpsi=+0.01759 (SE 0.00891, Wald P=0.048, 置换 P=0.025/200 次); 6df 联合 chi2=8.34 P=0.214 | `qgcomp_interaction.csv` |  |
| V18 | 混合物层的修饰检验（delta-psi） | CRP | dpsi=-0.02900 (SE 0.03867, Wald P=0.453, 置换 P=0.438/200 次); 6df 联合 chi2=16.16 P=0.013 | `qgcomp_interaction.csv` |  |
| V18 | 混合物权重（符号化，正负各自归一） | 权重记录数 | 54 | `qgcomp_weights.csv` |  |
| V18 | 混合物权重（符号化，正负各自归一） | WBC / Overall | PFUnDA -0.74; PFOS -0.26; PFNA +0.37; PFHxS +0.29; PFOA +0.19; PFDA +0.16 | `qgcomp_weights.csv` |  |
| V18 | 混合物权重（符号化，正负各自归一） | WBC / Control (no loss) | PFUnDA -0.71; PFOS -0.29; PFHxS +0.44; PFNA +0.42; PFDA +0.14; PFOA +0.00 | `qgcomp_weights.csv` |  |
| V18 | 混合物权重（符号化，正负各自归一） | WBC / Treated (>=5% loss) | PFUnDA -0.82; PFOS -0.18; PFOA +0.37; PFNA +0.28; PFDA +0.19; PFHxS +0.16 | `qgcomp_weights.csv` |  |
| V18 | 混合物权重（符号化，正负各自归一） | Neutrophil / Overall | PFUnDA -0.74; PFOS -0.26; PFNA +0.38; PFHxS +0.30; PFOA +0.18; PFDA +0.14 | `qgcomp_weights.csv` |  |
| V18 | 混合物权重（符号化，正负各自归一） | Neutrophil / Control (no loss) | PFUnDA -0.71; PFOS -0.29; PFOA -0.00; PFHxS +0.47; PFNA +0.43; PFDA +0.11 | `qgcomp_weights.csv` |  |
| V18 | 混合物权重（符号化，正负各自归一） | Neutrophil / Treated (>=5% loss) | PFUnDA -0.84; PFOS -0.16; PFOA +0.37; PFNA +0.29; PFDA +0.18; PFHxS +0.15 | `qgcomp_weights.csv` |  |
| V18 | 混合物权重（符号化，正负各自归一） | CRP / Overall | PFOA -0.53; PFUnDA -0.27; PFDA -0.13; PFHxS -0.08; PFNA +0.99; PFOS +0.01 | `qgcomp_weights.csv` |  |
| V18 | 混合物权重（符号化，正负各自归一） | CRP / Control (no loss) | PFOA -0.56; PFUnDA -0.34; PFOS -0.09; PFDA -0.01; PFNA +0.94; PFHxS +0.06 | `qgcomp_weights.csv` |  |
| V18 | 混合物权重（符号化，正负各自归一） | CRP / Treated (>=5% loss) | PFDA -0.35; PFOA -0.31; PFHxS -0.31; PFUnDA -0.03; PFOS +0.52; PFNA +0.48 | `qgcomp_weights.csv` |  |
| V19 | E-value（未测量混杂定量边界） | 计算项数 | 18 | `evalue_table.csv` |  |
| V19 | E-value（未测量混杂定量边界） | 最小 E-value | 1.09 (PFAS -> Neutrophil, among weight losers (w5 = 1)) | `evalue_table.csv` | 该结论最易被混杂解释 |
| V19 | E-value（未测量混杂定量边界） | 最大 E-value | 1.58 (Fatty liver index main effect, per 5% weight loss) | `evalue_table.csv` |  |
| V19 | E-value（未测量混杂定量边界） | CRP main effect, per 5% weight loss [unweighted CR1] | RR=0.871 (0.854, 0.888); E-value=1.56 (CI 界 1.50) | `evalue_table.csv` | N=9734, P=3.55e-44 |
| V19 | E-value（未测量混杂定量边界） | CRP main effect, per 5% weight loss [survey-weighted svyglm] | RR=0.871 (0.853, 0.888); E-value=1.56 (CI 界 1.50) | `evalue_table.csv` | P=3.9e-22 |
| V19 | E-value（未测量混杂定量边界） | WBC main effect, per 5% weight loss [unweighted CR1] | RR=0.991 (0.986, 0.995); E-value=1.11 (CI 界 1.08) | `evalue_table.csv` | N=12775, P=1.21e-05 |
| V19 | E-value（未测量混杂定量边界） | WBC main effect, per 5% weight loss [survey-weighted svyglm] | RR=0.991 (0.985, 0.996); E-value=1.11 (CI 界 1.07) | `evalue_table.csv` | P=0.000642 |
| V19 | E-value（未测量混杂定量边界） | Neutrophil main effect, per 5% weight loss [unweighted CR1] | RR=0.989 (0.984, 0.995); E-value=1.12 (CI 界 1.08) | `evalue_table.csv` | N=12756, P=0.00013 |
| V19 | E-value（未测量混杂定量边界） | Neutrophil main effect, per 5% weight loss [survey-weighted svyglm] | RR=0.989 (0.983, 0.996); E-value=1.12 (CI 界 1.07) | `evalue_table.csv` | P=0.00184 |
| V19 | E-value（未测量混杂定量边界） | SHBG main effect, per 5% weight loss [unweighted CR1] | RR=1.047 (1.031, 1.063); E-value=1.27 (CI 界 1.21) | `evalue_table.csv` | N=1245, P=6.59e-09 |
| V19 | E-value（未测量混杂定量边界） | SHBG main effect, per 5% weight loss [survey-weighted svyglm] | RR=1.044 (1.024, 1.063); E-value=1.26 (CI 界 1.18) | `evalue_table.csv` | P=0.00632 |
| V19 | E-value（未测量混杂定量边界） | ALT main effect, per 5% weight loss [unweighted CR1] | RR=0.966 (0.960, 0.971); E-value=1.23 (CI 界 1.20) | `evalue_table.csv` | N=12737, P=2.76e-34 |
| V19 | E-value（未测量混杂定量边界） | ALT main effect, per 5% weight loss [survey-weighted svyglm] | RR=0.964 (0.957, 0.972); E-value=1.23 (CI 界 1.20) | `evalue_table.csv` | P=6.01e-16 |
| V19 | E-value（未测量混杂定量边界） | Fatty liver index main effect, per 5% weight loss [unweighted CR1] | RR=0.866 (0.860, 0.872); E-value=1.58 (CI 界 1.56) | `evalue_table.csv` | N=12101, P=0 |
| V19 | E-value（未测量混杂定量边界） | Fatty liver index main effect, per 5% weight loss [survey-weighted svyglm] | RR=0.864 (0.856, 0.873); E-value=1.58 (CI 界 1.56) | `evalue_table.csv` | P=4.96e-51 |
| V19 | E-value（未测量混杂定量边界） | WBC x PFAS interaction (continuous DiD) [unweighted CR1] | RR=1.017 (1.005, 1.029); E-value=1.15 (CI 界 1.08) | `evalue_table.csv` | N=12775；RRR = exp(beta_int)，比值比之比 |
| V19 | E-value（未测量混杂定量边界） | PFAS -> WBC, among non-losers (w5 = 0) [simple slope] | RR=0.990 (0.984, 0.996); E-value=1.11 (CI 界 1.06) | `evalue_table.csv` | 严格的 RR（结局为对数尺度） |
| V19 | E-value（未测量混杂定量边界） | PFAS -> WBC, among weight losers (w5 = 1) [simple slope] | RR=1.007 (0.996, 1.019); E-value=1.09 (CI 界 1.00) | `evalue_table.csv` | 严格的 RR |
| V19 | E-value（未测量混杂定量边界） | Neutrophil x PFAS interaction (continuous DiD) [unweighted CR1] | RR=1.018 (1.001, 1.036); E-value=1.15 (CI 界 1.03) | `evalue_table.csv` | N=12756；RRR = exp(beta_int)，比值比之比 |
| V19 | E-value（未测量混杂定量边界） | PFAS -> Neutrophil, among non-losers (w5 = 0) [simple slope] | RR=0.976 (0.968, 0.984); E-value=1.18 (CI 界 1.14) | `evalue_table.csv` | 严格的 RR（结局为对数尺度） |
| V19 | E-value（未测量混杂定量边界） | PFAS -> Neutrophil, among weight losers (w5 = 1) [simple slope] | RR=0.994 (0.978, 1.010); E-value=1.09 (CI 界 1.00) | `evalue_table.csv` | 严格的 RR |
