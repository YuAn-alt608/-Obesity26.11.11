## ============================================================================
## build_evidence_pack.py — 把全文用到的支撑证据按编号复制到一个文件夹并建索引
## ----------------------------------------------------------------------------
## 输出：evidence_pack/
##         README.md            （人读总目录，含编号规则）
##         INDEX.csv            （机读索引：编号/幕/章节/论断/文件/来源脚本/关键数字）
##         01_feasibility/ … 07_code/  （按幕分文件夹，文件名前缀 = 编号）
## 原则：
##   1) 一律"复制"不改动原件；原件仍在 results/ verification/ figures/ code/
##   2) 文件名 = <编号>_<原文件名>，编号顺序 = 论证顺序，不是文件名顺序
##   3) 二进制中间产物（.rds / 大 .csv.gz）不复制，只在索引里登记原位路径
## ============================================================================
import os, re, csv, shutil, io, glob

BASE = os.environ.get("AGE22_ROOT", os.getcwd())
PACK = os.path.join(BASE, "evidence_pack")

## 幕（act）→ 文件夹名
ACT_DIR = {
    "1": "01_feasibility",     # 第一幕 能不能做：数据可得性与设计
    "2": "02_main_effects",    # 第二幕 三个稳健主效应
    "3": "03_modification",    # 第三幕 PFAS 修饰：三个突破口
    "4": "04_singlecell",      # 第四幕 单细胞机制
    "5": "05_verification",    # 第五幕 可信度：四层校验 + 四问
    "6": "06_figures",         # 图
    "7": "07_code",            # 代码
    "8": "08_reports",         # 正文与附录
}

## ---------------------------------------------------------------------------
## 证据清单：(编号, 幕, 章节, 支撑的论断, 源文件, 关键数字/备注)
## ---------------------------------------------------------------------------
M = [
## ===== 第一幕：能不能做 =====
("E1-01", "1", "§1", "数据可得性核查：8 周期 200 个 XPT 的结构与行数（DNAm 数据不存在、PCB/OCP 仅 2003-04）",
 "verification/xpt_structure.csv", "200 文件 / 1,549,371 行"),
("E1-02", "1", "§1", "变量字典：183 个 XPT、4,532 个变量的名称/标签/非缺失数",
 "data/meta/vars_all.csv", "任何变量映射可回查"),
("E1-03", "1", "§2.1", "各周期 N 与暴露可得性（ΣPFAS / PCB / OCP / PFOS 分列）",
 "results/qc_pop.csv", "C 周期 10,122 人；PFAS 子样本 2,094"),
("E1-04", "1", "§2.1", "分析样本 S0–S4 的构建质控（排除规则逐条计数）",
 "results/qc_build_core.csv", "S0 28,342 → S1 9,656"),
("E1-05", "1", "§3.1", "表 1：减重组 vs 未减重组基线（N=7,779 / 1,877）",
 "results/table1_baseline.csv", "全部主分析样本"),
("E1-06", "1", "附", "S1 全样本个体级数据（9,656 人 × 关键变量），可脱离脚本直接复算",
 "results/S1_subject_level.csv", "1.1 MB，审稿人复核入口"),

("E1-07", "1", "§1 / 第一幕", "★ 样本量漏斗：82,096 → 45,739 → 16,778 → 13,162 → 12,775（主分析域）",
 "results/cohort_flow.csv", "11 步；每步注明筛选条件，跨 branch 不算百分比"),

## ===== 第二幕：三个稳健主效应 =====
("E2-01", "2", "§3.2/3.9", "PhenoAge 对全因死亡的预测力：每加速 1 年 RR=1.032",
 "results/mort_phenoage.csv", "RR 1.032 (1.029–1.036), P=5.3e-69"),
("E2-02", "2", "§3.2", "PhenoAge 判别力 AUC=0.857，与实足年龄 0.856 等价",
 "logs/13_mort.log", "⚠ 仅日志有，未落盘 CSV（见 INDEX 备注）"),
("E2-03", "2", "§3.3–3.5", "主效应与交互总表（5 个模型：PFAS/PCB/OCP/睾酮/游离睾酮）",
 "results/table2_main.csv", "PhenoAge×PFAS b2=0.35, P=0.116"),
("E2-04", "2", "§3.4", "参数化剂量-反应：转折点、过零点、联合 Wald",
 "results/param_interaction.csv", "b2 联合检验 P=0.116（未确证）"),
("E2-05", "2", "§3.4", "剂量曲面网格（ΔW × PFAS z 的拟合面与置信带原料）",
 "results/surface_grid.csv", "417 KB，供图 3 复现"),
("E2-06", "2", "§3.10", "稳健性：分周期复现",
 "results/robust_bycycle.csv", "8 周期逐一同规格"),
("E2-07", "2", "§3.10", "稳健性：阈值/对照组重定义",
 "results/robust_threshold.csv", "≥5%/≥10%/稳定体重对照"),
("E2-08", "2", "§3.9", "硬终点：减重 × PFAS 对死亡率的放大（阴性）",
 "results/mort_w5_pfas.csv", "RR 1.233 (0.996–1.527) 跨零"),
("E2-09", "2", "§3.7.1", "雄激素轴八结局：每 5% 减重的效应",
 "results/androgen_per5.csv", "TT +4.0%、SHBG +4.8%、%游离 −2.3%、游离 +1.6% (P=0.17)"),
("E2-10", "2", "§3.7.3", "雄激素：PFAS 三分位交互（三种设定）",
 "results/androgen_tertile.csv", "八结局全部阴性"),
("E2-11", "2", "§3.7.2", "脱钩检验：总睾酮与 SHBG 的阻尼是否成比例（配对聚类自助法）",
 "results/androgen_decouple.csv", "阻尼成比例 → 不支持结合位点竞争"),
("E2-12", "2", "§3.7.4", "雄激素：单味同系物模型（PFOA/PFOS/PFHxS/PFNA）",
 "results/androgen_congener.csv", "无亲和力梯度"),
("E2-13", "2", "§3.7.4", "雄激素：多同系物互调（互斥共线性诊断）",
 "results/androgen_congener_multi.csv", "交互仍阴性"),
("E2-14", "2", "§3.8", "雄激素主表（已排除六类药物的主分析样本）",
 "results/androgen_axis.csv", "S3 2,244 → 清洁样本"),
("E2-15", "2", "§3.8", "同一模型的全样本版本（不排除药物）——用于展示药物混杂的量级",
 "results/androgen_axis_allsample.csv", "与 E2-14 并读，差异即药物效应"),
("E2-16", "2", "§3.7", "表 3（雄激素轴）/ 表 4（同系物）投稿版",
 "results/table3_androgen_axis.csv", "配套 table4_congener.csv（见 E2-17）"),
("E2-17", "2", "§3.7.4", "表 4：四味同系物的完整估计",
 "results/table4_congener.csv", ""),
("E2-18", "2", "§3.13.1", "肝-脂九结局：每 5% 减重（本研究最强信号）",
 "results/liver_per5.csv", "ALT −3.4%、GGT −3.1%、AST −1.3%"),
("E2-19", "2", "§3.13.1", "肝-脂：PFAS 三分位交互（表面显著的那一版）",
 "results/liver_tertile.csv", "拆解前"),
("E2-20", "2", "§3.13.1", "肝-脂：四个临床二分类终点（FLI≥60 的 OR=0.226 等）",
 "results/liver_binary.csv", "临床可解释性"),
("E2-21", "2", "§3.13.2", "肝-脂：二次剂量-反应与联合 Wald",
 "results/liver_quad.csv", ""),
("E2-22", "2", "§3.13.3", "肝-脂：分层分析",
 "results/liver_strat.csv", ""),
("E2-23", "2", "§3.13.2", "★ 评分拆解：纯肝脏部分 vs BMI 部分（把假交互拆开的关键一步）",
 "results/liver_decomp.csv", "HSI 肝部分 P=0.68、FLI 肝部分 P=0.14"),
("E2-24", "2", "§3.13.3", "★ 剂量构成诊断：PFAS 是否修饰减重幅度本身",
 "results/liver_weight_diag.csv", "ΔW1y 中位 9.9% vs 8.3%，交互 P=0.49–0.98"),
("E2-25", "2", "§3.14.1", "★ 炎症十结局主效应（修正截距后的权威规格）",
 "results/inflam_clean_primary.csv", "CRP −12.9%、WBC −0.9%、SII −1.0%"),
("E2-26", "2", "§3.14.2", "炎症：三分位交互（三档口径）",
 "results/inflam_clean_tertile.csv", "WBC −2.3%/+0.6%/+0.5%，联合 P=0.127"),
("E2-27", "2", "§3.14.1", "炎症：连续剂量（每 5%）与置换 P",
 "results/inflam_clean_dose.csv", ""),
("E2-28", "2", "§3.14.2", "炎症：总睾酮式脱钩检验（WBC vs 单核的阻尼是否成比例）",
 "results/inflam_decouple.csv", ""),
("E2-29", "2", "§3.14.2", "炎症：剂量构成诊断（PFAS 不修饰减重幅度）",
 "results/inflam_weight_diag.csv", "排除“剂量差异”这一非生物学解释"),
("E2-30", "2", "§3.14", "ΔBMI 实测剂量交互（与自报剂量对照）",
 "results/inflam_dbmi_interaction.csv", ""),
("E2-32", "2", "§3.3", "因果森林主输出：校准检验、ATE、十分位 CATE（ATE=−0.33 年的来源）",
 "logs/10_cf.log", "ATE −0.33 (SE 0.14)；CATE(P90)−CATE(P10)=+0.71 (P=0.082)"),
("E2-31", "2", "§3.14", "v3 旧规格的炎症结果（CRP/WBC 等，无截距误设前）——保留仅供对照，已被 E2-25 取代",
 "results/inflam_per5.csv", "⚠ 不得作为主报依据"),

("E2-33", "2", "§3.7 / 第四幕", "★ 雄激素轴：总量 vs 生物利用度（TT +4.0% 而 %游离 −2.3%）",
 "results/androgen_bioavailability.csv", "TT +4.02%、SHBG +4.78%、%游离 −2.31%、T:E2 +5.55%；白蛋白 +0.20% (ns)"),

## ===== 第三幕：PFAS 修饰，三个突破口 =====
("E3-01", "3", "§3.12", "突破口一：高敏亚组的实测效应（PFOS +26.8%、PFOA +24.2%，符号与假说相反）",
 "results/subgroup_highsens.csv", "n=56，处理组"),
("E3-02", "3", "§3.12", "★ 突破口一的功效量化：80% 功效下 MDE=±35.7%",
 "results/subgroup_power.csv", "任何小于 MDE 的点估计不得解释为“无效应”"),
("E3-03", "3", "§3.5", "★ 森林异质性检验的种子稳定性（只改种子重跑 5 次）",
 "results/forest_seed_stability.csv", "单一 P 值不足以支撑“异质性真实”"),
("E3-04", "3", "§3.5", "PhenoAge 侧森林的种子稳定性",
 "results/forest_seed_phenoage.csv", ""),
("E3-05", "3", "§3.14.2", "★ 突破口二判决：系数级随机化检验（PFAS 在“周期×性别”层内打乱）",
 "verification/permutation_pop_clean.csv", "WBC 置换 P=0.011/0.020；30 次检验唯一存活"),
("E3-06", "3", "§3.14.2", "★ 置换零分布诊断：随机暴露下估计量本身有偏（需位置校正）",
 "verification/permutation_pop_null.csv", "位置校正 P 的口径依据"),
("E3-07", "3", "§3.14.2", "系数稳定性诊断（10 结局 × 3 估计量的稳定性矩阵）",
 "verification/coef_stability_diag.csv", ""),
("E3-08", "3", "§3.13.2", "突破口三：肝酶反跳的表面交互与拆解后归零",
 "results/liver_decomp.csv", "与 E2-23 同一文件，两处引用"),
("E3-09", "3", "§3.13.2", "突破口三：BMI 项产生假交互的机制（剂量构成差异）",
 "results/liver_weight_diag.csv", "与 E2-24 同一文件"),
("E3-10", "3", "§3.6", "PCB/OCP 双轴：参数化交互（P10 −0.98 → P90 +0.45 年）",
 "results/param_interaction.csv", "与 E2-04 同一文件；两条轴互斥不可合并"),
("E3-11", "3", "§3.8", "规格敏感性：纳入/排除药物的两套结果并列",
 "verification/specification_sensitivity.csv", "不隐藏任何一套"),
("E3-12", "3", "§3.8", "药物标志的原始比对样本（含/不含用药者的激素数据切片）",
 "verification/_horm_sample.csv", "配套 _horm_full.csv（全量）"),
("E3-13", "3", "§3.14", "★ 混合暴露：qgcomp 联合效应（整体/未减重/减重三层）",
 "verification/qgcomp_mixture.csv", "WBC −1.62% vs +0.14%"),
("E3-14", "3", "§3.14", "★ 混合暴露：符号化权重（哪一味在驱动）",
 "verification/qgcomp_weights.csv", "WBC 由 PFNA/PFHxS 主导；CRP 由 PFNA +0.94 主导"),
("E3-15", "3", "§3.14", "★ 混合暴露：Δψ + 置换 + 6 df 联合 χ²（含不一致的披露）",
 "verification/qgcomp_interaction.csv", "Δψ=+0.0178，置换 P=0.020，联合 P=0.206"),

## ===== 第四幕：单细胞机制 =====
("E4-01", "4", "§4.2", "组成：巨噬细胞占比 25.3% → 13.5%",
 "results/geo_composition.csv", "18 池 × 7 类"),
("E4-02", "4", "§4.3", "脂肪细胞功能模块（脂解/脂滴/AHR/SASP…）",
 "results/geo_adipocyte_modules.csv", "配套检验表见 E4-03"),
("E4-03", "4", "§4.3", "模块检验（Welch t + BH）",
 "results/geo_adipocyte_module_tests.csv", ""),
("E4-04", "4", "§4.3/4.4", "★ 通路子评分（把脂解酶与脂滴结构蛋白分开，避免相互抵消）",
 "results/geo_adipocyte_subscores.csv", ""),
("E4-05", "4", "§4.4", "★ 子评分检验 + 脂解–外排耦合相关",
 "results/geo_adipocyte_subscore_tests.csv", "耦合不成立 → 假说核心前提不成立"),
("E4-06", "4", "§4.5", "关键基因 log₂FC（WL/Obese，BH 校正）",
 "results/geo_key_genes.csv", "LIPE +0.63、MGLL +0.62"),
("E4-07", "4", "§4.5", "关键基因的检验表",
 "results/geo_key_gene_tests.csv", ""),
("E4-08", "4", "§4.7.3", "M1/M2 极化（不变，P=0.47）",
 "results/sc_m1m2.csv", ""),
("E4-09", "4", "§4.7.2", "★ 净炎症负荷 = 占比 × 每细胞评分（不变，P=0.19）",
 "results/sc_burden.csv", "0.032 vs 0.039"),
("E4-10", "4", "§4.8", "配体-受体：三条关闭、一条打开",
 "results/sc_lr_pairs.csv", "IL6→IL6R −66.8%、CSF1→CSF1R +20.2%"),
("E4-11", "4", "§4.8", "配体/受体覆盖率（通讯分析的可检出性边界）",
 "results/sc_lig_rec_pct.csv", ""),
("E4-12", "4", "§4.2", "细胞类型均值矩阵（fig7 与假体 bulk 的原料）",
 "results/cell_means.csv", ""),
("E4-13", "4", "§4.6", "★ 外部背靠背：GSE199063 个体配对纵向（50/49/38 + 28 对照）",
 "results/ext_gse199063_paired.csv", "mac −1.18/−1.58 SD，P≤2e-14"),
("E4-14", "4", "§4.6", "外部：与对照组比较",
 "results/ext_gse199063_vs_control.csv", ""),
("E4-15", "4", "§4.6", "外部：逐样本模块评分",
 "results/ext_gse199063_scores.csv", ""),

("E4-16", "4", "§4.7.1 / 第二幕", "★ 选择性耗竭：7 类细胞逐一比较（只有巨噬显著下降，其余 6 类 ns）",
 "results/geo_composition_selectivity.csv", "macrophage −46.5% (BH-P=0.035, d=−2.07)；adipocyte +40.7% 属分母效应"),

("E4-17", "4", "§5.3.1 / 第四幕", "★ 口袋电荷对比：SHBG 净 −1（仅 Asp65）vs 白蛋白 +14（6 Arg + 1 Lys）",
 "struct_pfas/t3_pocket_charge.csv", "解释 PFAS 头基（−1）在 SHBG 被排斥、在白蛋白被吸引"),
("E4-18", "4", "§5.3.1 / 第四幕", "配体形状：PFOA vdW 体积 221 Å³ < DHT 237 Å³，排除「位阻装不下」",
 "struct_pfas/t1_ligand_shape.csv", "PFOA 最粗半径 3.64 Å、PFOS 3.59 Å、DHT 3.90 Å"),
("E4-19", "4", "§5.3.1 / 第四幕", "SHBG 口袋沿长轴 20 个切片的可用半径剖面（最窄 4.17 Å）",
 "struct_pfas/t2_pocket_profile.csv", "配体最粗半径全部低于最窄处"),
("E4-20", "4", "§5.3.1 / 第四幕", "★ 占据率敏感性：56 种组合，需 Kd ≤ 30 nM 才够 4% 占据",
 "struct_pfas/t4_occupancy.csv", "Kd=3nM→1.6%；100nM→0.05%；1μM→0.005%"),
("E4-21", "4", "§5.3.1 / 第四幕", "结构分析方法学与 AF3 作业包的判读标准（含三条硬限制）",
 "struct_pfas/README.md", "先定判据再看结果；AF3 只能给 pose，给不了 Kd"),

## ===== 第五幕：可信度 =====
## L1 源头
("E5-01", "5", "C.1", "L1 全部原始文件的 SHA256 指纹清单",
 "verification/sha256_manifest.csv", "208 文件"),
("E5-02", "5", "C.1", "★ L1 与 CDC 官方逐字节比对：200/200 完全一致、0 字节差异",
 "verification/official_bytewise_comparison.csv", "机读可核"),
("E5-03", "5", "C.1", "L1 跨周期单位一致性（量级跳跃检测）",
 "verification/unit_distribution.csv", ""),
("E5-04", "5", "C.1", "L1 缺失值哨兵拦截（9999/7777/999 等）",
 "verification/sentinel_audit.csv", ""),
("E5-05", "5", "C.1", "★ L1 CRP 两代变量单位裁定（mg/dL vs mg/L 的 10 倍断裂）",
 "verification/crp_unit_adjudication.csv", "E-1 勘误的证据"),
("E5-06", "5", "C.1", "L1 自报体重史哨兵拦截（4,135 个哨兵）",
 "verification/whq_sentinel_interception.csv", ""),
("E5-07", "5", "C.1", "L1 人体测量学边界（32–227 kg）",
 "verification/anthropometric_bounds.csv", ""),
("E5-08", "5", "C.1", "L1 SEQN 跨文件链接（孤儿 SEQN 数）",
 "verification/seqn_linkage.csv", ""),
## L2 基准
("E5-09", "5", "C.2", "L2 主键唯一性（各数据集 SEQN 重复数）",
 "verification/key_integrity.csv", ""),
("E5-10", "5", "C.2", "L2 连接扇出（无行数膨胀）",
 "verification/join_fanout.csv", ""),
("E5-11", "5", "C.2", "★ L2 金标准锚点 9/9 命中（文献可查的量级基准）",
 "verification/anchor_benchmarks.csv", ""),
("E5-12", "5", "C.2", "L2 金标准锚点明细",
 "verification/anchor_detail.csv", ""),
("E5-13", "5", "C.2", "★ L2 Vermeulen 游离睾酮公式独立复算（最大偏差）",
 "verification/vermeulen_check.csv", "逐人复算，非抽样"),
("E5-14", "5", "C.2", "L2 跨脚本时序审计（产物是否早于脚本 = 是否伪造）",
 "verification/audit_trail.csv", "撕裂数 0"),
("E5-15", "5", "C.2", "L2 CDC 官方单位表比对",
 "verification/official_units.csv", ""),
## L3 防伪
("E5-16", "5", "C.3", "★ L3 系数级随机化检验（人群侧，修正规格后）",
 "verification/permutation_pop_clean.csv", "与 E3-05 同一文件"),
("E5-17", "5", "C.3", "★ L3 单细胞精确置换（12 池，C(12,6)=924 全枚举）",
 "verification/permutation_sc_exact12.csv", "p=0.0065 等"),
("E5-18", "5", "C.3", "★ L3 单细胞精确置换（18 池，C(18,6)=18,564 全枚举）",
 "verification/permutation_sc.csv", "p=0.080（如实报告不再显著）"),
("E5-19", "5", "C.3", "★ L3 R vs Python 双语言复算（最大 |Δβ|/|ΔSE|/ΔP）",
 "verification/python_vs_r.csv", "机器精度内一致"),
("E5-20", "5", "C.3", "L3 双语言复算的原始两侧结果",
 "verification/python_crosscheck.csv", "配套 python_crosscheck_r.csv"),
## L4 单细胞防伪
("E5-21", "5", "C.4", "★ L4 Leave-One-Pool-Out（18 池，0 次方向逆转）",
 "verification/lopo_pool_level.csv", ""),
("E5-22", "5", "C.4", "L4 单池驱动力（rel_drop>0.5 的仅有本已不显著项）",
 "verification/lopo_influence.csv", ""),
("E5-23", "5", "C.4", "★ L4 无监督聚类 vs 标记基因分群（ARI）",
 "verification/cluster_ari.csv", "KMeans 0.406 / Louvain 0.240"),
("E5-24", "5", "C.4", "L4 聚类纯度（巨噬 0.847）",
 "verification/cluster_purity.csv", ""),
("E5-25", "5", "C.4", "L4 子样本代表性（126 单元，max 0.19 pp）",
 "verification/label_reproduction.csv", ""),
("E5-26", "5", "C.4", "L4 Marker DotPlot 特异度（38/38）",
 "verification/dotplot_specificity.csv", ""),
("E5-27", "5", "C.5", "外部数据完整性（54 个 gzip，0 损坏）",
 "verification/geo_extracted_integrity.csv", ""),
## 四问
("E5-28", "5", "C.7", "① 加权：8 周期合并设计框架（119 层 × 242 PSU，223 百万）",
 "verification/survey_design_summary.csv", ""),
("E5-29", "5", "C.7.3", "① 加权：手算 Kish 设计效应（DEFF 2.04–4.61）",
 "verification/survey_deff.csv", "有效 N 仅名义的 1/3–1/2"),
("E5-30", "5", "C.7.4", "★ ① 加权：31 项逐项比对（0 方向翻转；WBC 交互 P 0.0056→0.072）",
 "verification/survey_weighted_sensitivity.csv", "Supplementary Table S1"),
("E5-31", "5", "C.8", "② LOD：逐周期 %<LOD + 支链段删失（双口径，最大偏差 1.72 pp）",
 "verification/lod_audit.csv", "PFUnDA 全周期 >30%"),
("E5-32", "5", "C.8.3", "② LOD：实测 <2LOD 取值（据此反推 LOD）",
 "verification/lod_imputation_values.csv", ""),
("E5-33", "5", "C.8.5", "★ ② LOD：三方案替换 + 逐同系物（点估计偏移 ≤1.1%）",
 "verification/lod_substitution_sensitivity.csv", "PFUnDA 单独建模失效（P=0.35）"),
("E5-34", "5", "C.9", "③ 混合物：ψ（见 E3-13）",
 "verification/qgcomp_mixture.csv", "同 E3-13"),
("E5-35", "5", "C.9", "③ 混合物：权重（见 E3-14）",
 "verification/qgcomp_weights.csv", "同 E3-14"),
("E5-36", "5", "C.9", "③ 混合物：Δψ（见 E3-15）",
 "verification/qgcomp_interaction.csv", "同 E3-15"),
("E5-37", "5", "C.10", "★ ④ E-value：18 项（WBC 交互 1.15、CRP 1.56、FLI 1.58）",
 "verification/evalue_table.csv", ""),
("E5-38", "5", "C.0.1", "★ 总台账：256 条记录，每行 = 一个数字 + 一个来源 CSV",
 "verification/verification_ledger.csv", "layer 字段 L1–L4 / V16–V19 / EXT"),
("E5-39", "5", "C.0.1", "总台账的人读版（Markdown 表格）",
 "verification/VERIFICATION_SUMMARY.md", ""),

## ===== 图 =====
("E6-01", "6", "图 1", "暴露与体重变化分布（含“整数化”堆积）", "figures/fig1_exposure.png", ""),
("E6-02", "6", "图 2", "主效应 ATE 森林", "figures/fig2_ate.png", ""),
("E6-03", "6", "图 3", "CATE 曲面与剂量-反应", "figures/fig3_cate.png", ""),
("E6-04", "6", "图 4", "睾酮轴剂量曲面", "figures/fig4_testosterone_surface.png", ""),
("E6-05", "6", "图 5", "死亡随访与 PhenoAge 效度", "figures/fig5_mortality.png", ""),
("E6-06", "6", "图 6", "稳健性（分周期/阈值/置换）", "figures/fig6_robustness.png", ""),
("E6-07", "6", "图 7", "单细胞四面板（组成/通路/关键基因/耦合）", "figures/fig7_singlecell.png", ""),
("E6-08", "6", "图 8", "雄激素轴四面板（含药物筛查 8d）", "figures/fig8_androgen_axis.png", ""),
("E6-09", "6", "图 9", "肝-脂轴（含拆解与剂量构成）", "figures/fig9_liver_axis.png", ""),
("E6-10", "6", "图 10", "高敏亚组功效演示", "figures/fig10_subgroup_power.png", ""),
("E6-11", "6", "图 11", "免疫重塑与细胞通讯", "figures/fig11_scrna_remodel.png", ""),
("E6-12", "6", "图 12", "系统炎症与髓系钝化（含加权/混合物对照）", "figures/fig12_inflammation.png", ""),
("E6-13", "6", "图 13", "Marker DotPlot（38 基因 × 7 类）", "figures/fig13_dotplot_markers.png", ""),
("E6-14", "6", "图 14", "PFAS × SHBG 结构不相容性：四面板（形状 / 口袋剖面 / 电荷 / 占据率）", "figures/fig14_pfas_shbg_structure.png", "最窄 4.17 Å vs 配体 3.59–3.64 Å"),

## ===== 报告 =====
("E8-01", "8", "全文", "主报告（v4，含 §0–§8 与可核查性附录）", "report/Full_Report.md", ""),
("E8-02", "8", "附录", "附录 C：可复现性与数据真实性声明（C.0–C.11）",
 "report/Appendix_C_Reproducibility.md", ""),
]

## PDF 版图（若存在）同编号加后缀 _pdf
PDF_EXTRA = [("E6-01","figures/fig1_exposure.pdf"),("E6-02","figures/fig2_ate.pdf"),
             ("E6-03","figures/fig3_cate.pdf"),("E6-04","figures/fig4_testosterone_surface.pdf"),
             ("E6-05","figures/fig5_mortality.pdf"),("E6-06","figures/fig6_robustness.pdf"),
             ("E6-07","figures/fig7_singlecell.pdf"),("E6-08","figures/fig8_androgen_axis.pdf"),
             ("E6-09","figures/fig9_liver_axis.pdf"),("E6-10","figures/fig10_subgroup_power.pdf"),
             ("E6-11","figures/fig11_scrna_remodel.pdf"),("E6-12","figures/fig12_inflammation.pdf"),
             ("E6-13","figures/fig13_dotplot_markers.pdf"),
             ("E6-14","figures/fig14_pfas_shbg_structure.pdf")]

## 关键脚本（单独编号，供"从原始数据到该数字"的一跳）
SCRIPTS = [
 ("E7-01", "§2/§3.1", "人群构建与样本定义（S0–S4）", "code/06_build_pop.R"),
 ("E7-02", "§2.1", "分析数据集与个体级导出", "code/09_prep_analysis.R"),
 ("E7-03", "§3.2–3.5", "因果森林主流程", "code/10_causal_forest.R"),
 ("E7-04", "§3.4", "参数化剂量-反应交互", "code/11_param.R"),
 ("E7-05", "§3.9", "死亡随访（泊松率模型）", "code/13_mortality.R"),
 ("E7-06", "§3.10", "稳健性全套", "code/15_robust.R"),
 ("E7-07", "§3.7–3.8", "雄激素轴全套（八结局/三分位/自助法/同系物）", "code/22_androgen_axis.R"),
 ("E7-08", "§3.8", "处方药标志构建（rx_flags）", "code/23_rx_flags.R"),
 ("E7-09", "§3.13", "肝-脂轴九结局", "code/30_liver_analysis.R"),
 ("E7-10", "§3.13.2", "★ 评分拆解与剂量诊断（本报告最关键的一步）", "code/31_liver_decomp.R"),
 ("E7-11", "§3.12", "亚组功效量化", "code/32_subgroup_highsens.R"),
 ("E7-12", "§4.7–4.8", "单细胞极化与配体-受体", "code/33_singlecell_lr.R"),
 ("E7-13", "§3.14", "★ 系统炎症修正规格（权威交互表）", "code/36_inflam_clean.R"),
 ("E7-14", "§4", "单细胞主流程", "code/14_geo.R"),
 ("E7-15", "§4.4", "通路子评分与耦合检验", "code/20_geo_subscores.R"),
 ("E7-16", "C.1–C.6", "四层校验脚本 V01–V12（目录整体，含 V01b/V10a/V10b 变体）", "code/verify/"),
 ("E7-17", "C.7–C.10", "四问补检脚本 V16–V19（目录整体）", "code/verify/"),
 ("E7-18", "C.0.1", "台账生成器 V15", "code/verify/V15_ledger.R"),
]

## ---------------------------------------------------------------------------
def producer_index():
    """把每个 CSV 反查生成它的脚本（同名 grep），用于索引的 source_script 列"""
    files = glob.glob(os.path.join(BASE, "code", "*.R")) + \
            glob.glob(os.path.join(BASE, "code", "*.py")) + \
            glob.glob(os.path.join(BASE, "code", "verify", "*.*"))
    idx = {}
    for f in files:
        if os.path.basename(f) == "build_evidence_pack.py":
            continue                      # 本脚本自己会提到全部文件名，必须排除
        try:
            txt = io.open(f, encoding="utf-8", errors="ignore").read()
        except Exception:
            continue
        for m in re.finditer(r"([A-Za-z0-9_\-\.]+\.(?:csv|png|pdf|md|rds))", txt):
            idx.setdefault(m.group(1), set()).add(
                os.path.relpath(f, BASE).replace("\\", "/"))
    return idx

def main():
    os.makedirs(PACK, exist_ok=True)
    for d in ACT_DIR.values():
        os.makedirs(os.path.join(PACK, d), exist_ok=True)
    prod = producer_index()
    rows, missing, seen = [], [], {}
    for eid, act, sec, claim, src, note in M:
        p = os.path.join(BASE, src.replace("/", os.sep))
        if not os.path.exists(p):
            missing.append((eid, src)); continue
        base = os.path.basename(src.rstrip("/"))
        if os.path.isdir(p):
            continue
        newname = "%s_%s" % (eid, base)
        dst = os.path.join(PACK, ACT_DIR[act], newname)
        shutil.copy2(p, dst)
        seen.setdefault(eid, 0)
        seen[eid] += 1
        if seen[eid] > 1:                        # 同一编号被多处引用 → 加后缀
            stem, ext = os.path.splitext(newname)
            newname = "%s_%d%s" % (stem, seen[eid], ext)
            os.rename(dst, os.path.join(PACK, ACT_DIR[act], newname))
        ps = sorted(prod.get(base, []))
        if not ps and src.startswith("report/"):
            ps = ["（作者撰写，非脚本生成）"]
        if not ps and re.match(r"^\d+_.*\.log$", base):      # 日志 → 同名编号脚本
            num = base.split("_")[0]
            ps = sorted(os.path.relpath(x, BASE).replace("\\", "/")
                        for x in glob.glob(os.path.join(BASE, "code", num + "_*.R")))
        rows.append(dict(id=eid, act=act, section=sec, claim=claim,
                         file="%s/%s" % (ACT_DIR[act], newname),
                         origin=src, source_script="; ".join(ps) if ps else "（历史/中间产物，无单一生成脚本）",
                         key_number=note))
    ## 图 PDF
    for eid, src in PDF_EXTRA:
        p = os.path.join(BASE, src.replace("/", os.sep))
        if os.path.exists(p):
            base = os.path.basename(src)
            newname = "%s_%s" % (eid, base)
            shutil.copy2(p, os.path.join(PACK, ACT_DIR["6"], newname))
            rows.append(dict(id=eid, act="6", section="图（矢量版）", claim="同编号 PNG 的 PDF 矢量版",
                             file="06_figures/" + newname, origin=src,
                             source_script="code/16_figures.R 等", key_number="投稿排版用"))
    ## 脚本
    for eid, sec, claim, src in SCRIPTS:
        p = os.path.join(BASE, src.replace("/", os.sep))
        if os.path.isdir(p):
            import re as _re
            pat = _re.compile(r"V(0[1-9]|1[0-5])") if eid == "E7-16" else _re.compile(r"V(1[6-9])")
            for f in sorted(glob.glob(os.path.join(p, "*"))):
                if os.path.isfile(f) and pat.search(os.path.basename(f)):
                    nm = "%s_%s" % (eid, os.path.basename(f))
                    shutil.copy2(f, os.path.join(PACK, ACT_DIR["7"], nm))
            rows.append(dict(id=eid, act="7", section=sec, claim=claim,
                             file="07_code/（目录内全部文件，前缀 %s_）" % eid,
                             origin=src, source_script=src, key_number=""))
        elif os.path.exists(p):
            nm = "%s_%s" % (eid, os.path.basename(src))
            shutil.copy2(p, os.path.join(PACK, ACT_DIR["7"], nm))
            rows.append(dict(id=eid, act="7", section=sec, claim=claim,
                             file="07_code/" + nm, origin=src,
                             source_script=src, key_number=""))
        else:
            missing.append((eid, src))

    rows.sort(key=lambda r: (r["act"], r["id"]))
    with io.open(os.path.join(PACK, "INDEX.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id","act","section","claim","file","origin","source_script","key_number"])
        w.writeheader(); w.writerows(rows)
    print("索引条目:", len(rows), " 缺失:", missing)
    return rows

if __name__ == "__main__":
    main()
