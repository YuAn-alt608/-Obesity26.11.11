# 投稿包清单与合规校验

> 目标刊：**Obesity — Brief Cutting Edge Report**
> 初版校验 2026-11-21；本轮全文数据复核与格式修正在同日完成。
> 数字均为实测（文档脚本 `code/build_docx.py`，图脚本 `code/make_figures_bc.py`）。

---

## 1. 文件清单

| 类别 | 文件 | 说明 |
|---|---|---|
| 主文件 | `manuscript/Main_Manuscript.docx` | 标题页 + 摘要 + Study Importance + 正文 + Table 1 + 图注 + 18 条文献；双倍行距、连续行号、页码、真上标引文 |
| 源文本 | `manuscript/main_manuscript.md` | 与 docx 同源，改这一份再重跑脚本，避免版本漂移 |
| 投稿信 | `manuscript/Cover_Letter.docx` | 单作者署名，含"与完整报告关系"的主动说明（防 salami 质疑） |
| 图 1 | `figures/fig1_artefacts.png` / `.pdf` | 三面板：混杂 / 处方药 / 公式 |
| 图 2 | `figures/fig2_androgen.png` / `.pdf` | 两面板：激素森林图 / 三柱对照 |
| 补充材料 | `supplement/Supplementary_Material.docx` | S1–S9 |
| 补充表 | `supplement/TableS1/S2/S3/S7/S8/S9_*.csv` | 由项目结果直接导出，未手工录入 |
| STROBE | `supplement/STROBE_Checklist.docx` + `.csv` | 22 条逐项对应到章节 |
| 图脚本 | `code/make_figures_bc.py` | 图中每个数字的来源写在脚本注释里 |
| 文档脚本 | `code/build_docx.py` | 从 md 生成四份 docx |
| 领域审计 | `02_field_audit_notes.md` + `02_field_audit.csv` | 12 篇 NHANES 研究逐篇核对 |

---

## 2. 硬指标实测

| 项目 | Obesity 上限 | 本稿实测 | 状态 |
|---|---|---|---|
| 正文词数（Intro→Discussion，不含摘要/文献/图表） | 1,500 | **1,386** | ✅ 余 114 |
| 摘要词数（结构化，含 Objective/Methods/Results/Conclusions 四个标签） | 200 | **196** | ✅ |
| 参考文献 | 20 | **18** | ✅ |
| 图 + 表 | 3 | **3**（2 图 + 1 表） | ✅ |
| 标题字符数 | 125 | **121** | ✅ |
| 标题形式 | 禁问句 | 陈述句 | ✅ |
| Study Importance 三问 | 必需 | 已写（3×34–37 词） | ✅ |
| 正文文件内是否含图 | 不得含 | 不含 | ✅ |
| person-first 语言 | 强制 | 全文无 "obese people/subjects/the obese"；唯一 "obese individuals" 出现在文献 1 的原题中（不得改写引文） | ✅ |
| STROBE 清单 | 必需 | 已附 | ✅ |
| 标题页必备项 | 作者/单位/通讯/字数/图表数/基金/利益冲突 | 均已填 | ✅ |

---

## 3. 关键数字的最终口径（与源文件一致）

| 正文表述 | 数值 | 源 |
|---|---|---|
| 粗估 PhenoAge 差 | +1.87 年（= 1.35 − 0.52） | 报告 §3.3 |
| 调整后 | −0.33（森林 ATE）/ −0.52（参数） | `table2_main.csv` / 报告 §3.3 |
| 位移 | 森林 **2.20 年**（1.87→−0.33）；参数 **2.39 年**（1.87→−0.52） | 相减 |
| 每 5% 减重 PhenoAge | −0.302 年，P=3.9×10⁻⁷ | 报告 §3.4 |
| 用药者 | 167 / 2,244 = 7.45% | 报告 §3.8 |
| CATE 梯度 | −3.4~−4.6 pp → +0.1~+0.8 pp | 报告 §3.8 |
| 交互 P | 0.012–0.091（5 个种子）→ 0.383–0.670 | `forest_seed_stability.csv` |
| HSI / FLI 拆解 | 1.2×10⁻⁶→0.68；0.017→0.14；log ALT 0.61 | `liver_decomp.csv` |
| 30 次检验 | 5→3→0 | `inflam_clean_primary.csv` 实算 |
| 激素面板（N=1,232） | TT +4.02 (1.70, 6.39)；SHBG +4.78 (3.15, 6.43)；%free −2.31 (−3.21, −1.40)；FT +1.62 (−0.65, 3.94)；Alb +0.20；T:E2 +5.55；E2 −1.45 | `table3_androgen_axis.csv` |
| 脱钩检验 | TT–SHBG P=0.83；FAI–SHBG P=0.91；SHBG–albumin P=0.52 | `androgen_decouple.csv` |
| 雄激素×PFAS 交互 | 0.165–0.661 | `table3_androgen_axis.csv` |
| 混合暴露 | 关节项 P=0.261（TT）/ 0.668（FT） | `table2_main.csv` |
| 死亡 | RR 1.23（最低暴露）→1.10（最高） | `mort_w5_pfas.csv` |
| 结构 | 口袋净电荷 −1；白蛋白 +14；Kd=3 nM 时占据率 PFOA 1.6% / PFOS 4.3%（均 <5%） | `t3_pocket_charge.csv` / `t4_occupancy.csv` |

---

## 4. 本轮全文复核：改掉的 12 处（含 5 处会伤及结论的）

| # | 位置 | 原来的写法 | 问题 | 现在 |
|---|---|---|---|---|
| 1 | Results 伪影一 | "moved the estimate 2.20 years … −0.33（森林）and −0.52（参数）" | 2.20 只对应森林；1.87→−0.52 是 2.39 | 两个位移分别写明 2.20 与 2.39 |
| 2 | Results 首段 / 摘要 / cover letter | 把 P=0.0049 称为 "interaction" | 该值是**森林异质性** P；参数交互项 P=0.261。全篇批的就是口径不严，这里不能含糊 | 改为 "forest heterogeneity P=0.0049"，摘要作 "(forest P=0.0049)" |
| 3 | Results 伪影一 | "…and reproduced within cycles (Table S2)" | 报告 §3.10 明写**分周期不能稳定重复**（仅 2/6 周期显著）；且 Table S2 是交互审计表，不是分周期表 | 删除该分句 |
| 4 | Methods 数据 | "Eight NHANES cycles (2003–2023)" | 实际 8 个周期**不含 2007–2008**（`DEMO_*` 只有 C/D/F/G/H/I/J/L） | 改为 "between 2003 and 2023 (2007–2008 was not included)" |
| 5 | Discussion 局限 | "lack drug names for 2017–2018" | 2017–2020（J）的 RXQ_RX **有** RXDDRUG；真正缺药品名的是 **2021–2023（L）**，而该周期本就无睾酮 | 改为 2021–2023，并注明该周期无睾酮测量 |
| 6 | Methods / Results | 1,232 与 2,244 并存无解释 | 审稿人必问"男性样本到底多少人" | Methods 写明 2,244 人为雄激素森林样本，其中 1,232 人有 SHBG/白蛋白/雌二醇全panel |
| 7 | Discussion 首段 | "Each of the three mechanisms…"；Table 1 却有四行 | 三 vs 四自相矛盾 | 表题改 "Three artefacts, one specification error…"，末行标注 "(not an exposure artefact)"，正文补一句"规范错误同样偏向有害" |
| 8 | S3 抽样效率 | "effective sample sizes of 38–57%" | 按 `n_eff/n = 1/deff` 实为 **22–33%**（979–4,452 人）；57% 是仅权重的 Kish 值 | 改为 22–33%，并补一句胆红素跨过显著性阈值的披露 |
| 9 | S4 占据率 | "4% occupancy would require Kd of 30 nM" | Kd=30 nM 时 PFOA 仅 0.17%；4% 需 **≈1 nM** | 改为 ≈1 nM，并补 PFOS 的 4.3%/1.3%/0.14% |
| 10 | S5 组成 | "only macrophages changed significantly" 后又写脂肪细胞 +40.7% 显著 | 自相矛盾（脂肪细胞 BH-P=0.019） | 改为 "only the macrophage fraction **fell** significantly"，脂肪细胞单列并注明分母效应 |
| 11 | S7 种子 | "across three random seeds"（雄激素） | 实际跑了 **5 个**种子；PhenoAge 三个种子的范围是 1.2×10⁻⁷–1.0×10⁻⁴ | 改为 5 个种子，并改正 PhenoAge 范围 |
| 12 | S6 台账 | "Every number … is linked"；PhenoAge "within the analytic sample" | 绝对化表述难兑现；死亡验证是 N=37,853 的全成人样本 | 改为 "Each reported estimate…"；验证样本改为全成人样本 |

另外两处体例改动：全文第一人称由 **we/our 改为 I**（单作者）；引文上标改为 Word 真上标（原先是用 Unicode 上标字符拼的，排版会随字体漂移）。

---

## 5. 还需要你补的三项

1. **数据/代码仓库 DOI**。现写"可向通讯作者索取，接收后存入公共仓库"。若你能开 Zenodo/OSF，把它换成 DOI 更强——这篇的方法学主张靠可复现性支撑。
2. **3 位推荐审稿人**（Obesity 投稿系统第 5 步要求；cover letter 里已不写，留系统填）。
3. **是否挂预印本**（若挂，须在 cover letter 写明服务器与编号）。

---

## 6. 两处需要你确认的表述选择

**（1）标题**。现用 A 案（121 字符）。若编辑偏好更临床，B 案 99 字符：`Total testosterone rose 4% while the free fraction fell 2.3%: what "toxin rebound" alarms get wrong`。两案均不含问句。

**（2）"临床该警惕"的分寸**。正文落点是"总睾酮高估了获益"，并明写游离睾酮绝对值无变化（+1.6%，P=0.16）。**不建议**改成"生物利用度下降造成危害"——那会从读数问题滑成因果主张，审稿人一问"游离睾酮不是升了吗"就站不住。现稿用 300→312 ng/dL 的例子落到临床，是能守住的最强版本。

---

## 7. 尚未做、但可加

- 图 2b 只画三根柱，恒等式写在图注里；如要更直观，可加第三面板（不增加 display item 数）。
- 补充材料现为 Word 正文 + 独立 CSV；若期刊要求单文件上传，可把 CSV 并入 docx（脚本可改）。
- 主报告 `figures/fig14` 面板 d 的 "Kd ≤ 30 nM" 标注与 S4 同源，属同一处高估，需同步改为 ≈1 nM 后重绘（不影响本次投稿文件）。
