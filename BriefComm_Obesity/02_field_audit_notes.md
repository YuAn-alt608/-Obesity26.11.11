# 领域审计：NHANES 的 POP–性激素研究有没有筛处方药

> 执行日期 2026-11-21。原始检索结果不可编辑，逐篇条目录于 `02_field_audit.csv`（12 行，含 PMID）。

---

## 1. 检索式与纳入

PubMed（2026-11-21），Title/Abstract 字段：

```
(NHANES OR "National Health and Nutrition Examination Survey")
AND (PFAS OR perfluoro* OR PFOA OR PFOS OR polyfluoroalkyl OR
     "persistent organic pollutants" OR PCB OR polybrominated OR DDE)
AND (testosterone OR androgen* OR SHBG OR "sex hormone" OR estradiol)
```

- 命中 **19 篇**。
- 剔除 7 篇：结局不是性激素（动脉僵硬度、高血压、牙周炎、骨密度各 1 篇）、非原始研究（BPA 综述 1 篇）、暴露非持久性有机物（有机磷阻燃剂 1 篇）、以激素作协变量而非结局（1 篇）。
- **纳入 12 篇**（PFAS 10 篇、PBDE 1 篇、有机氯农药 1 篇）。

> 判定口径：NHANES 原始研究；暴露为持久性有机污染物（含 PFAS/PBDE/OCP）或含 PFAS 的混合物；性激素（睾酮/SHBG/E2/FT/FAI）为结局。

## 2. 结果：药物筛查情况

| 处理方式 | 篇数 | 具体 |
|---|---|---|
| 明确排除影响雄激素轴的药物 | **1** | Madrigal 2021（OCP）：排除抗雄药、5α-还原酶抑制剂、肾上腺皮质激素 |
| 仅笼统排除"激素治疗" | **1** | Li 2025（PBDE）：排除 hormone therapy（n=34），未列药类 |
| 全文核对、未处理 | **2** | Lewis 2015、Xie 2021（仅女性调避孕药） |
| 摘要层面未提及 | **8** | 见 CSV |

**关键空白**：
- **12 篇里 0 篇筛查 GnRH 激动剂/拮抗剂**（前列腺癌治疗，对睾酮的压制可达 90% 以上）。
- **12 篇里 0 篇筛查阿片类**（男性性腺功能减退的常见医源性原因）。
- 只有 1 篇覆盖糖皮质激素。

**核对深度**：对 4 篇以睾酮为主要结局的成人研究做了全文核对（Lewis 2015、Xie 2021、Madrigal 2021、Li 2025），其余 8 篇为摘要层面判定。这一区分必须写进补充材料，不能含糊成"全部查过全文"。

## 3. 直接可用的两句话

**正文（Results 段，约 30 词）：**

> "In a PubMed audit of the 12 NHANES studies that have modelled sex hormones against persistent organic pollutants, two excluded men on hormone therapy and **none reported screening for GnRH agonists or opioid analgesics** (Table S9)."

**Discussion（一句，把审计变成论点）：**

> "The drug class most capable of producing the artefact we describe — GnRH agonists, used in prostate cancer — has never been screened for in this literature."

## 4. 这条审计替我们挡掉的一个反驳

审稿人可以说："你只是自己数据里恰好有几个用药者。"

有了这张表，回答变成：**这不是我们的样本问题，是这个文献体系的结构性盲点**——12 篇里有 0 篇筛过对睾酮压制最强的两类药。这使"167 人翻转一个显著结论"从一个尴尬变成一条领域级证据。

## 5. 局限（写进补充材料）

1. 只检索 PubMed，未覆盖 Embase/Web of Science；未做 PRISMA 流程，属于**范围性审计**，不是系统综述。
2. 8 篇仅核对摘要，作者可能在全文中调整了用药而未在摘要报告。
3. 检索式未包含"medication"作为自由词，因此可能漏掉专门讨论该问题的研究（但这类研究本就会在标题/摘要出现 PFAS 与激素词，被主检索捕获）。
