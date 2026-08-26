# Zenodo 13736891 SecAct 配对治疗动态方案 v0.1

**冻结日期：** 2026-08-24  
**冻结时点：** 已核实数据结构、论文总体设计及 Figure 2B 存在 CR/non-CR 注释；尚未人工重建 9 位 bulk 患者的逐例疗效，也未运行 SecAct  
**性质：** 小样本、患者级、探索性治疗动态分析；不作为 HRA006718 或化疗反应预测的替代验证

## 一、研究问题

在接受 camrelizumab、apatinib 和 TPC 联合诱导治疗的鼻咽癌患者中，一周期治疗后的肿瘤/转移淋巴结 SecAct 分泌蛋白活动是否相对治疗前发生患者内一致变化？

## 二、数据和统计单位

- 数据文件：Zenodo DOI 10.5281/zenodo.13736891 中的 `BulkRNAseq.rar/BulkRNAseq.xls`。
- 9 位患者，每位患者一份 BT 和一份 OT 样本，共 18 个表达列。
- 统计单位为患者，n=9；18 个样本不作为 18 个独立个体。
- 本版本不使用 CR/non-CR 标签，不开展疗效分层。

## 三、预处理

1. 读取实际为 tab-separated 文本的 `BulkRNAseq.xls`。
2. 使用 `gene_symbol`，删除缺失符号；只保留 SecAct signature 中出现的人类基因。
3. 若保留基因符号重复，按患者样本列取均值；当前审计显示非空 `gene_symbol` 唯一。
4. 所有表达值必须为非负且无缺失；转换为 `log2(x + 1)`。
5. 由列名解析 patient ID、BT/OT，并要求 9 位患者各有且仅有一对样本。

## 四、SecAct 和推断参数

- SecActPy 版本：运行时记录。
- `is_differential=False`。
- `n_rand=1000`。
- RNG：GSL。
- seed：20260823。
- gene sorting：启用。
- 对相关系数 >=0.9 的签名按 SecActPy `group_signatures` 合并，预期检验 657 个独立签名组。

对每个签名组计算 OT - BT 的患者内差值，并报告：

- mean/median paired difference；
- 配对 Cohen dz 及小样本校正后的 Hedges gz；
- 配对 t 检验；
- Wilcoxon signed-rank 检验；
- 512 种完整 sign-flip 的双侧精确 P 值；
- 精确 P 的 BH FDR；
- 2,000 次按患者对重采样的 mean-difference 95% CI 和同向率。

## 五、探索性动态组门槛

同时满足：

- sign-flip BH FDR <0.10；
- |Hedges gz| >=1.0；
- bootstrap 同向率 >=0.90。

该门槛只定义“值得后续机制核实的治疗动态组”，不定义临床生物标志物。若 0 个组通过，按阴性结果报告，不改用 CR/non-CR、单侧检验或原始 P 值救结果。

## 六、敏感性和边界

- 对排名靠前的组同时报告配对 t、Wilcoxon 和 exact sign-flip，不选择最有利检验。
- 不因某位患者影响大而事后删除；可报告 leave-one-patient-out 方向稳定性，但不改变正式门槛。
- 不使用单细胞的 87,191 个细胞扩大样本量。
- 本分析评估的是联合免疫化疗后的动态，不能单独归因于 PD-1、apatinib、紫杉类、顺铂或卡培他滨中的任一成分。
- 与 GP 诱导化疗不是同一治疗情境，不能用于验证 GSE132112 或 HRA006718 的 GP 疗效预测。

## 七、停止规则

- 0 个动态组：结束全 SecAct 筛选，仅保留完整结果和总体描述。
- 有动态组：只允许进入预先说明的患者级单细胞来源/受体定位，不开展结果后全库追靶。
- 无论结果如何，不在本数据中建立分类器、选截点或报告交叉验证 AUC。
