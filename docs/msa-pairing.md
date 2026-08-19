# Paired MSA 与 Unpaired MSA

面向复合物结构预测的 MSA 组织方式说明。结论部分均已在真实数据上核对，
数据来源见文末「证据」一节。

## 1. 为什么需要 MSA

AlphaFold 一类模型从多序列比对里读的是**共进化信号**：如果第 i 列和第 j 列的
突变总是同时发生，说明这两个残基在三维空间中彼此接触——一个位置发生了替换，
另一个位置必须跟着补偿，否则结构或功能被破坏。单体折叠依靠的就是这种链内相关性。

## 2. 到了复合物，问题变了

预测两条链能否结合、以什么姿态结合，需要的是**跨链**的共进化：A 链第 i 位与
B 链第 j 位一起变化，提示这两个残基构成界面接触。

但要计算跨链相关性，同一行里的 A 段和 B 段**必须来自同一个物种**。如果这一行
A 段取自大肠杆菌、B 段取自人，把它们并排比较得到的只是噪音——两者从未在同一个
细胞里共同演化过。

于是 MSA 被拆成两块：

```
             ← A 的列 →      ← B 的列 →
paired    [ A@小鼠      ][ B@小鼠      ]   ← 同物种，携带跨链信号
          [ A@斑马鱼    ][ B@斑马鱼    ]
          [ A@酵母      ][ B@酵母      ]

unpaired  [ A的某同源物 ][ ---gap---   ]   ← 块对角，只有链内信号
          [ A的某同源物 ][ ---gap---   ]
          [ ---gap---   ][ B的某同源物 ]
          [ ---gap---   ][ B的某同源物 ]
```

- **paired（配对）**：一行内两段来自同一物种，携带界面信息，决定"能不能对上"。
- **unpaired（非配对）**：块对角排列，A 的命中在 B 的列上全填 gap，反之亦然。
  只携带链内信息，决定"每条链自己折得对不对"。

两块拼进同一个特征矩阵一起喂给模型，缺一不可。

## 3. 关键认识：配的是同源物的物种，不是查询序列的物种

这是最容易误解的一点。查询的两条蛋白自己属于什么物种**不重要**，重要的是
它们各自的同源物能否在同一个物种里同时找到。

| 场景 | 配对结果 | 原因 |
|---|---|---|
| 人蛋白 A × 人蛋白 B | 配对很深，可达上千行 | 两者在小鼠、斑马鱼、果蝇、酵母都有直系同源物 |
| 人蛋白 × 病毒蛋白 | 配对近乎为空 | A 的同源物在脊椎动物，B 的同源物在病毒株，没有物种同时含两者 |

**这正是宿主—病原互作预测格外困难的根本原因**：不是算法不行，是配对 MSA
本身就是空的，模型只能靠先验去猜界面。

### 对全配对筛选的直接后果

跨物种对与同物种对的 ipTM / ipSAE **不可直接放在一张表里比大小**。同物种对
天然配对更深、分数更高，会被误读为"更可能互作"。

因此筛选结果表必须附带一列 **paired MSA 行数（配对深度）**，用于分层解读。
这一列的计算成本几乎为零，但缺了它，整个 N×N 筛选的结论都可能是配对深度的
假象而非真实的互作信号。

## 4. AlphaFold2-Multimer 的配对算法

来源：`alphafold/data/msa_pairing.py`

1. 按物种标识分组（`msa_species_identifiers_all_seq`，取自 UniProt 那一路 MSA）
2. **只在一条链中出现的物种直接丢弃**（原文：跳过只出现于单条链的物种）
3. 同一物种内若两条链各有多个命中，按「与各自查询序列的相似度」降序排列，
   再按名次逐位配对——第 1 名配第 1 名，第 2 名配第 2 名
4. 配对必须使用**带物种标注的数据库（UniProt）**；unpaired 可以使用任意库，
   包括没有可靠物种归属的宏基因组库（BFD / MGnify / colabfold_envdb）

第 4 条是两块 MSA 必须分开检索、分开缓存的根本原因——不是实现上的取舍，
而是数据来源本身不同。

ColabFold 提供两种策略，对应上面第 2 步的松紧：

- `pairgreedy`：允许部分行（某条链在该物种缺失时留 gap）
- `paircomplete`：只保留所有链都有命中的物种（AF2 的原始行为）

两条链时差别不大；三条链以上差别显著。

## 5. 证据：真实数据中的头部格式

在 A100 上已有的猪蛋白组 MSA（`AniVirusInteractome`，UP000008227 / taxid 9823，
22,804 条序列）上核对：

**Protenix MSA 服务返回的 `pairing.a3m` —— 头部带 taxid**

```
>query
FPAAARHRRGLQTCSRYQTLELEKEFQCNPYLTCKRWIEVSHALGLTERQIKIWFQNRRMKWKKREQ...
>UniRef100_UPI0008F9D4CD_7038/	102	0.351	2.228E-22	0	90	93	112	201	286
>UniRef100_K0A1Z2_1902835/	100	0.355	7.911E-22	2	90	93	105	194	311
>UniRef100_A0A7M5VAS7_252671/	99	0.303	1.491E-21	0	91	93	14	112	224
```

格式为 `UniRef100_<accession>_<taxid>/` 后接制表符分隔的比对统计量
（score、seqid、evalue、qstart、qend、qlen、tstart、tend、tlen）。

抽样统计（条目 268）：7428 条头部中 7427 条带 taxid（唯一的例外是 `>query`），
涉及 1694 个不同物种；同一物种最多出现 50 次——这正是 AF2 需要「组内按相似度
排序」这一步的原因。换条目复核（条目 145）：1443 / 1444，模式一致。

**同一条序列的 `non_pairing.a3m` —— 不带 taxid**

```
>query
FPAAARHRRGLQTCSRYQTLELEKEFQCNPYLTCKRWIEVSHALGLTERQIKIWFQNRRMKWKKREQ...
>ERR1719370_563322	87	0.311	3.413E-17	1	90	93	25	111	114
```

`ERR` 开头是 ENA 测序运行编号，属宏基因组来源，没有可靠物种归属。与第 4 条
完全吻合。

**ColabFold API 的 `uniref.a3m`（`mode=env`，非配对）—— 不带 taxid**

```
>UniRef100_A0A8D1GZC5	923	0.991	1.444E-294	2	738	739	25	761	762
```

只有 accession，没有 taxid 后缀。

## 6. 对工具设计的意义

`pairing.a3m` 是**逐条序列**产出的，且自带 taxid。这意味着：

- 配对 MSA 不必按复合物向服务器提交。**每条唯一序列提交一次**，把带 taxid 的
  命中缓存下来，之后任意两条序列的配对都可以在本地按 taxid 完成。
- N 条序列的全配对筛选（N(N+1)/2 个组合）所需的 MSA 检索次数是 **O(N)**，
  不是 O(N²)。100 条序列的 5050 个组合，只需 100 次检索。

注意：ColabFold 客户端自带的 `pair_sequences()`（`colabfold/input.py`）是按
**行号**做拼接的，不做 taxid 匹配——它依赖服务器端已经把各链区块按物种对齐。
因此若要走上面的 O(N) 路线，需要自行实现第 4 节那套按 taxid 分组配对的逻辑，
而不能直接复用 ColabFold 的客户端函数。

## 7. 数据来源速查

| 内容 | 来源库 | 带 taxid | 用途 |
|---|---|---|---|
| `pairing.a3m` | UniRef100 / UniProt | 是 | 配对 MSA，跨链信号 |
| `non_pairing.a3m` | 宏基因组（ENA 等） | 否 | 非配对 MSA，链内信号 |
| `uniref.a3m`（ColabFold `mode=env`） | UniRef30 | 否 | 非配对 MSA |
| `bfd.mgnify30.metaeuk30.smag30.a3m` | BFD / MGnify / MetaEuk / SMAG | 否 | 非配对 MSA |

---

参考：

- `alphafold/data/msa_pairing.py`（AF2-Multimer 配对实现）
- `colabfold/input.py` 中的 `pair_sequences` / `pad_sequences` / `msa_to_str`
- `colabfold/mmseqs/search.py` 中的 `mmseqs_search_pair`（本地 `pairaln` 路径）
