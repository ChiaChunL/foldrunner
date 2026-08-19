# 各引擎输入格式核对表

每一条都对照上游源码或官方规范核对过,不是转述文档散文。**不要**为了让格式看起来
整齐而"修正"它们——这些不一致是真实存在的。

## 真机验证状态

生成的文件已经在 A100/A40 上喂给各引擎**自己的解析器**跑通,不是只对着我们自己的
测试跑。测试面板:3 条真实序列 → 6 个复合物,MSA 来自缓存。

| 引擎 | 验证器 | 结果 |
|---|---|---|
| AlphaFold 3 | `folding_input.Input.from_json` | 6/6,`unpairedMsaPath`/`pairedMsaPath` 正确解析 |
| AlphaFold Server | `from_alphafoldserver_fold_job`(骨架)+ 官方 schema | 骨架 6/6,`count: 2` → 链 A/B |
| Boltz-2 | `boltz.main.check_inputs` | 6/6 |
| Protenix | `SampleDictToFeatures.get_feature_dict` | 6/6,`count: 2` → 64 token 两条链 |
| Chai-1 | `inference_dataset.read_inputs` | 6/6,同源二聚体 → 2 条记录 |
| ColabFold | `get_queries` + `unserialize_msa` | 目录 6、CSV 6 带 a3m;`#len\tcard` 头正确还原 |
| AlphaFold2-Multimer | `alphafold.data.parsers.parse_fasta` | 2 条记录,链名 B/C,无冒号 |
| SeedFold | 与官方样例结构比对 | 顶层键、entity 键、类型、取值全部一致 |

### 陷阱:不要用 AF3 的转换器去"修正" AlphaFold Server 输出

AF3 本地带一个 `Input.from_alphafoldserver_fold_job`,看起来是 server dialect 的
权威校验器,**它不是**。它只支持服务端能力的一个子集,会拒绝
`unpairedMsa`、`glycans`、`maxTemplateDate`,而这三个在 DeepMind 自己的
`server/README.md` 里都有明确文档(`unpairedMsa` 在第 111 行有说明、第 161 行有示例)。

照着这个转换器去删 `unpairedMsa`,会把唯一能复用 MSA 的网页端引擎变成不能复用。
骨架可以用它验,MSA 字段只能对文档验。

## 速查

| 引擎 | 载体 | 多拷贝表达 | 保留输入链名 | 外部 MSA 接口 |
|---|---|---|---|---|
| AlphaFold 3 | 一案例一 JSON | `id: ["A","B"]` | 是 | `unpairedMsaPath` / `pairedMsaPath`(或内嵌字符串) |
| AlphaFold Server | JSON 数组 | `count: 2` | 否(强制 A/B) | `unpairedMsa`(只能内嵌) |
| Boltz-2 | 一案例一 YAML | `id: [A, B]` | 是 | `msa: <path>.a3m` |
| Protenix | JSON 数组 | `count: 2` | 否(强制 A/B) | `unpairedMsaPath` / `pairedMsaPath` / `templatesPath` |
| Chai-1 | FASTA(每链一条记录) | 重复记录 | — | `msa_directory` 下的 `<sha256>.aligned.pqt` |
| ColabFold | FASTA(`:` 连接)或 CSV | a3m 头部 cardinality | — | `.a3m` 输入,或 CSV 的 `a3mpath` 列 |
| AlphaFold2-Multimer | FASTA(每链一条记录) | 重复记录 | 否(从 B 开始) | `--use-precomputed-msas` |
| SeedFold | JSON(单个/批量两种) | `copies: 2` | — | **无**,网页端自算 |

## 逐项要点

### AlphaFold 3
`dialect` 必须是 `"alphafold3"`,`version` 当前支持 1–4。`id` 可以是字符串或列表,
给列表即表示同一实体的多份拷贝。给了 MSA 之后 `templates` 仍可不给或设 `null`,
AF3 会自行搜模板。自定义 MSA 的第一条序列必须与 query 完全相同,且移除所有小写
插入后每条长度必须与 query 一致。`--max_template_date` 默认 2021-09-30。

**不要放任它自建 MSA**:数据流水线按 job 运行,N×N 面板会把同一条序列搜 N 次。

### AlphaFold Server
顶层必须是列表。实体键是 `proteinChain` 而非 `protein`。`modelSeeds` 是**字符串**
数组。配体只接受 19 个 CCD 白名单代号,**不接受自定义 SMILES**;离子白名单 10 个。
`maxTemplateDate` 对训练截止后的目标很重要。学术账号每天 20–30 个 job,writer 必须
按配额分批出文件。

官方规范说 `modelSeeds` 支持多个 seed,但网页端上传实测只接受单个——writer 默认
一 job 一 seed,两种情况都能过。

### Boltz-2
`msa` 字段三种取值:不给(=自动搜索)、`"empty"`(=单序列模式)、路径(=复用)。
面板场景两个默认值都不要:前者按对重复搜索,后者丢弃比对。

**同一输入目录里不能混放带 `affinity` 和不带 `affinity` 的条目。** Boltz 会对目录里
每条都跑 affinity 阶段,然后在缺 `pre_affinity_*.npz` 的条目上崩溃。带亲和力的案例
必须单独放一个目录。

### Protenix
`msa: {precomputed_msa_dir, pairing_db}` 是**已废弃**格式,现在用 proteinChain 下的
`pairedMsaPath` / `unpairedMsaPath` / `templatesPath`,建议绝对路径。配体三种写法:
裸 SMILES、`CCD_ATP`、`FILE_/path/to.sdf`。`protenix pred -s 2066,318,1030` 一次可
给多个 seed。

### Chai-1
FASTA 头部分隔符是 `|`,实体类型 `protein` / `ligand` / `rna` / `dna` / `glycan`,
`>protein|name=xxx` 和 `>protein|xxx` 两种写法都接受。

MSA 走 `msa_directory` 下的 `<sha256>.aligned.pqt`,hash 是序列大写后的 SHA-256,
一个唯一链序列一个文件。四列:`sequence` / `source_database` / `pairing_key` /
`comment`,跨链配对靠 `pairing_key` 匹配。

**CLI `chai-lab fold` 会永久丢失 PAE**:`scores.model_idx_*.npz` 只有 ranking 汇总,
完整 PAE 只在 Python API `run_inference` 返回的 `StructureCandidates.pae` 里。

### ColabFold
多链用 `:` 连在**同一条**序列里,与 AF2-Multimer 的多条记录写法相反。

复合物 a3m 的头部格式(`colabfold/input.py` 的 `msa_to_str`):

```
#len1,len2	card1,card2
>101
<拼接后的 query 序列>
<配对块,然后是各链的非配对块>
```

第 0 行是 `#` + 逗号分隔的各链长度 + TAB + 逗号分隔的拷贝数。错一个字符就会被当作
单体处理。

`get_queries` 还接受 `.csv` / `.tsv`,列名 `id` / `sequence`,可选 `a3mpath` /
`templatepath` ——比往目录里扔 a3m 干净得多。

### AlphaFold2-Multimer
一个复合物一个文件,每条链一条记录。没有 seed 选项;
`--num-multimer-predictions-per-model N` 得到的是**同一 MSA 下的多次采样**,不是
随机种子,命名上要区分。纯 ipTM 不在 `ranking_debug.json` 里(那里只有
0.8·ipTM+0.2·pTM 的组合分),要从 `result_model_*.pkl` 读。

### SeedFold
两种输入。单个预测是裸的实体列表:

```json
[{"entity": "Protein", "copies": 4, "sequence": "QLEDSEVEAVAKGLEEM..."}]
```

批量提交是 job 列表:

```json
[{"job_name": "example_job_001", "model": "SeedFold-Linear_v1.0.0",
  "entities": [{"entity": "Protein", "copies": 4, "sequence": "..."}]}]
```

实体类型五种:`Protein` / `DNA` / `RNA` / `Ligand/Ion-CCD` / `Ligand-Smiles`,
全部用同一个 `sequence` 键装内容(CCD 代号和 SMILES 也放这里)。
`model` 可选 `SeedFold_v1.0.0` 与 `SeedFold-Linear_v1.0.0`。

**没有 MSA 字段,也没有 seed 字段。** 它是八家里唯一完全无法复用外部比对的引擎,
采样随机性也不可复现。任何跨引擎比较都要把它单独标注,否则会把"自算 MSA"的差异
误读成模型差异。


## AlphaFold2-Multimer 的 MSA 复用只能做一半

`--use-precomputed-msas` 是**逐文件**判断的(`pipeline.py`:
`if not use_precomputed_msas or not os.path.exists(msa_out_path)`),文件在就读,
不在就自己搜。它找的文件是:

| 文件 | 格式 | 我们能不能给 |
|---|---|---|
| `bfd_uniref_hits.a3m` | a3m | **能**,直接放缓存里的 unpaired |
| `uniref90_hits.sto` | Stockholm | 不给 |
| `mgnify_hits.sto` | Stockholm | 不给 |
| `uniprot_hits.sto` | Stockholm | 不给 |

只写 a3m 那一个。Stockholm 那三个理论上可以转换,但 `uniprot_hits.sto` 是多聚体
**配对**的来源,配对靠的是序列标识里的物种信息;转换过程中如果标识形态变了,配对会
静默变空——那不是复用 MSA,是悄悄换了一个 MSA。宁可让 AF2 自己去搜。

因此 AF2 在 manifest 里虽然 `msa_used=1`,但 `extra.msa_partial` 会列出仍由它自己
搜索的那三个文件。任何跨引擎比较都要把这一点算进去。
