# 各引擎的输出布局

写给下游读取方(主要是 `foldmetrics`)。**每一条都标注了来源**:

- **实测** —— 在 A100/A40 上由 `foldrunner run` 真实跑出、我逐个 `find` 看过的路径
- **未验证** —— 尚未真跑,不要据此实现

不要把"未验证"当事实用。这份文档里凡是没标"实测"的,都可能是错的。

## foldrunner 的外层约定

`foldrunner run -o <results>` 建立的外层目录,是我们能控制的唯一一层:

```
<results>/<engine>/<stem>/    ← 到这里为止由 foldrunner 决定
```

`<stem>` 取决于引擎的**运行粒度**:

| 粒度 | `<stem>` | 引擎 |
|---|---|---|
| 每个复合物一次 | job 名,如 `P0__P1` | af3、chai1、af2_multimer |
| 整个面板一次 | `panel` | boltz2、protenix、colabfold |
| 多 seed 拆分 | `<job>_seed<N>` 或 `panel_seed<N>` | 命令模板含单个 `{seed}` 的引擎 |

**这一层之下的结构完全由引擎自己决定,foldrunner 不干预。** 下游能否识别 target,
取决于引擎自己建不建 job 目录——这正是差异的来源。

## 逐引擎

### Boltz-2 —— 实测

```
<results>/boltz2/panel_seed2066/          多 seed 拆成平级目录
└── boltz_results_boltz2/
    ├── lightning_logs/version_0/hparams.yaml
    ├── msa/
    ├── processed/records/P0__P0.json
    └── predictions/P0__P0/
        ├── P0__P0_model_0.cif
        ├── confidence_P0__P0_model_0.json
        ├── pae_P0__P0_model_0.npz
        ├── pde_P0__P0_model_0.npz
        └── plddt_P0__P0_model_0.npz
<results>/boltz2/panel_seed318/           同上
```

132 个文件里只有 12 个是结构。**输出文件名里不带 seed**,所以 foldrunner 必须把
不同 seed 拆进不同目录,否则第二次运行会覆盖第一次。

Boltz 首次运行会下载 CCD 数据(`mols.tar`)。离线节点上必须用 `--cache` 指向已有
的本地副本,否则卡在网络超时。

面板级运行,但 Boltz 自己会建 `predictions/<job>/`,所以 target 可从目录取得。

### Protenix —— 实测

```
<results>/protenix/panel/
└── P0__P1/
    └── seed_2066/
        └── predictions/
            ├── P0__P1_sample_0.cif
            ├── P0__P1_summary_confidence_sample_0.json
            └── P0__P1_full_data_sample_0.json     ← 需要 --need_atom_confidence
```

面板级运行,job 目录在最外层,中间隔着 `seed_<N>/predictions/`。

**`full_data` 里才有 `token_pair_pae`。** 不加 `--need_atom_confidence True` 就只有
`.cif` 和 `summary_confidence`,所有 PAE 类指标(ipSAE / pDockQ2 / LIS)全部拿不到,
且事后无法补救。foldrunner 已把这个 flag 写进默认命令。

### ColabFold —— 实测

```
<results>/colabfold/panel/
├── P0__P0_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb
├── P0__P0_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json
├── P0__P0_predicted_aligned_error_v1.json
├── P0__P0.a3m                    输入的副本
├── P0__P0_coverage.png
├── P0__P0_pae.png
├── P0__P0_plddt.png
├── P0__P0.done.txt
├── config.json / log.txt / cite.bibtex
└── ...（其余复合物同样平铺）
```

51 个文件里只有 6 个是结构。

**平铺,不建 job 目录。** 这是唯一一个下游无法从目录路径取得 target 的引擎——
job 名只存在于文件名前缀里。

命名规则(由 ColabFold 自己保证,稳定):

```
{jobname}_(un)relaxed_rank_{NNN}_{model_type}_model_{N}_seed_{NNN}.(pdb|cif)
{jobname}_scores_rank_{NNN}_{model_type}_model_{N}_seed_{NNN}.json
```

`{jobname}` 经过 ColabFold 的 `safe_filename()`(非 `[A-Za-z0-9_.-]` 变 `_`),
foldrunner 的 `A__B` 形式原样保留。

改成每复合物一次运行可以解决,但 ColabFold 每次调用都要重载 AlphaFold2 参数,
面板一大就全耗在加载上,所以维持面板级运行。

### AlphaFold 3 —— 实测

```
<results>/af3/P0__P0/          ← foldrunner 建的
└── P0__P0/                    ← AF3 建的，用 JSON 里的 name，无时间戳
    ├── P0__P0_model.cif                 排名最佳
    ├── P0__P0_confidences.json
    ├── P0__P0_summary_confidences.json
    ├── P0__P0_ranking_scores.csv
    ├── P0__P0_data.json                 解析后的输入，不是结果
    ├── seed-2066_sample-{0..4}/
    │   ├── P0__P0_seed-2066_sample-0_model.cif
    │   ├── ..._confidences.json
    │   └── ..._summary_confidences.json
    └── seed-318_sample-{0..4}/          多 seed 在同一目录下并列
```

一个复合物 216 个文件里 11 个 cif,但**只有 10 个是独立模型**:顶层
`P0__P0_model.cif` 是排名第一那个 sample 的副本,递归扫会重复计一次,
`ranking_scores.csv` 能对出是哪一个。

**目录名不带时间戳。** 时间戳形式(`P0__P0_20260818_145715`)只在同一输出目录
重复运行时出现——AF3 不覆盖,第二次起加时间戳后缀。拿重跑结果当原生布局会得出
错误结论。

容器化运行两个必须注意的点:

- job JSON 里的 `unpairedMsaPath` / `pairedMsaPath` 是**绝对路径**,MSA 缓存必须
  以相同路径挂进容器,否则 `--norun_data_pipeline` 找不到比对
- 给了 MSA 就等于跳过 data pipeline,而跳过的 pipeline 也填不了 templates。
  **`templates` 必须显式给出**(空列表即可),否则模型加载完才报
  `Protein chain N is missing Templates`

### Chai-1 —— 实测

```
<results>/chai1/P0__P0/
├── pred.model_idx_{0..4}.cif
├── scores.model_idx_{0..4}.npz
└── pae.npz                     ← foldrunner 的脚本额外保存，形状 (5, N, N)
```

foldrunner 通过生成的 `run_chai.py` 驱动而非 `chai-lab fold`:CLI 只写 ranking
汇总,完整 PAE 只存在于 Python API 返回对象上,CLI 跑完就永久丢失。

**比对不由 foldrunner 转成 parquet。** parquet 只能向后兼容:新版 pyarrow 写的
文件在旧版上会以 `Repetition level histogram size mismatch` 失败,而读它的是引擎
环境,版本不受本包控制(实测写 25.0.1 / 读 19.0.0 即失败,五种写法参数全无效)。
所以 foldrunner 只写 a3m,由跑在 Chai 环境里的脚本用 Chai 自己的函数转换。

**每条记录的 name 必须唯一。** 同源二聚体不能写两条同名记录,foldrunner 用链名
后缀区分(`P0_A` / `P0_B`)。注意 `read_inputs()` 会接受重名文件,拒绝发生在更
后面的特征构建阶段。

### AlphaFold2-Multimer —— 实测

```
<results>/af2_multimer/P0__P0/     ← foldrunner 建的
└── P0__P0/                        ← AF2 按 FASTA basename 建的
    ├── ranked_{0..4}.pdb + .cif                   按分数排序
    ├── unrelaxed_model_{1..5}_multimer_v3_pred_0.pdb + .cif
    ├── confidence_model_{1..5}_..._pred_0.json
    ├── pae_model_{1..5}_..._pred_0.json
    ├── result_model_{1..5}_..._pred_0.pkl         大文件
    ├── ranking_debug.json                         键: iptm+ptm, order
    ├── features.pkl / timings.json / .done
    └── msas/
        ├── A/{bfd_uniref_hits.a3m, mgnify_hits.sto, pdb_hits.sto, uniref90_hits.sto}
        └── chain_id_map.json
```

20 个结构文件里**只有 5 个独立模型**:每个模型有 `.pdb` 和 `.cif` 两份,
`ranked_*` 又是同样五个的副本(按分数排序)。

**两套链命名,不要混。**

```json
chain_id_map.json:  "A" → description "P0__P0_B"
                    "B" → description "P0__P0_C"
```

- `msas/<X>/` 从 **A** 开始,按**唯一序列**编号——一个序列搜一次,
  所以同源二聚体只有一个目录,不是两个
- 预测结构里的链从 **B** 开始
- `chain_id_map.json` 是连接两者的唯一凭据

**纯 ipTM 不在 `ranking_debug.json` 里**,那里只有 `iptm+ptm` 组合分,
要从 `result_model_*.pkl` 读。

运行时必须显式要求多聚体预设。AF2 自己和常见 wrapper 都默认单体流水线,
而单体流水线会以 `More than one input sequence found` 拒绝多条记录的 FASTA
——报错发生在任务已经启动之后。参数拼写按入口不同:`run_alphafold.py` 用
`--model_preset=multimer`。

### AlphaFold Server / SeedFold —— 不适用

网页表单,没有命令行。结果由人工下载,目录结构取决于用户怎么放。
`scripts/MANUAL.md` 会列出要上传的文件。

## 跨引擎对齐靠什么

job 名是唯一贯穿所有引擎的键,形式 `<A>__<B>`,双下划线分隔,实体名清洗到
`[A-Za-z0-9._-]`。同一个 job 在八家里名字完全一致。

`manifest.tsv` 记录了每家实际分配的链名,因为 Protenix 和 AlphaFold Server 会把
链名强制改成 A/B/C,AF2-Multimer 从 B 开始——不看 manifest 就对不上号。
