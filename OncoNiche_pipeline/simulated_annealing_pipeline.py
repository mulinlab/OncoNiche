"""OncoNiche SA:给定 uKIN 输出的种子相关子网络,用模拟退火寻找 conductance
最小、组织特异基因富集的紧密子网络。

"""
from pathlib import Path
from multiprocessing import Pool
import argparse, functools, math, random
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import fisher_exact


## ============ 基础网络操作 ============

def load_network(path):
    return nx.read_edgelist(path)


def ukin_seed_subnet(tissue, seed, G_bg, top_k, work_dir):
    """读 uKIN 输出取 top-K 基因,返回包含 seed 的连通分量节点集。"""
    path = Path(work_dir) / f"OncoNiche_pipeline/uKIN_pipeline/output" \
                            f"/uKIN_seed_SIA_CGC_{tissue}/output_{seed}_results.txt"
    top = pd.read_csv(path, sep=" ", header=None, nrows=top_k)[0].tolist()
    for comp in nx.connected_components(G_bg.subgraph(top)):
        if seed in comp:
            return set(comp)
    return set()


def cc_with_seed(G, S, seed):
    """G 限制在节点集 S 上,取含 seed 的连通分量。"""
    if seed not in S:
        return set()
    for comp in nx.connected_components(G.subgraph(S)):
        if seed in comp:
            return set(comp)
    return set()


## ============ 评分与检验 ============

def conductance(G, S):
    """S 相对 G 的 conductance = cut / (2*internal + cut)。"""
    if not S:
        return 1.0
    cut = nx.cut_size(G, S)                 # 边界边数,networkx 底层实现
    vol = sum(dict(G.degree(S)).values())   # ∑ deg_G(v) for v in S
    denom = vol                              # = 2*internal + cut
    return cut / denom if denom > 0 else 1.0


def fisher_pvalue(S_nodes, ts_set, all_nodes):
    """2×2 表: [子网络∩ts | 非子网络∩ts; 子网络非ts | 非子网络非ts]。"""
    S = set(S_nodes)
    a = len(S & ts_set)
    b = len(ts_set - S)
    c = len(S - ts_set)
    d = len(all_nodes) - a - b - c
    _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return p


## ============ SA 主循环 ============

def cooling(n_iter, T0, alpha):
    return T0 * alpha ** np.arange(1, n_iter + 1)


def sa_optimize(
    seed, S_init, G_bg, ts_gene, all_nodes,
    n_iter=5000, T0=10, alpha=0.95, p_cutoff=0.01,
    init_frac=0.5, max_init_retry=None, rng_seed=None,
):
    """给定初始子网 S_init(节点集),模拟退火迭代。"""
    if seed not in S_init:
        raise ValueError(f"seed {seed!r} not in S_init")

    rng      = random.Random(rng_seed)
    T        = cooling(n_iter, T0, alpha)
    N_all    = list(S_init)                     # 候选池 = 初始子网所有节点
    non_seed = [n for n in N_all if n != seed]
    ts_set   = set(ts_gene) & set(all_nodes)
    max_init_retry = max_init_retry or max(1, len(N_all) // 4)

    ## ---- 初始化:随机丢 init_frac 节点,重试直到富集显著 ----
    def _try_init():
        drop = set(rng.sample(non_seed, k=int(len(N_all) * init_frac)))
        S    = cc_with_seed(G_bg, set(N_all) - drop, seed)
        return S, fisher_pvalue(S, ts_set, all_nodes)

    S, p_curr = _try_init()
    for _ in range(max_init_retry):
        if len(S) >= len(N_all) / 4 and p_curr <= p_cutoff:
            break
        S, p_curr = _try_init()

    cond_curr = conductance(G_bg, S)

    ## ---- SA 主循环 ----
    history, members = [], []
    for t in range(n_iter):
        members.append({
            "seed": seed,
            "Subnetwork member genes": "|".join(sorted(S)),
            "Times": t,
        })

        ## 提议:翻转一个非种子节点
        v     = rng.choice(non_seed)
        S_new = cc_with_seed(G_bg, (S - {v}) if v in S else (S | {v}), seed)

        if not S_new:                          # 翻转后种子被隔离 → 拒绝
            decision, cond_new, p_new, diff, metro, r = "refuse", np.nan, np.nan, np.nan, np.nan, np.nan
        else:
            cond_new = conductance(G_bg, S_new)
            p_new    = fisher_pvalue(S_new, ts_set, all_nodes)
            diff     = cond_new - cond_curr
            metro    = min(1.0, math.exp(-diff / T[t])) if diff > 0 else 1.0
            r        = rng.random()
            decision = "apply" if (r < metro and p_new < p_cutoff) else "refuse"

        history.append({
            "Times": t,
            "Conductance score": cond_new,
            "Conductance score difference": diff,
            "Probability": metro,
            "Random number": r,
            "Temperature": T[t],
            "P value": p_new,
            "Number of subnetwork members": len(S_new if decision == "apply" else S),
            "Rotation decision": decision,
        })

        if decision == "apply":
            S, cond_curr = S_new, cond_new

    return pd.DataFrame(history), pd.DataFrame(members)


## ============ 端到端 driver ============

def process_seed(seed, *, tissue, G_bg, all_nodes, ts_gene, work_dir,
                 top_k=100, **sa_kwargs):
    out_dir = Path(work_dir) / f"OncoNiche_pipeline/OncoNiche_output/{tissue}"
    out_dir.mkdir(parents=True, exist_ok=True)

    S_init = ukin_seed_subnet(tissue, seed, G_bg, top_k, work_dir)
    if not S_init or seed not in S_init:
        return seed, "skip:no_seed_in_ukin_subnet"

    hist, members = sa_optimize(seed, S_init, G_bg, ts_gene, all_nodes, **sa_kwargs)
    hist["tissue"] = tissue;  hist["ts_seed"] = seed
    members["tissue"] = tissue

    ## 命名统一 —— 与后处理 filter 函数对齐(修复原代码 bug)
    hist.to_csv(   out_dir / f"{tissue}_{seed}_subnetwork_argument.txt",     sep="\t", index=False)
    members.to_csv(out_dir / f"{tissue}_{seed}_subnetwork_member_genes.txt", sep="\t", index=False)
    return seed, "ok"


def filter_final_subnetworks(tissue, work_dir, min_size=10):
    """扫描所有 seed 的 SA 结果,取最后 apply 的子网络成员,合并输出。"""
    out_dir = Path(work_dir) / f"OncoNiche_pipeline/OncoNiche_output/{tissue}"
    all_rows = []

    for arg_file in out_dir.glob(f"{tissue}_*_subnetwork_argument.txt"):
        seed        = arg_file.name.replace(f"{tissue}_", "").replace("_subnetwork_argument.txt", "")
        member_file = out_dir / f"{tissue}_{seed}_subnetwork_member_genes.txt"
        if not member_file.exists():
            continue

        hist    = pd.read_csv(arg_file, sep="\t")
        applied = hist[hist["Rotation decision"] == "apply"]
        if applied.empty:
            continue
        last_t = applied["Times"].iloc[-1]

        members = pd.read_csv(member_file, sep="\t")
        row     = members.loc[members["Times"] == last_t].iloc[0]
        genes   = row["Subnetwork member genes"].split("|")
        if len(genes) < min_size:
            continue

        df = pd.DataFrame({
            "Subnetwork member genes":       genes,
            "Tissue-specific mutated genes": seed,
            "P value":                       applied["P value"].iloc[-1],
            "Conductance score":             applied["Conductance score"].iloc[-1],
            "Number of subnetwork members":  applied["Number of subnetwork members"].iloc[-1],
            "Tissue":                        tissue,
        })
        all_rows.append(df)

    if not all_rows:
        return pd.DataFrame()
    result = pd.concat(all_rows, ignore_index=True).drop_duplicates()
    result.to_csv(out_dir / f"All_tissue_specific_subnetworks_in_{tissue}.txt",
                  sep="\t", index=False)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue",     required=True)
    ap.add_argument("--work_dir",   required=True, type=Path)
    ap.add_argument("--n_proc",     type=int, default=1)
    ap.add_argument("--n_iter",     type=int, default=5000)
    ap.add_argument("--top_k_ukin", type=int, default=100)
    args = ap.parse_args()

    G_bg = load_network(
        args.work_dir /
        "OncoNiche_pipeline/uKIN_pipeline/Global_Gini_CGC_driver/background_network/Nested Systems in Tumors network.tsv"
    )
    all_nodes = set(G_bg.nodes())

    base    = args.work_dir / f"OncoNiche_pipeline/uKIN_pipeline/Global_Gini_CGC_driver/{args.tissue}_global_Gini"
    ts_mut  = pd.read_csv(base / "tissue_mut_score", sep="\t", header=None)[1].tolist()
    ts_exp  = pd.read_csv(base / "tissue_exp_score", sep="\t", header=None)[1].tolist()
    ts_gene = set(ts_mut) | set(ts_exp)

    ukin_out = args.work_dir / f"OncoNiche_pipeline/uKIN_pipeline/output/uKIN_seed_SIA_CGC_{args.tissue}"
    seeds    = [f.name.split("output_")[1].split("_results.txt")[0]
                for f in ukin_out.glob("output_*_results.txt")]
    print(f"[{args.tissue}] {len(seeds)} seeds")

    fn = functools.partial(
        process_seed,
        tissue    = args.tissue,
        G_bg      = G_bg,
        all_nodes = all_nodes,
        ts_gene   = ts_gene,
        work_dir  = args.work_dir,
        top_k     = args.top_k_ukin,
        n_iter    = args.n_iter,
    )

    if args.n_proc > 1:
        with Pool(args.n_proc) as pool:
            for seed, status in pool.imap_unordered(fn, seeds):
                print(f"[{status}] {seed}", flush=True)
    else:
        for seed in seeds:
            _, status = fn(seed); print(f"[{status}] {seed}", flush=True)

    filter_final_subnetworks(args.tissue, args.work_dir)


if __name__ == "__main__":
    main()
