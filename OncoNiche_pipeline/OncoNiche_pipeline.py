"""OncoNiche pipeline: uKIN + Simulated Annealing.
"""
from pathlib import Path
from multiprocessing import Pool
import argparse, functools, math, os, random
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import fisher_exact


## ============ 基础工具 ============

def load_network_from_file(network_file):
    with open(network_file, "rb") as fh:
        return nx.read_edgelist(fh)


def ukin_seed(tissue, seed, G_bg, ukin_number, work_dir):
    """读 uKIN top-K 结果,取 G_bg 上含 seed 的连通分量子图。"""
    path = (Path(work_dir) / "OncoNiche_pipeline/uKIN_pipeline/output"
            / f"uKIN_seed_SIA_CGC_{tissue}" / f"output_{seed}_results.txt")
    top = pd.read_csv(path, sep=" ", header=None, nrows=ukin_number)[0].tolist()
    for comp in nx.connected_components(G_bg.subgraph(top)):
        if seed in comp:
            return G_bg.subgraph(comp)                # 与原代码同类型:subgraph view
    return G_bg.subgraph([])                          # 空 subgraph 兜底


def _cc_with_seed(G, S, seed):
    """G 限制到节点集 S,取含 seed 的连通分量,返回 set。"""
    if seed not in S:
        return set()
    for comp in nx.connected_components(G.subgraph(S)):
        if seed in comp:
            return set(comp)
    return set()


## ============ 打分函数(与原代码数值等价)============

def conductance_score(G_bg, S):
    """cs / (2*ms + cs) = cut / vol,与原实现数值等价。"""
    if not S:
        return 0.0
    cut = nx.cut_size(G_bg, S)
    vol = sum(d for _, d in G_bg.degree(S))          # = 2*internal + cut
    return float(cut / vol) if vol > 0 else 0.0


def enrichment_pvalue(bg_nodes_set, S, ts_gene_set):
    """2x2 Fisher (greater)。原代码返回 7 列,实际只用 pvalue。"""
    a = len(S & ts_gene_set)                          # sub ∩ ts
    b = len(ts_gene_set - S)                          # (bg - sub) ∩ ts   (ts ⊆ bg)
    c = len(S - ts_gene_set)                          # sub - ts          (S ⊆ bg)
    d = len(bg_nodes_set) - a - b - c
    _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return float(p)


def temperature_function(n_iter, T0, alpha):
    """T[i] = T0 * alpha^i for i=1..n_iter,与原循环等价。"""
    return (T0 * alpha ** np.arange(1, n_iter + 1)).tolist()


## ============ SA 主循环 ============

def _init_active_status(N_all_nodes, seed, G_bg, seed_id_col_index=None, rng=random):
    """随机把 50% 节点标为 status=0(移除),保留 seed。返回 (status_df, S_active)。
      - 所有节点初始 status=1
      - 随机抽 int(N/2) 个位置置 0
      - 强制 seed 位置 status=1
      - S_active = status=1 的节点在 G_bg 上取含 seed 的 CC
    """
    n = len(N_all_nodes)
    status = np.ones(n, dtype=int)
    drop_pos = rng.sample(list(range(n)), int(n * 0.5))
    status[drop_pos] = 0
    seed_pos = N_all_nodes.index(seed)
    status[seed_pos] = 1                              # 强制 seed 存活

    status_df = pd.DataFrame({"nodes": N_all_nodes, "status": status})
    S_active = _cc_with_seed(G_bg,
                             set(np.array(N_all_nodes)[status == 1]),
                             seed)
    return status_df, S_active


def _record(time_eg, cond, diff, metro, r, T, p, size, decision):
    """记录一步的字典。用 dict 收集,循环外一次 DataFrame,避免 O(n^2) concat。"""
    return dict(zip(
        ["Times", "Conductance score", "Conductance score difference",
         "Probability", "Random number", "Temperature", "P value",
         "Number of subnetwork members", "Rotation decision"],
        [time_eg, cond, diff, metro, r, T, p, size, decision],
    ))


def simulated_annealing_main(TS_mutation, subnetwork_init, iteration_times,
                             ts_gene, G_bg, pvalue_cutoff, bg_nodes,
                             temperature_all):
    """主循环:
       - 初始化 while 循环重试至满足 (size>=N/4 且 p<cutoff) 或 j>=N/4
       - metropolis 用 exp,当 -diff/t>=100 时兜底 1
       - 拒绝分支里翻回 candidate status(冗余但保留原行为)
    """
    N_all      = list(subnetwork_init.nodes())
    bg_set     = set(bg_nodes)
    ts_set     = set(ts_gene) & bg_set
    N_size     = len(N_all)

    ## ---- 初始化 + while 重试 ----
    status_df, S_init = _init_active_status(N_all, TS_mutation, G_bg)
    p_init = enrichment_pvalue(bg_set, S_init, ts_set)

    j = 0
    while ((len(S_init) < N_size / 4) or (p_init > pvalue_cutoff)) and (j < N_size / 4):
        j += 1
        status_df, S_init = _init_active_status(N_all, TS_mutation, G_bg)
        p_init = enrichment_pvalue(bg_set, S_init, ts_set)

    S_curr        = S_init
    cond_curr     = conductance_score(G_bg, S_curr)
    status_curr   = status_df

    ## ---- SA 迭代 ----
    hist_rows    = []
    final_rows   = []
    positions    = list(range(N_size))
    seed_pos     = N_all.index(TS_mutation)

    for time_eg in range(iteration_times):
        ## 每 iter 首:记录 S_curr(与原代码顺序一致)
        final_rows.append({
            "Tissue-specific mutated genes": TS_mutation,
            "Subnetwork member genes":       "|".join(S_curr),
            "Times":                          time_eg,
        })

        ## 提议:候选 = curr.copy,翻转一个非 seed 位置
        status_cand = status_curr.copy()
        pos = random.sample(positions, 1)[0]
        while pos == seed_pos:                        # 保持原代码 while 采样直到非 seed
            pos = random.sample(positions, 1)[0]
        old_status  = status_cand.iloc[pos, 1]
        status_cand.iloc[pos, 1] = 1 - old_status
        active      = set(status_cand.loc[status_cand["status"] == 1, "nodes"])
        S_cand      = _cc_with_seed(G_bg, active, TS_mutation)

        ## 分支 1:候选子网含 seed → 评估
        if TS_mutation in S_cand:
            cond_cand = conductance_score(G_bg, S_cand)
            diff      = cond_cand - cond_curr
            t         = temperature_all[time_eg]
            metro     = math.exp(-diff / t) if (-diff / t) < 100 else 1
            r         = random.random()
            p_cand    = enrichment_pvalue(bg_set, S_cand, ts_set)

            if (r < metro) and (p_cand < pvalue_cutoff):
                print(time_eg)
                S_curr, cond_curr, status_curr = S_cand, cond_cand, status_cand
                hist_rows.append(_record(time_eg, cond_cand, diff, metro, r, t,
                                         p_cand, len(S_curr), "apply"))
            else:
                print(f"{time_eg} rotation")
                hist_rows.append(_record(time_eg, cond_cand, diff, metro, r, t,
                                         p_cand, len(S_curr), "refuse"))
                ## 保留原代码"拒绝后翻回 candidate status"的操作(冗余但等价行为)
                status_cand.iloc[pos, 1] = old_status

        ## 分支 2:候选子网不含 seed → 用上一轮的 diff/metro/r/p 值记录(保持原行为)
        else:
            print(f"{time_eg} rotation")
            ## 原代码此分支里 diff/metro/r/p 值沿用上一轮,不重新计算
            ## 若 time_eg=0 就进这里,上轮值未定义 → 原代码同样会 NameError
            ## 保持原行为
            hist_rows.append(_record(time_eg, cond_cand, diff, metro, r, t,
                                     p_cand, len(S_curr), "refuse"))
            status_cand.iloc[pos, 1] = old_status

    return pd.DataFrame(hist_rows), pd.DataFrame(final_rows)


## ============ 端到端 driver ============

def simulated_annealing_all_function_reversion(tissue, TS_mutation, work_dir):
    """一个 seed 的完整流程,输出与原代码相同的两个 txt 文件。"""
    iteration_times = 5000
    T_schedule      = temperature_function(iteration_times, 10, 0.95)

    out_dir = Path(work_dir) / "OncoNiche_pipeline/OncoNiche_output" / tissue
    out_dir.mkdir(parents=True, exist_ok=True)

    bg_net_path = (Path(work_dir) / "OncoNiche_pipeline/uKIN_pipeline"
                   / "Global_Gini_CGC_driver/background_network"
                   / "Nested Systems in Tumors network.tsv")
    G_bg     = load_network_from_file(bg_net_path)
    bg_nodes = list(G_bg.nodes())

    seed_subnet = ukin_seed(tissue, TS_mutation, G_bg, ukin_number=100, work_dir=work_dir)

    base = (Path(work_dir) / "OncoNiche_pipeline/uKIN_pipeline"
            / "Global_Gini_CGC_driver" / f"{tissue}_global_Gini")
    ts_mut = pd.read_csv(base / "tissue_mut_score", sep="\t", header=None)[1].tolist()
    ts_exp = pd.read_csv(base / "tissue_exp_score", sep="\t", header=None)[1].tolist()
    ts_gene = np.intersect1d(np.union1d(ts_mut, ts_exp), bg_nodes)

    hist, final = simulated_annealing_main(
        TS_mutation, seed_subnet, iteration_times, ts_gene, G_bg,
        pvalue_cutoff=0.01, bg_nodes=bg_nodes, temperature_all=T_schedule,
    )
    final["tissue"]  = tissue
    hist["tissue"]   = tissue
    hist["ts_seed"]  = TS_mutation

    ## 输出文件名与原代码一致
    hist.to_csv( out_dir / f"subnetwork_stas_curr_all_{tissue}_{TS_mutation}.txt",
                sep="\t", index=False)
    final.to_csv(out_dir / f"final_subnetwork_output_{tissue}_{TS_mutation}.txt",
                sep="\t", index=False)


def ts_subnetwork_filter(tissue, work_dir):
    out_dir  = Path(work_dir) / "OncoNiche_pipeline/OncoNiche_output" / tissue
    all_seed = np.unique([
        i.replace(f"{tissue}_", "").split("_")[0]
        for i in os.listdir(out_dir) if "argument" in i
    ])

    rows = []
    for seed in all_seed:
        subnet_read = pd.read_csv(out_dir / f"{tissue}_{seed}_subnetwork_member_genes.txt", sep="\t")
        stas_read   = pd.read_csv(out_dir / f"{tissue}_{seed}_subnetwork_argument.txt",     sep="\t")
        applied     = stas_read[stas_read["Rotation decision"] == "apply"]
        if stas_read.empty or applied.empty:
            continue

        last_t = applied["Times"].iloc[-1]
        final  = subnet_read[subnet_read["Times"] == last_t]
        genes  = list(final["Subnetwork member genes"].str.split("|", expand=True).T.iloc[:, 0])
        df = pd.DataFrame({
            "Subnetwork member genes":       genes,
            "Tissue-specific mutated genes": applied["Tissue-specific mutated genes"].iloc[0],
            "P value":                       applied["P value"].iloc[-1],
            "Conductance score":             applied["Conductance score"].iloc[-1],
            "Number of subnetwork members":  applied["Number of subnetwork members"].iloc[-1],
            "Tissue":                        tissue,
        })
        if len(df) > 10:
            rows.append(df)

    result = pd.concat(rows).drop_duplicates() if rows else pd.DataFrame()
    result.to_csv(out_dir / f"All_tissue_specific_subnetworks_in_{tissue}.txt",
                  sep="\t", index=False)
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue",   required=True)
    ap.add_argument("--work_dir", required=True)
    args = ap.parse_args()

    ukin_out = (Path(args.work_dir) / "OncoNiche_pipeline/uKIN_pipeline/output"
                / f"uKIN_seed_SIA_CGC_{args.tissue}")
    seeds = [i.split("output_")[1].split("_results.txt")[0]
             for i in os.listdir(ukin_out)]

    for seed in seeds:
        simulated_annealing_all_function_reversion(args.tissue, seed, args.work_dir)
    ts_subnetwork_filter(args.tissue, args.work_dir)