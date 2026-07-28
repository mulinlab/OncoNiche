"""uKIN wrapper:对每个 seed 并行调用 Ruby 版 uKIN。"""
from pathlib import Path
from multiprocessing import Pool
import subprocess, functools, argparse


def run_ukin(seed, *, network_path, prior_path, seed_dir, matlab_bin, output_dir, timeout=3600):
    """一个 seed 一次调用。list argv,不走 shell,避免路径含空格问题。"""
    cmd = [
        "ruby", "uKIN-master/uKIN.rb",
        str(network_path),
        str(prior_path),
        str(seed_dir / f"{seed}_mut_seed.txt"),
        f"matlab={matlab_bin}",
        f"output_prefix={output_dir / f'output_{seed}'}",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
        return seed, "ok"
    except subprocess.CalledProcessError as e:
        return seed, f"fail:{e.returncode}:{e.stderr[:200]}"
    except subprocess.TimeoutExpired:
        return seed, "timeout"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tissue",   required=True)
    p.add_argument("--work_dir", required=True, type=Path)
    p.add_argument("--n_proc",   type=int, default=10)
    args = p.parse_args()

    base       = args.work_dir / "uKIN_pipeline/Global_Gini_CGC_driver"
    seed_dir   = base / f"{args.tissue}_global_Gini"
    output_dir = args.work_dir / "OncoNiche_pipeline/uKIN_pipeline/output" / f"uKIN_seed_SIA_CGC_{args.tissue}"
    output_dir.mkdir(parents=True, exist_ok=True)

    seeds = [f.name.split("_mut_seed.txt")[0]
             for f in seed_dir.glob("*_mut_seed.txt")]
    print(f"[{args.tissue}] {len(seeds)} seeds to process")

    fn = functools.partial(
        run_ukin,
        network_path = args.work_dir / "sci_paper_plot_stas/github/Nested Systems in Tumors network.tsv",
        prior_path   = seed_dir / "prior_knowledge.txt",
        seed_dir     = seed_dir,
        matlab_bin   = "matlab_run/bin",
        output_dir   = output_dir,
    )
    with Pool(args.n_proc) as pool:
        for seed, status in pool.imap_unordered(fn, seeds):
            print(f"[{status}] {seed}", flush=True)


if __name__ == "__main__":
    main()
