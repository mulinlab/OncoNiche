from multiprocessing import Pool
import functools
import os
import argparse
import subprocess
def run_seed_ukin(gene,tissue,work_dir):
    ruby_cmd=f'''ruby uKIN-master/uKIN.rb \
    /h/tianyi/TS_datasets_reversion/sci_paper_plot_stas/github/Nested Systems in Tumors network.tsv \
    {work_dir}/uKIN_pipeline/Global_Gini_CGC_driver/{tissue}_global_Gini/prior_knowledge.txt \
    {work_dir}/github/uKIN_pipeline/Global_Gini_CGC_driver/{tissue}_global_Gini/{gene}_mut_seed.txt  \
    matlab=matlab_run/bin \
    output_prefix=output_{gene} '''
    process_completed = subprocess.run(
    ruby_cmd, shell =True, encoding='utf-8', stdout = subprocess.PIPE,
    stderr = subprocess.PIPE
    )
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", help="Please input tissue type", type=str)
    parser.add_argument("--work_dir", help="Please the working directory path where the code is located", type=str)
    args = parser.parse_args()
    tissue=args.tissue
    work_dir=args.work_dir
    run_seed_ukin_partial=functools.partial(run_seed_ukin,tissue=tissue,work_dir=work_dir)
    ts_mutated_genes=[i.split('_mut_seed.txt')[0] for i in os.listdir(f'{work_dir}/uKIN_pipeline/Global_Gini_CGC_driver/{tissue}_global_Gini') if '_mut_seed.txt' in i]
    pool = Pool(10)
    run_cmd =pool.map(run_seed_ukin_partial, ts_mutated_genes) 
    pool.close()
    pool.join()