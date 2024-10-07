import networkx as nx
import numpy as np
import argparse
import pandas as pd
import random
import math
from multiprocessing import Pool
import scipy
import os
import argparse


def ukin_seed(tissue,seed,background_network,ukin_number,work_dir):
    ukin_read=pd.read_csv(f'{work_dir}/OncoNiche_pipeline/uKIN_pipeline/output/uKIN_seed_SIA_CGC_{tissue}/output_{seed}_results.txt',sep=' ',header=None)
    ukin_read_top=ukin_read.iloc[0:ukin_number,:]
    subnetwork_filter=background_network.subgraph(ukin_read_top[0].tolist())
    subnetwork_filter_component= list(nx.connected_components(subnetwork_filter))
    subnetwork_filter_component_seed_nodes=[]
    for subnetwork_remove_componet_nodes_eg in subnetwork_filter_component:
        if seed in subnetwork_remove_componet_nodes_eg:
            subnetwork_filter_component_seed_nodes=subnetwork_filter_component_seed_nodes+list(subnetwork_remove_componet_nodes_eg)
    subnetwork_filter_output=background_network.subgraph(subnetwork_filter_component_seed_nodes)
    return(subnetwork_filter_output)
def load_network_from_file(network_file):
    fh = open(network_file, "rb")
    G = nx.read_edgelist(fh)
    fh.close()
    return G
def conductance_score_function(subnetwork,background_network):
    subnetwork_nodes=subnetwork.nodes()
    conn_edges=[]
    for sub_nodes_eg in subnetwork_nodes:
        sub_nodes_eg_edge=nx.edges(background_network, [sub_nodes_eg])
        conn_edges=conn_edges+list(sub_nodes_eg_edge)
    ms=len(subnetwork.edges())
    cs=len(conn_edges) - 2*ms
    return(float(cs/(2*ms+cs)))
def remove_network(subnetwork,remove_node,TS_mutation):
    subnetwork.remove_nodes_from(remove_node)
    subnetwork_remove_componet_nodes= list(nx.connected_components(subnetwork))
    subnetwork_remove_componet_nodes_contain_ts_mut=[]
    for subnetwork_remove_componet_nodes_eg in subnetwork_remove_componet_nodes:
        if TS_mutation in subnetwork_remove_componet_nodes_eg:
            subnetwork_remove_componet_nodes_contain_ts_mut=subnetwork_remove_componet_nodes_contain_ts_mut+list(subnetwork_remove_componet_nodes_eg)
    remove_subnetwork_output=subnetwork.subgraph(subnetwork_remove_componet_nodes_contain_ts_mut)
    return(remove_subnetwork_output)
def enrichment_analysis_tsgenes(background_network_nodes,subnetwork,ts_gene):
    subnetwork_nodes=subnetwork.nodes()
    table_sub_ts=np.intersect1d(subnetwork_nodes,np.intersect1d(background_network_nodes,ts_gene))
    table_nosub_ts=np.intersect1d(np.setdiff1d(background_network_nodes,subnetwork_nodes),ts_gene)
    table_sub_nots=np.intersect1d(np.setdiff1d(background_network_nodes,ts_gene),subnetwork_nodes)
    table_nosub_nots=np.intersect1d(np.setdiff1d(background_network_nodes,subnetwork_nodes),np.setdiff1d(background_network_nodes,ts_gene))
    oddsratio, pvalue = scipy.stats.fisher_exact([[len(table_sub_ts),len(table_nosub_ts)],[len(table_sub_nots),len(table_nosub_nots)]],alternative='greater')
    fisher_output=pd.DataFrame()
    fisher_output['table_sub_ts']=[len(table_sub_ts)]
    fisher_output['table_nosub_ts']=[len(table_nosub_ts)]
    fisher_output['table_sub_nots']=[len(table_sub_nots)]
    fisher_output['table_nosub_nots']=[len(table_nosub_nots)]
    fisher_output['overlap_gene']=['|'.join(table_sub_ts)]
    fisher_output['oddsratio']=[oddsratio]
    fisher_output['pvalue']=[pvalue]
    return(fisher_output)
def temperatureFunction(iteration_times,initial_temp,alpha):
    temp_iteration=[]
    for i in np.arange(1,iteration_times+1):
        initial_temp=initial_temp*alpha
        temp_iteration=temp_iteration+[initial_temp]
    return(temp_iteration)
def simulated_annealing_main_pvalue_score_rotation(TS_mutation,subnetwork,iteration_times,ts_gene,background_network,pvalue_cutoff,background_network_nodes,temperatureFunction_num):
    subnetwork_original=nx.Graph(subnetwork)
    background_network=nx.Graph(background_network)
    subnetwork_nodes=subnetwork_original.nodes()
    subnetwork_nodes_status=[1]*len(subnetwork_nodes)
    subnetwork_nodes_status_df=pd.DataFrame()
    subnetwork_nodes_status_df['nodes']=list(subnetwork_nodes)
    subnetwork_nodes_status_df['status']=subnetwork_nodes_status
    subnetwork_nodes_status_df.reset_index(drop=True,inplace=True)
    position_all=list(np.arange(0,len(subnetwork_nodes)))
    initation_pos=random.sample(list(np.arange(0,subnetwork_nodes_status_df.shape[0])),int(len(list(subnetwork.nodes()))*0.5))
    subnetwork_nodes_status_df.iloc[initation_pos,1]=0
    subnetwork_nodes_status_df.loc[subnetwork_nodes_status_df['nodes']==TS_mutation,'status']=1
    remove_nodes_initialize=subnetwork_nodes_status_df.loc[subnetwork_nodes_status_df['status']==0,'nodes'].tolist()
    subnetwork=nx.Graph(subnetwork_original)
    subnetwork_initialize=nx.Graph(remove_network(subnetwork,remove_nodes_initialize,TS_mutation))
    enrichment_testing_subnetwork_initialize=enrichment_analysis_tsgenes(background_network_nodes,subnetwork_initialize,ts_gene)
    enrichment_testing_subnetwork_initialize_pvalue=float(enrichment_testing_subnetwork_initialize['pvalue'])
    j=0
    while ((len(subnetwork_initialize.nodes()) <len(subnetwork_original)/4)  | (enrichment_testing_subnetwork_initialize_pvalue > pvalue_cutoff)) & (j<len(subnetwork_original)/4):
        j=j+1
        subnetwork_nodes=subnetwork_original.nodes()
        subnetwork_nodes_status=[1]*len(subnetwork_nodes)
        subnetwork_nodes_status_df=pd.DataFrame()
        subnetwork_nodes_status_df['nodes']=subnetwork_nodes
        subnetwork_nodes_status_df['status']=subnetwork_nodes_status
        subnetwork_nodes_status_df.reset_index(drop=True,inplace=True)
        position_all=list(np.arange(0,len(subnetwork_nodes)))
        initation_pos=random.sample(list(np.arange(0,subnetwork_nodes_status_df.shape[0])),int(len(list(subnetwork.nodes()))*0.5))
        subnetwork_nodes_status_df.iloc[initation_pos,1]=0
        subnetwork_nodes_status_df.loc[subnetwork_nodes_status_df['nodes']==TS_mutation,'status']=1
        remove_nodes_initialize=subnetwork_nodes_status_df.loc[subnetwork_nodes_status_df['status']==0,'nodes'].tolist()
        subnetwork_initialize=nx.Graph(remove_network(subnetwork,remove_nodes_initialize,TS_mutation))
        subnetwork=nx.Graph(subnetwork_original)
        enrichment_testing_subnetwork_initialize=enrichment_analysis_tsgenes(background_network_nodes,subnetwork_initialize,ts_gene)
        enrichment_testing_subnetwork_initialize_pvalue=float(enrichment_testing_subnetwork_initialize['pvalue'])
    conductance_score_curr,subnetwork_curr,subnetwork_nodes_status_df_curr=conductance_score_function(subnetwork_initialize,background_network),nx.Graph(subnetwork_initialize),subnetwork_nodes_status_df
    subnetwork_stas_curr_all=pd.DataFrame()
    final_subnetwork_output_all=pd.DataFrame()
    temperature_all=temperatureFunction_num
    for time_eg in np.arange(0,iteration_times):
        final_subnetwork_output=pd.DataFrame()
        final_subnetwork_output['Tissue-specific mutated genes']=[TS_mutation]
        final_subnetwork_output['Subnetwork member genes']='|'.join(subnetwork_curr.nodes())
        final_subnetwork_output['Times']=time_eg
        final_subnetwork_output_all=pd.concat([final_subnetwork_output_all,final_subnetwork_output])
        subnetwork_nodes_status_df_candidate=subnetwork_nodes_status_df_curr.copy()
        random_position=random.sample(position_all,1)
        while subnetwork_nodes_status_df_candidate.iloc[random_position[0],:]['nodes']==TS_mutation:
            random_position=random.sample(position_all,1)
        candidate_status=subnetwork_nodes_status_df_candidate.iloc[random_position[0],:]['status']
        if candidate_status==0:
            subnetwork_nodes_status_df_candidate.iloc[random_position[0],1]=1
        else:
            subnetwork_nodes_status_df_candidate.iloc[random_position[0],1]=0
        remove_nodes=subnetwork_nodes_status_df_candidate.loc[subnetwork_nodes_status_df_candidate['status']==0,'nodes'].tolist()
        subnetwork_candidate=nx.Graph(remove_network(subnetwork,remove_nodes,TS_mutation))
        subnetwork=nx.Graph(subnetwork_original)
        if len(np.intersect1d(subnetwork_candidate.nodes(),TS_mutation))!=0:
            conductance_score_candidate=conductance_score_function(subnetwork_candidate,background_network)
            diff = conductance_score_candidate-conductance_score_curr
            t = temperature_all[time_eg]
            if (-diff / t)<100:
                metropolis = math.exp(-diff / t)
            else:
                metropolis=1
            random_num=random.random()
            enrichment_testing_subnetwork=enrichment_analysis_tsgenes(background_network_nodes,subnetwork_candidate,ts_gene)
            enrichment_testing_subnetwork_pvalue=float(enrichment_testing_subnetwork['pvalue'])
            if  (random_num< metropolis) & (enrichment_testing_subnetwork_pvalue<pvalue_cutoff) :
                print(time_eg)
                subnetwork_curr,conductance_score_curr, subnetwork_nodes_status_df_curr= nx.Graph(subnetwork_candidate),conductance_score_candidate,subnetwork_nodes_status_df_candidate
                subnetwork_stas_curr_eg=pd.DataFrame()
                subnetwork_stas_curr_eg['Times']=[time_eg]
                subnetwork_stas_curr_eg['Conductance score']=[conductance_score_candidate]
                subnetwork_stas_curr_eg['Conductance score difference']=[diff]
                subnetwork_stas_curr_eg['Probability']=[metropolis]
                subnetwork_stas_curr_eg['Random number']=[random_num]
                subnetwork_stas_curr_eg['Temperature']=[t]
                subnetwork_stas_curr_eg['P value']=[enrichment_testing_subnetwork_pvalue]
                subnetwork_stas_curr_eg['Number of subnetwork members']=[len(subnetwork_curr.nodes())]
                subnetwork_stas_curr_eg['Rotation decision']=['apply']
            else:
                print(f'{time_eg} rotation')
                subnetwork_stas_curr_eg=pd.DataFrame()
                subnetwork_stas_curr_eg['Times']=[time_eg]
                subnetwork_stas_curr_eg['Conductance score']=[conductance_score_candidate]
                subnetwork_stas_curr_eg['Conductance score difference']=[diff]
                subnetwork_stas_curr_eg['Probability']=[metropolis]
                subnetwork_stas_curr_eg['Random number']=[random_num]
                subnetwork_stas_curr_eg['Temperature']=[t]
                subnetwork_stas_curr_eg['P value']=[enrichment_testing_subnetwork_pvalue]
                subnetwork_stas_curr_eg['Number of subnetwork members']=[len(subnetwork_curr.nodes())]
                subnetwork_stas_curr_eg['Rotation decision']=['refuse']
                candidate_status=subnetwork_nodes_status_df_candidate.iloc[random_position[0],:]['status']
                if candidate_status==0:
                    subnetwork_nodes_status_df_candidate.iloc[random_position[0],1]=1
                else:
                    subnetwork_nodes_status_df_candidate.iloc[random_position[0],1]=0
        else:
            print(f'{time_eg} rotation')
            subnetwork_stas_curr_eg=pd.DataFrame()
            subnetwork_stas_curr_eg['Times']=[time_eg]
            subnetwork_stas_curr_eg['Conductance score']=[conductance_score_candidate]
            subnetwork_stas_curr_eg['Conductance score difference']=[diff]
            subnetwork_stas_curr_eg['Probability']=[metropolis]
            subnetwork_stas_curr_eg['Random number']=[random_num]
            subnetwork_stas_curr_eg['Temperature']=[t]
            subnetwork_stas_curr_eg['P value']=[enrichment_testing_subnetwork_pvalue]
            subnetwork_stas_curr_eg['Number of subnetwork members']=[len(subnetwork_curr.nodes())]
            subnetwork_stas_curr_eg['Rotation decision']=['refuse']
            candidate_status=subnetwork_nodes_status_df_candidate.iloc[random_position[0],:]['status']
            if candidate_status==0:
                subnetwork_nodes_status_df_candidate.iloc[random_position[0],1]=1
            else:
                subnetwork_nodes_status_df_candidate.iloc[random_position[0],1]=0
        subnetwork_stas_curr_all=pd.concat([subnetwork_stas_curr_all,subnetwork_stas_curr_eg])
    return(subnetwork_stas_curr_all,final_subnetwork_output_all)


def simulated_annealing_all_function_reversion(tissue,TS_mutation,work_dir):
    iteration_times=5000
    temperatureFunction_num=temperatureFunction(iteration_times,10,0.95)
    try:
        os.mkdir(f'{work_dir}/OncoNiche_pipeline/OncoNiche_output/{tissue}')
    except:
        pass
    background_network=nx.Graph(load_network_from_file(f'{work_dir}/OncoNiche_pipeline/uKIN_pipeline/Global_Gini_CGC_driver/background_network/Nested Systems in Tumors network.tsv'))
    ukin_number=100
    cancer_related_subnetwork=ukin_seed(tissue,TS_mutation,background_network,ukin_number,work_dir)
    background_network_nodes=list(background_network.nodes())
    ts_mut_gene=pd.read_csv(f'{work_dir}/OncoNiche_pipeline/uKIN_pipeline/Global_Gini_CGC_driver/{tissue}_global_Gini/tissue_mut_score',sep='\t',header=None)[1].tolist()
    ts_exp_gene=pd.read_csv(f'{work_dir}/OncoNiche_pipeline/uKIN_pipeline/Global_Gini_CGC_driver/{tissue}_global_Gini/tissue_exp_score',sep='\t',header=None)[1].tolist()
    ts_exp_mut_gene=np.union1d(ts_mut_gene,ts_exp_gene)
    ts_gene=np.intersect1d(ts_exp_mut_gene,background_network_nodes)
    enrich_cutoff=0.01
    simulated_annealing_main_output=simulated_annealing_main_pvalue_score_rotation(TS_mutation,nx.Graph(cancer_related_subnetwork),iteration_times,ts_gene,nx.Graph(background_network),enrich_cutoff,background_network_nodes,temperatureFunction_num)
    subnetwork_stas_curr_all=simulated_annealing_main_output[0]
    final_subnetwork_output=simulated_annealing_main_output[1]
    final_subnetwork_output['tissue']=tissue
    subnetwork_stas_curr_all['tissue']=tissue
    subnetwork_stas_curr_all['ts_seed']=TS_mutation
    subnetwork_stas_curr_all.to_csv(f'{work_dir}/OncoNiche_pipeline/OncoNiche_output/{tissue}/subnetwork_stas_curr_all_{tissue}_{TS_mutation}.txt',sep='\t',index=False)
    final_subnetwork_output.to_csv(f'{work_dir}/OncoNiche_pipeline/OncoNiche_output/{tissue}/final_subnetwork_output_{tissue}_{TS_mutation}.txt',sep='\t',index=False)

def ts_subnetwork_filter(tissue,work_dir):
    all_seed=np.unique([i.replace(f'{tissue}_','').split('_')[0] for i in os.listdir(f'{work_dir}/OncoNiche_pipeline/OncoNiche_output/{tissue}') if 'argument' in i] )
    sa_seed_network_all=pd.DataFrame()
    for seed_gene in all_seed:
        subnet_read=pd.read_csv(f'{work_dir}/OncoNiche_pipeline/OncoNiche_output/{tissue}/{tissue}_{seed_gene}_subnetwork_member_genes.txt',sep='\t')
        stas_read=pd.read_csv(f'{work_dir}/OncoNiche_pipeline/OncoNiche_output/{tissue}/{tissue}_{seed_gene}_subnetwork_argument.txt',sep='\t')
        stas_read_apply=stas_read[stas_read['Rotation decision']=='apply']
        if (stas_read.shape[0]!=0) & (stas_read_apply.shape[0]!=0):
            subnet_read_final=subnet_read[subnet_read['Times']==stas_read_apply['Times'].iloc[-1]]
            sa_seed_network=pd.DataFrame()
            sa_seed_network['Subnetwork member genes']=list(subnet_read_final['Subnetwork member genes'].str.split('|',expand=True).T.iloc[:,0])
            sa_seed_network['Tissue-specific mutated genes']=stas_read_apply['Tissue-specific mutated genes'].iloc[0]
            sa_seed_network['P value']=stas_read_apply['P value'].iloc[-1]
            sa_seed_network['Conductance score']=stas_read_apply['Conductance score'].iloc[-1]
            sa_seed_network['Number of subnetwork members']=stas_read_apply['Number of subnetwork members'].iloc[-1]
            sa_seed_network['Tissue']=tissue
            if sa_seed_network.shape[0]>10:
                sa_seed_network_all=pd.concat([sa_seed_network_all,sa_seed_network])
    sa_seed_network_all=sa_seed_network_all.drop_duplicates()
    sa_seed_network_all.to_csv(f'{work_dir}/OncoNiche_pipeline/OncoNiche_output/{tissue}/All_tissue_specific_subnetworks_in_{tissue}.txt',sep='\t',index=False)
    return(sa_seed_network_all)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", help="Please input tissue types", type=str)
    parser.add_argument("--work_dir", help="Please the working directory path where the code is located", type=str)
    args = parser.parse_args()
    tissue=args.tissue
    work_dir=args.work_dir
    TS_mutation_all=[ i.split('output_')[1].split('_results.txt')[0] for i in os.listdir(f'{work_dir}/OncoNiche_pipeline/uKIN_pipeline/output/uKIN_seed_SIA_CGC_{tissue}') ]
    for TS_mutation_eg in TS_mutation_all:
        simulated_annealing_all_function_reversion_run=simulated_annealing_all_function_reversion(tissue,TS_mutation_eg,work_dir)
    ts_subnetwork_process_output=ts_subnetwork_filter(tissue,work_dir)

