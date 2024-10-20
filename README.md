# OncoNiche
The OncoNiche model identifies contextual information associated with tissue-specific mutated genes. Initially, using the uKIN method, OncoNiche performs guided random walks and heat diffusion in a network propagation framework to detect cancer-risk genes influenced by these mutations, identifying these as cancer-related subnetworks. The connectivity of each subnetwork is quantified through conductance scores. By applying simulated annealing, OncoNiche identifies the subnetwork with the optimal global conductance score. At each iteration of the annealing process, we ensured significant enrichment of tissue-specific mutated and expressed genes within the refined subnetworks. Ultimately, OncoNiche identifies subnetworks closely connected with each tissue-specific mutated gene, which likely reflect the tissue context associated with each tissue-specific mutations.
![image](https://github.com/mulinlab/OncoNiche/blob/main/OncoNiche_pipeline/img/OncoNiche.png)
## Identifying cancer-related subnetworks (uKIN method)
### 1. Input
There are three required input files:
1) a background network file 
2) a prior knowledge file containing a list of cancer driver genes, tissue-specific mutated and expressed genes
3) a file of tissue-specific mutated genes, each weighted by the Gini index.

Using skin cancer data as an example, we create a working directory with the structure outlined below and place the uKIN input data files within it.

```
work_dir          
└─OncoNiche_pipeline          
   └─uKIN_pipeline
      └─Global_Gini_CGC_driver
         ├─background_network
         │      Nested Systems in Tumors network.tsv    # background network file
         │      
         └─Skin_global_Gini
                 tissue_specific_genes.txt    # a file of tissue specific genes	
                 prior_knowledge.txt	# a prior knowledge file
```

### 2. Output
output_tissue_specific_genes_results.txt is written in the uKIN/output directory. The file contains a list of candidate genes ranked by how frequently they are visited as the guided walks reach the stationary distribution.
### 3. How to run
Run command:
```
python uKIN_pipeline.py --tissue tissue_type --work_dir work_dir
```
In uKIN_pipeline/output folder, you can find the uKIN_seed_SIA_CGC_tissue folder, which contains the candidate genes and visit frequencies associated with all tissue-specific mutated genes for the input tissue type.
## Identification of densely connected subnetworks (simulated annealing)
### 1. Input
There are three required input files:
 1) a background network file 
 2) a gene visit frequency file output from uKIN
 3) a list file of tissue-specific mutated and expressed genes
Using skin cancer data as an example, we create a working directory with the structure outlined below and place the simulated annealing input data files within it.
```
work_dir          
└─OncoNiche_pipeline
   └─uKIN_pipeline
       ├─Global_Gini_CGC_driver
       │  ├─background_network
       │  │      Nested Systems in Tumors network.tsv    # background network file
       │  │      
       │  └─Skin_global_Gini
       │          tissue_exp_score    # a list of tissue-specific expressed genes
       │          tissue_mut_score    # a list of tissue-specific mutated genes
       └─output
           └─uKIN_seed_SIA_CGC_Skin
                   output_tissue_specific_mutated genes_results.txt    # a list of candidate genes ranked by visit frequencies
```
### 2. Output
The file “Tissue_gene_subnetwork_member_genes.txt”, “Tissue_gene_subnetwork_argument.txt” and  "All_tissue_specific_subnetworks_in_tissue.txt" were written in the OncoNiche_output/tissue directory. 


The file "Tissue_gene_subnetwork_member_genes.txt" contains the genes of subnetwork members for each iteration, organized as follows.
| Tissue-specific mutated genes  | Subnetwork member genes | Times |
| ------------- | ------------- | ------------- |

The file "Tissue_gene_subnetwork_argument.txt" contains the parameter data calculated by the model for each iteration, organized as follows.
| Times  | Conductance score | Conductance score difference | Probability | Random number | Temperature | P value | Number of subnetwork members | Rotation decision |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |

The file "All_tissue_specific_subnetworks_in_tissue.txt" contains all tissue-specific subnetworks identified by OncoNiche, organized as follows.
| Subnetwork member genes  | Tissue-specific mutated genes | P value | Conductance score | Number of subnetwork members | Tissue |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |

```
work_dir          
└─OncoNiche_pipeline
   └─OncoNiche_output
       └─Skin
             Skin_BRAF_subnetwork_member_genes.txt  # the genes of subnetwork members for each iteration
             Skin_BRAF_subnetwork_argument.txt  # the parameter data calculated by the model for each iteration
             All_tissue_specific_subnetworks_in_tissue.txt  # all tissue-specific subnetworks identified by OncoNiche
```
### 3. How to run
Run command:
```
python simulated_annealing_pipeline.py --tissue tissue_type --work_dir work_dir
```
In the OncoNiche_output/tissue directory, you can find the files "Tissue_gene_subnetwork_member_genes.txt", "Tissue_gene_subnetwork_argument.txt" and "All_tissue_specific_subnetworks_in_tissue.txt", which contain the subnetwork member genes, the parameter data calculated by the model during each iteration and all tissue-specific subnetworks identified by OncoNiche in this tissue.
# Copyright
Copyright (c) Mulinlab@Tianjin Medical University 2021-2024.
# Citation
A network-based methodology for identifying molecular networks associated with tissue-specific mutations.
