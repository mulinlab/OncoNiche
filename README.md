# OncoNiche
The OncoNiche model identifies contextual information associated with tissue-specific mutated genes. Initially, using the uKIN method, OncoNiche performs guided random walks and heat diffusion in a network propagation framework to detect cancer-risk genes influenced by these mutations, identifying these as cancer-related subnetworks. The connectivity of each subnetwork is quantified through conductance scores. By applying simulated annealing, OncoNiche identifies the subnetwork with the optimal global conductance score. At each iteration of the annealing process, we ensured significant enrichment of tissue-specific mutated and expressed genes within the refined subnetworks. Ultimately, OncoNiche identifies subnetworks closely connected with each tissue-specific mutated gene, which likely reflect the tissue context associated with each tissue-specific mutations.
## uKIN method
### 1. Input
There are three required input files:
1) a background network file 
2) a prior knowledge file containing a list of cancer driver genes, tissue-specific mutated and expressed genes
3) a file of tissue-specific mutated genes, each weighted by the Gini index.

Using skin cancer data as an example, we create a working directory with the structure outlined below and place the uKIN input data files within it.
```
work_dir          
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
```
python uKIN_pipeline.py --tissue tissue_type
```
## Simulated annealing
### 1. Input
### 2. Output
### 3. How to run
