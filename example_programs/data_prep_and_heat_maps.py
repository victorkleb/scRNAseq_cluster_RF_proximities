

################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program:  data_prep_and_heat_maps.py                                                                                ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
#
# invoke FUNCTION  data_prep_and_heat_maps
#
################################################################################################################################

import pandas as pd
import numpy  as np

 
from pathlib import Path
import os


import sys
sys.path.append("D:/scRNA-seq/Spec_clust_RFproximities_scRNAseq")


from  utilities import *
from FUNCTIONS_Spec_clust_RFproximities_scRNAseq import *



pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 
#########################################################################################

cell_cluster_order_pkl = "cell_cluster_order.pkl"
gene_cluster_order_pkl = "gene_cluster_order.pkl"
null_residuals_gene_cell_ordered_pkl = "null_residuals_gene_cell_ordered.pkl"


standardized_null_residuals_pkl = 'standardized_null_residuals.pkl'
consensus_clusterings_pkl =  "consensus_clusterings.pkl"
gene_mean_proximities_pkl = "mean proximity - genes.pkl"
cell_mean_proximities_pkl = "mean proximity - cells.pkl"


data_directory  = "Zhengmix4eq"
data_folder = r"D:/scRNA-seq/" + data_directory
data_path = Path( data_folder )

gene_trees_data_folder = data_folder + "/gene_data/1000 trees" 
gene_trees_data_path = Path( gene_trees_data_folder )

cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 
cell_trees_data_folder = cell_data_folder + "/2000 trees" 
cell_trees_data_path = Path( cell_trees_data_folder )




heatmap_data_folder = data_folder + "/heat_map 5 gene, 5 cell clusters" 
heatmap_data_path = Path( heatmap_data_folder )

### Create folder if it doesn't exist
if not os.path.exists( heatmap_data_path ):
    os.mkdir( heatmap_data_path )
		

 
#  log output
logoutfile = heatmap_data_path / "data_prep_and_heat_maps.txt"
logfile = open ( logoutfile,'w')

# pickle outputs 
cell_cluster_order_dsn = heatmap_data_path / cell_cluster_order_pkl
gene_cluster_order_dsn = heatmap_data_path / gene_cluster_order_pkl
null_residuals_gene_cell_ordered_dsn = heatmap_data_path / null_residuals_gene_cell_ordered_pkl


# pickle inputs
standardized_null_residuals_dsn = data_path / standardized_null_residuals_pkl	
gene_consensus_clusterings_dsn = gene_trees_data_path / consensus_clusterings_pkl 
cell_consensus_clusterings_dsn = cell_trees_data_path / consensus_clusterings_pkl 
gene_mean_proximities_dsn = gene_trees_data_path / gene_mean_proximities_pkl 
cell_mean_proximities_dsn = cell_trees_data_path / cell_mean_proximities_pkl 

#######################################################################################

plt.ioff()

     
df_null_residuals_standardized = pd.read_pickle ( standardized_null_residuals_dsn )
plog ( logfile, '\n\n df_null_residuals_standardized: \n\n', df_null_residuals_standardized ) 

   
  
df_cell_consensus_clusterings= pd.read_pickle ( cell_consensus_clusterings_dsn )     
plog ( logfile,  '\n\n df_cell_consensus_clusterings: \n\n', df_cell_consensus_clusterings )

df_cell_mean_proximities = pd.read_pickle ( cell_mean_proximities_dsn )
plog ( logfile,  '\n\n df_cell_mean_proximities: \n\n', df_cell_mean_proximities )  

df_gene_consensus_clusterings= pd.read_pickle ( gene_consensus_clusterings_dsn )     
plog ( logfile,  '\n\n df_gene_consensus_clusterings: \n\n', df_gene_consensus_clusterings )

df_gene_mean_proximities = pd.read_pickle ( gene_mean_proximities_dsn )
plog ( logfile,  '\n\n df_gene_mean_proximities: \n\n', df_gene_mean_proximities )  
pdline( logfile )


######################################################################################
  
  
result_list = data_prep_and_heat_maps  ( 'Zhengmix4eq',  5, 5, df_null_residuals_standardized, df_gene_consensus_clusterings, df_cell_consensus_clusterings, df_gene_mean_proximities, df_cell_mean_proximities, heatmap_data_path  )
 
 
pa_cluster_mean_list =  result_list[0]
permutation_list = result_list[1]
dataframe_list = result_list[2]

pa_cluster_genes_mean = pa_cluster_mean_list[0]
pa_cluster_cells_mean = pa_cluster_mean_list[1] 

gene_permutation = permutation_list[0]
cell_permutation = permutation_list[1]

df_gene_cluster_order = dataframe_list[0]
df_cell_cluster_order = dataframe_list[1]
df_null_residuals_gene_cell_ordered = dataframe_list[2]
 

plog ( logfile,  '\n\n\n proximity_array - gene cluster means: \n\n', pa_cluster_genes_mean )
plog ( logfile,  '\n\n\n proximity_array - cell cluster means: \n\n', pa_cluster_cells_mean ) 
 
plog ( logfile,  '\n\n\n gene_permutation \n\n', gene_permutation  )  
plog ( logfile,  '\n\n\n cell_permutation \n\n', cell_permutation  )  

plog ( logfile,  '\n\n\n df_gene_cluster_order \n\n', df_gene_cluster_order  )  
plog ( logfile,  '\n\n\n df_cell_cluster_order \n\n', df_cell_cluster_order  )  
plog ( logfile,  '\n\n\n df_null_residuals_gene_cell_ordered \n\n', df_null_residuals_gene_cell_ordered   )  



df_cell_cluster_order.to_pickle ( cell_cluster_order_dsn )
df_gene_cluster_order.to_pickle ( gene_cluster_order_dsn )
df_null_residuals_gene_cell_ordered.to_pickle  ( null_residuals_gene_cell_ordered_dsn ) 




logfile.close()  


################################################################################################################################ 
##                                                                                                                            ##
##  end program:  data_prep_and_heat_maps.py                                                                                  ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################


