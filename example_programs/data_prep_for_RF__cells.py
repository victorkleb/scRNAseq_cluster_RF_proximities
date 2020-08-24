

################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program: data_prep_for_RF__cells.py                                                                                 ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION data_prep_for_RF_cells
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
 
standardized_null_residuals_pkl = 'standardized_null_residuals.pkl'
consensus_clusterings_pkl =  "consensus_clusterings.pkl"

gene_cluster_means_pkl = "gene_cluster_means.pkl"


data_directory  = "Zhengmix4eq"
data_folder = r"D:/scRNA-seq/" + data_directory
data_path = Path ( data_folder )

trees_data_folder = data_folder + "/gene_data/1000 trees" 
trees_data_path = Path( trees_data_folder )

cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 
cell_data_path = Path( cell_data_folder )

# Create folder if necessary	 
if not os.path.exists( cell_data_path ):
    os.mkdir( cell_data_path )

	

####  log output
logfile_dsn  =  cell_data_path /  "data_prep_for_RF__cells.txt"
logfile = open ( logfile_dsn,'w')

# pickle output
gene_cluster_means_dsn = cell_data_path / gene_cluster_means_pkl


#### pickle inputs 
standardized_null_residuals_dsn = data_path / standardized_null_residuals_pkl
consensus_clusterings_dsn = trees_data_path / consensus_clusterings_pkl 


#############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

df_null_residuals_standardized = pd.read_pickle ( standardized_null_residuals_dsn )
plog ( logfile, '\n\n df_null_residuals_standardized :\n\n', df_null_residuals_standardized )

df_consensus_clusterings= pd.read_pickle ( consensus_clusterings_dsn )    
plog ( logfile, '\n\n df_consensus_clusterings: \n\n', df_consensus_clusterings )




result_list = data_prep_for_RF_cells  ( df_null_residuals_standardized, df_consensus_clusterings, 5, cell_data_path )

df_gene_cluster_means = result_list[0]
plog ( logfile, '\n\n df_gene_cluster_means \n\n', df_gene_cluster_means )




df_gene_cluster_means.to_pickle ( gene_cluster_means_dsn )

 


logfile.close()



################################################################################################################################ 
##                                                                                                                            ##
##  end program: data_prep_for_RF__cells.py                                                                                   ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################

