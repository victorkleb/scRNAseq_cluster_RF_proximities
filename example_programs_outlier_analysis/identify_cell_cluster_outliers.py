



################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program:  identify_cell_cluster_outliers.py                                                                         ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
import pandas as pd
import numpy  as np


import sys
sys.path.append("D:/scRNA-seq/Spec_clust_RFproximities_scRNAseq")

from  utilities import *

from pathlib import Path



pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 
#########################################################################################

cell_clusters = 8



logfile_txt = "cell_cluster_" + str ( cell_clusters ) + "_outliers_MME_0p1.txt"
cell_outliers_pkl = "cell_cluster_" + str ( cell_clusters ) + "_outliers_MME_0p1.pkl"


cell_ME_pkl = "cell_misclassification_error_individual_clusterings.pkl"


data_directory  = "Zhengmix8eq"



data_folder = r"D:/scRNA-seq/" + data_directory
data_path = Path( data_folder )

	
cell_trees_data_folder = data_folder + "/cell_data - 5 gene clusters/2000 trees" 
cell_trees_data_path = Path( cell_trees_data_folder )
	

###  log output
logfile_dsn  =  cell_trees_data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# pickle output 
cell_outliers_dsn = cell_trees_data_path / cell_outliers_pkl


# pickle input
cell_ME_dsn = cell_trees_data_path / cell_ME_pkl 

#######################################################################################

pctl_list = [.01,.05, .10, .25, .5, .75, .90, .95, .99 ]

########################################################################################


df_cell_misclassification_error_individual_clusterings= pd.read_pickle ( cell_ME_dsn )   

plog ( logfile, '\n\n\n df_cell_misclassification_error_individual_clusterings : \n\n', df_cell_misclassification_error_individual_clusterings )   
 

df_cell_select_cluster = df_cell_misclassification_error_individual_clusterings[[ cell_clusters ]]
plog ( logfile, '\n\n\n df_cell_select_cluster.describe : \n\n', df_cell_select_cluster.describe( percentiles=pctl_list )  )

df_cell_select_outliers = df_cell_select_cluster.loc [ df_cell_select_cluster[ cell_clusters ] > 0.1 ] . rename ( columns= { cell_clusters: 'outliers' } )
plog ( logfile, '\n\n\n df_cell_select_outliers.describe : \n\n', df_cell_select_outliers.describe( percentiles=pctl_list )  )
plog ( logfile, '\n\n\n df_cell_select_outliers : \n\n', df_cell_select_outliers )




df_cell_select_outliers.to_pickle ( cell_outliers_dsn )




logfile.close()  



################################################################################################################################ 
##                                                                                                                            ##
##  end program:  identify_cell_cluster_outliers.py                                                                           ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################

