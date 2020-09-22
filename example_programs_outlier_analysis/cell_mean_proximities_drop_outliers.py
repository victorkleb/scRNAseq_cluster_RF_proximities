

################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program: cell_mean_proximities_drop_outliers.py                                                                     ##                                                                        
##                                                                                                                            ## 
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

logfile_txt = "cell_mean_proximities_drop_8_outliers_MME_0p1.txt"

mean_proximity_pkl = "mean proximity - cells.pkl"
cell_outliers_pkl = "cell_cluster_8_outliers_MME_0p1.pkl"

	

data_directory  = "Zhengmix8eq"
data_folder = r"D:/scRNA-seq/" + data_directory
data_path = Path( data_folder )



cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 

trees_data_folder = cell_data_folder + "/2000 trees" 
trees_data_path = Path( trees_data_folder )

drop_outliers_trees_data_folder = cell_data_folder  +  "/2000 trees drop outliers 8 clusters" 
drop_outliers_trees_data_path = Path ( drop_outliers_trees_data_folder )

	

###  log output
logfile_dsn  =  trees_data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

#### pickle output
mean_proximity_out_dsn = drop_outliers_trees_data_path / mean_proximity_pkl


#### pickle inputs
cell_outliers_dsn = trees_data_path / cell_outliers_pkl
mean_proximity_in_dsn = trees_data_path / mean_proximity_pkl

#############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

df_cell_select_outliers = pd.read_pickle ( cell_outliers_dsn )
plog ( logfile, '\n\n\n df_cell_select_outliers : \n\n', df_cell_select_outliers )

df_cell_select_outliers['outliers'] = True
plog ( logfile, '\n\n\n df_cell_select_outliers : \n\n', df_cell_select_outliers )

df_mean_proximity = pd.read_pickle ( mean_proximity_in_dsn )
plog ( logfile, '\n\n\n df_mean_proximity : \n\n', df_mean_proximity )


df_cell = pd.DataFrame ( index = range ( df_mean_proximity.shape[0] ),  data = df_mean_proximity.index.values, columns=['cell'] ) 
df_cell_append_outliers = df_cell.merge ( df_cell_select_outliers, how='left', left_on = 'cell', right_index=True )
df_cell_append_outliers ['outliers'] = df_cell_append_outliers ['outliers'].fillna( False )

df_cell_drop_outliers = df_cell_append_outliers.loc [ ~ df_cell_append_outliers ['outliers'] ]. drop ( columns = ['outliers'] ).sort_index()
plog ( logfile, '\n\n\n df_cell_drop_outliers : \n\n', df_cell_drop_outliers )

list_cell_drop_outliers = df_cell_drop_outliers['cell'].values.tolist()

df_mean_proximity_drop_outliers = df_mean_proximity[ list_cell_drop_outliers ].loc [ list_cell_drop_outliers ] 
plog ( logfile, '\n\n\n df_mean_proximity_drop_outliers : \n\n', df_mean_proximity_drop_outliers )

df_mean_proximity_drop_outliers.to_pickle ( mean_proximity_out_dsn )
  

  

  
  
logfile.close()



################################################################################################################################ 
##                                                                                                                            ##
##  end program: cell_mean_proximities_drop_outliers.py                                                                       ##                                                                               
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################

