

################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program: proximities_cells_drop_outliers.py                                                                         ##                                                                                
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

logfile_txt = "proximities_cells_drop_outliers.txt"
cell_outliers_pkl = "cell_cluster_8_outliers_MME_0p1.pkl"

prox_ds_base = "cell_proximities_"


data_directory  = "Zhengmix8eq"

data_folder = r"D:/scRNA-seq/" + data_directory


cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 
cell_data_path = Path( cell_data_folder )

trees_data_folder =  cell_data_folder +  "/2000 trees" 
trees_data_path = Path ( trees_data_folder )

drop_outliers_trees_data_folder =  cell_data_folder +  "/2000 trees drop outliers 8 clusters" 
drop_outliers_trees_data_path = Path ( drop_outliers_trees_data_folder )


	

###  log output
logfile_dsn  =  cell_data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')



# Create folders if necessary	 
if not os.path.exists( drop_outliers_trees_data_path ):
    os.mkdir( drop_outliers_trees_data_path )
	


#### pickle input
cell_outliers_dsn = trees_data_path / cell_outliers_pkl


#############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

df_cell_select_outliers = pd.read_pickle ( cell_outliers_dsn )
plog ( logfile, '\n\n\n df_cell_select_outliers : \n\n', df_cell_select_outliers )

df_cell_select_outliers['outliers'] = True
plog ( logfile, '\n\n\n df_cell_select_outliers : \n\n', df_cell_select_outliers )



for synth_copy in [0]:
  print ( 'synth_copy: ', synth_copy )
  
  prox_dsn = prox_ds_base + str ( synth_copy ) + ".csv"  
  prox_in_path_dsn	= trees_data_path / prox_dsn 
  prox_out_path_dsn	= drop_outliers_trees_data_path / prox_dsn  
  
  df_prox = pd.read_csv ( prox_in_path_dsn  )
  plog ( logfile, 'df_prox: \n\n', df_prox )

  
  df_cell = df_prox[[ 'cell' ]]
  df_cell_append_outliers = df_cell.merge ( df_cell_select_outliers, how='left', left_on = 'cell', right_index=True )
  df_cell_append_outliers ['outliers'] = df_cell_append_outliers ['outliers'].fillna( False )
 
  df_cell_drop_outliers = df_cell_append_outliers.loc [ ~ df_cell_append_outliers ['outliers'] ]. drop ( columns = ['outliers'] ).sort_index()
  plog ( logfile, '\n\n\n df_cell_drop_outliers : \n\n', df_cell_drop_outliers )

  list_cell_drop_outliers = df_cell_drop_outliers['cell'].values.tolist()

  df_prox_cell_index = df_prox . set_index ( ['cell'] )
  df_prox_drop_outliers = df_prox_cell_index[ list_cell_drop_outliers ].loc [ list_cell_drop_outliers ] 
  plog ( logfile, '\n\n\n df_prox_drop_outliers : \n\n', df_prox_drop_outliers )

  df_prox_drop_outliers.to_csv ( prox_out_path_dsn )
  


for synth_copy in range ( 1, 100 ): 
  print ( 'synth_copy: ', synth_copy ) 

  prox_dsn = prox_ds_base + str ( synth_copy ) + ".csv"  
  prox_in_path_dsn	= trees_data_path / prox_dsn 
  prox_out_path_dsn	= drop_outliers_trees_data_path / prox_dsn  
  
  df_prox = pd.read_csv ( prox_in_path_dsn  )
  plog ( logfile, 'df_prox: \n\n', df_prox )

  df_prox_cell_index = df_prox . set_index ( ['cell'] )
  df_prox_drop_outliers = df_prox_cell_index[ list_cell_drop_outliers ].loc [ list_cell_drop_outliers ] 

  df_prox_drop_outliers.to_csv ( prox_out_path_dsn )

  

  
  
logfile.close()



################################################################################################################################ 
##                                                                                                                            ##
##  end program: proximities_cells_drop_outliers.py                                                                           ##                                                                               
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################

