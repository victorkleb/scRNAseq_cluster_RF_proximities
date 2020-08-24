

################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program:  calculate_mean_proximity_matrix__cells.py                                                                 ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION   cells_mean_proximities
#
################################################################################################################################

import pandas as pd
import numpy  as np

 


  

import sys
sys.path.append("D:/scRNA-seq/Spec_clust_RFproximities_scRNAseq")


from  utilities import *
from FUNCTIONS_Spec_clust_RFproximities_scRNAseq import *


 
from pathlib import Path




pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 
#########################################################################################
########################################################################################

mean_proximity_pkl =  "mean proximity - cells.pkl"



data_directory  = "Zhengmix4eq"
data_folder = r"D:/scRNA-seq/" + data_directory

cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 

trees_data_folder = cell_data_folder + "/2000 trees" 
trees_data_path = Path( trees_data_folder )

	

####  log output
logfile_dsn  =  trees_data_path /  "calculate_mean_proximity_matrix__cells.txt"
logfile = open ( logfile_dsn,'w')

# pickle outputs
mean_proximity_dsn = trees_data_path / mean_proximity_pkl

######################################################################################################################################
 

result_list =  cells_mean_proximities  ( trees_data_path )  

df_mean_proximity = result_list[0]   
plog ( logfile, '\n df_mean_proximity : \n\n', df_mean_proximity )  
 

 
df_mean_proximity.to_pickle ( mean_proximity_dsn )  
  
 


logfile.close()



################################################################################################################################ 
##                                                                                                                            ##
##  end program:  calculate_mean_proximity_matrix__cells.py                                                                   ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################

