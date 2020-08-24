
################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program: calculate_mean_proximity_matrix__genes.py                                                                  ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION   genes_mean_proximities
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

mean_proximity_pkl =  "mean proximity - genes.pkl"


data_folder = r"D:/scRNA-seq/Zhengmix4eq"

trees_data_folder = data_folder + "/gene_data/1000 trees" 
trees_data_path = Path( trees_data_folder )

	

####  log output
logfile_dsn  =  trees_data_path /  "calculate_mean_proximity_matrix__genes.txt"
logfile = open ( logfile_dsn,'w')

# pickle outputs
mean_proximity_dsn = trees_data_path / mean_proximity_pkl

######################################################################################################################################
 

result_list =  genes_mean_proximities  ( trees_data_path )  

df_mean_proximity = result_list[0]   
plog ( logfile, '\n df_mean_proximity : \n\n', df_mean_proximity )  
 

 
df_mean_proximity.to_pickle ( mean_proximity_dsn )  
  
 


logfile.close()

################################################################################################################################ 
##                                                                                                                            ##
##  end program: calculate_mean_proximity_matrix__genes.py                                                                    ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################
