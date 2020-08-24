
################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program:  cell_clusters__individual_and_consensus.py                                                                ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION  cell_clusters__individual_and_consensus
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

individual_clusterings_pkl =  "individual_clusterings.pkl"
consensus_clusterings_pkl =  "consensus_clusterings.pkl"
consensus_clusterings_H_pkl =  "consensus_clusterings_H.pkl"
consensus_clusterings_Q_pkl =  "consensus_clusterings_Q.pkl"


data_folder = r"D:/scRNA-seq/Zhengmix4eq"


cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 

trees_data_folder = cell_data_folder + "/2000 trees" 
trees_data_path = Path( trees_data_folder )

 

	
	

####  log output
logfile_dsn  =  trees_data_path /  "cell_clusters__individual_and_consensus.txt"
logfile = open ( logfile_dsn,'w')

# pickle outputs
individual_clusterings_dsn = trees_data_path / individual_clusterings_pkl
consensus_clusterings_dsn = trees_data_path / consensus_clusterings_pkl 
consensus_clusterings_H_dsn = trees_data_path / consensus_clusterings_H_pkl  
consensus_clusterings_Q_dsn = trees_data_path / consensus_clusterings_Q_pkl  

######################################################################################################################################

result_list =  cell_clusters__individual_and_consensus  ( 20, trees_data_path )
  
df_individual_clusterings = result_list[0]
df_consensus_clusterings = result_list[1]
df_consensus_clusterings_H = result_list[2]
df_consensus_clusterings_Q = result_list[3]  
  
  
  

plog ( logfile, '\n\n\n df_individual_clusterings : \n\n', df_individual_clusterings )  
plog ( logfile, '\n\n\n df_consensus_clusterings_Q: \n\n', df_consensus_clusterings_Q )
plog ( logfile, '\n\n\n df_consensus_clusterings_H: \n\n', df_consensus_clusterings_H )
plog ( logfile, '\n\n\n df_consensus_clusterings: \n\n', df_consensus_clusterings ) 


df_individual_clusterings.to_pickle ( individual_clusterings_dsn )  
df_consensus_clusterings_Q.to_pickle ( consensus_clusterings_Q_dsn )
df_consensus_clusterings_H.to_pickle ( consensus_clusterings_H_dsn )
df_consensus_clusterings.to_pickle ( consensus_clusterings_dsn )    
 

 


logfile.close()

################################################################################################################################ 
##                                                                                                                            ##
##  end program:  cell_clusters__individual_and_consensus.py                                                                  ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################


