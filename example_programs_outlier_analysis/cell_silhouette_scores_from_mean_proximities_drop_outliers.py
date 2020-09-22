

#########################################################################################
#########################################################################################
#                                                                                       #   
#  start program  cell_silhouette_scores_from_mean_proximities_drop_outliers.py         #
#                                                                                       #  
#########################################################################################



import pandas as pd
import numpy as np



 
from pathlib import Path
import os


import sys
sys.path.append("D:/scRNA-seq/Spec_clust_RFproximities_scRNAseq")


from  utilities import *
from FUNCTIONS_Spec_clust_RFproximities_scRNAseq import *


from sklearn.metrics import silhouette_samples, silhouette_score

 


pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 

#########################################################################################

pctl_list = [.01,.05, .10, .25, .5, .75, .90, .95, .99 ]

########################################################################################

   
consensus_clusterings_pkl =  "consensus_clusterings.pkl"
mean_proximity_pkl =  "mean proximity - cells.pkl"

logfile_txt = "cell_silhouette_scores_from_mean_proximities.txt"
cell_silhouette_scores_pkl = "cell_silhouette_scores_from_mean_proximities.pkl"



data_directory  = "Zhengmix8eq"
data_folder = r"D:/scRNA-seq/" + data_directory



cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 

trees_data_folder = cell_data_folder  +  "/2000 trees drop outliers 8 clusters" 
trees_data_path = Path( trees_data_folder )


	

####  log output
logfile_dsn  =  trees_data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

 

# pickle output
cell_silhouette_scores_dsn = trees_data_path / cell_silhouette_scores_pkl


# pickle inputs 
consensus_clusterings_dsn = trees_data_path / consensus_clusterings_pkl
mean_proximity_dsn = trees_data_path / mean_proximity_pkl
 
########################################################################################

df_consensus_clusterings = pd.read_pickle ( consensus_clusterings_dsn )    
plog ( logfile, '\n df_consensus_clusterings : \n\n', df_consensus_clusterings )  
 
df_mean_proximity = pd.read_pickle ( mean_proximity_dsn )  
plog ( logfile, '\n df_mean_proximity : \n\n', df_mean_proximity )  
 
arr_distance = 1 - df_mean_proximity.values
np.fill_diagonal ( arr_distance, 0. ) 

del ( df_mean_proximity )


df_silhouette_scores_list = []

indx = df_consensus_clusterings.index

for clusters in  df_consensus_clusterings.columns.values.tolist():

  cluster_labels = df_consensus_clusterings[ clusters ].values

  sample_silhouette_values = silhouette_samples( arr_distance, cluster_labels, metric = "precomputed" )
  df_cell_silhouette_scores = pd.DataFrame ( index=indx, data=sample_silhouette_values, columns = [ clusters ] )
  df_silhouette_scores_list.append ( df_cell_silhouette_scores )
  
 
df_silhouette_scores = pd.concat ( df_silhouette_scores_list, axis=1, sort=False )
 

plog ( logfile, '\n\n df_silhouette_scores: \n\n', df_silhouette_scores )
plog ( logfile, '\n\n distribution of silhouette scores: \n\n', df_silhouette_scores.describe ( percentiles = pctl_list ). transpose()  )

 
df_silhouette_scores.to_pickle ( cell_silhouette_scores_dsn )
 
 
 
logfile.close()  



#########################################################################################
#                                                                                       #   
#  end program  cell_silhouette_scores_from_mean_proximities_drop_outliers.py           #
#                                                                                       #  
#########################################################################################
#########################################################################################


