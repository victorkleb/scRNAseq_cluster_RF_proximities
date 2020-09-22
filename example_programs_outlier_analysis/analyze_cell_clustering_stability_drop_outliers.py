



################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program:  analyze_cell_clustering_stability_drop_outliers.py                                                        ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTIONS		ARI_and_Misclassification_Error_for_individual_clusteringss
#						ARI_and_Misclassification_Error_for_individual_clusterings_wrt_consensus_clusters
#						ARI_and_Misclassification_Error_for_subset_consensus_clusterings
#						
#
################################################################################################################################

import pandas as pd
import numpy  as np


import sys
sys.path.append("D:/scRNA-seq/Spec_clust_RFproximities_scRNAseq")

from  utilities import *
from FUNCTIONS_Spec_clust_RFproximities_scRNAseq import *

from pathlib import Path



import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 
#########################################################################################

pctl_list = [.01,.05, .10, .25, .5, .75, .90, .95, .99 ]

########################################################################################

individual_clusterings_pkl =  "individual_clusterings.pkl"
consensus_clusterings_pkl =  "consensus_clusterings.pkl"
consensus_clusterings_H_pkl =  "consensus_clusterings_H.pkl"
consensus_clusterings_Q_pkl =  "consensus_clusterings_Q.pkl"

ARI_pkl =  "adjusted_Rand_index_individual_clusterings.pkl"
ME_pkl =   "misclassification_error_individual_clusterings.pkl"
cell_ME_pkl = "cell_misclassification_error_individual_clusterings.pkl"

ARI_CC_pkl =  "adjusted_Rand_index_individual_clusterings_wrt_consensus.pkl"
ME_CC_pkl =   "misclassification_error_individual_clusterings._wrt_consensus.pkl"
cell_ME_CC_pkl = "cell_misclassification_error_individual_clusterings_wrt_consensus.pkl"

ARI_ME_pdf =  "cell clustering stability.pdf"



	

data_directory  = "Zhengmix8eq"
data_folder = r"D:/scRNA-seq/" + data_directory


cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 


trees_data_folder =  cell_data_folder  +  "/2000 trees drop outliers 8 clusters" 
trees_data_path = Path ( trees_data_folder )


	

###  log output
logfile_dsn  =  trees_data_path /  "analyze_cell_clustering_stability.txt"
logfile = open ( logfile_dsn,'w')

###  PDF output
plot_dsn = trees_data_path / ARI_ME_pdf
pdf_pages = PdfPages( plot_dsn )

# pickle outputs
ARI_dsn = trees_data_path / ARI_pkl
ME_dsn = trees_data_path / ME_pkl 
cell_ME_dsn = trees_data_path / cell_ME_pkl  
ARI_CC_dsn = trees_data_path / ARI_CC_pkl
ME_CC_dsn = trees_data_path / ME_CC_pkl 
cell_ME_CC_dsn = trees_data_path / cell_ME_CC_pkl  



# pickle inputs
individual_clusterings_dsn = trees_data_path / individual_clusterings_pkl
consensus_clusterings_dsn = trees_data_path / consensus_clusterings_pkl 
consensus_clusterings_H_dsn = trees_data_path / consensus_clusterings_H_pkl  
consensus_clusterings_Q_dsn = trees_data_path / consensus_clusterings_Q_pkl  

#######################################################################################

df_individual_clusterings= pd.read_pickle ( individual_clusterings_dsn )  
df_consensus_clusterings_Q= pd.read_pickle ( consensus_clusterings_Q_dsn )
df_consensus_clusterings_H= pd.read_pickle ( consensus_clusterings_H_dsn )
df_consensus_clusterings= pd.read_pickle ( consensus_clusterings_dsn )    

plog ( logfile, '\n\n\n df_individual_clusterings : \n\n', df_individual_clusterings )  
plog ( logfile, '\n\n\n df_consensus_clusterings_Q: \n\n', df_consensus_clusterings_Q )
plog ( logfile, '\n\n\n df_consensus_clusterings_H: \n\n', df_consensus_clusterings_H )
plog ( logfile, '\n\n\n df_consensus_clusterings: \n\n', df_consensus_clusterings ) 

pdline( logfile ) 


  
result_list =  ARI_and_Misclassification_Error_for_individual_clusterings ( df_individual_clusterings )
  
df_Adjusted_Rand_Index_individual_clusterings = result_list[0]
df_Misclassification_Error_individual_clusterings = result_list[1]  
df_cell_misclassification_error_individual_clusterings = result_list[2]  
 
plog ( logfile, '\n\n df_Adjusted_Rand_Index_individual_clusterings: \n\n', df_Adjusted_Rand_Index_individual_clusterings )  
plog ( logfile, '\n\n\n df_Misclassification_Error_individual_clusterings: \n\n', df_Misclassification_Error_individual_clusterings )
plog ( logfile, '\n\n\n df_cell_misclassification_error_individual_clusterings: \n\n', df_cell_misclassification_error_individual_clusterings ) 

plog ( logfile, '\n\n\n distribution of Adjusted Rand Index: \n\n', df_Adjusted_Rand_Index_individual_clusterings.describe ( percentiles = pctl_list ).transpose() )
plog ( logfile, '\n\n\n distribution of Misclassification Error: \n\n', df_Misclassification_Error_individual_clusterings.describe ( percentiles = pctl_list ).transpose() )  
plog ( logfile, '\n\n\n distribution of CELL misclassification error: \n\n', df_cell_misclassification_error_individual_clusterings.describe ( percentiles = pctl_list ).transpose() )    
  
pdline( logfile )    

  
 
 
result_list =  ARI_and_Misclassification_Error_for_individual_clusterings_wrt_consensus_clusters ( df_individual_clusterings, df_consensus_clusterings )

df_Adjusted_Rand_Index_to_CC = result_list[0]
df_Misclassification_Error_to_CC = result_list[1]    
df_cell_misclassification_error_to_CC = result_list[2]  

  
plog ( logfile, '\n\n df_Adjusted_Rand_Index_to_CC: \n\n', df_Adjusted_Rand_Index_to_CC )
plog ( logfile, '\n\n\n df_Misclassification_Error_to_CC: \n\n', df_Misclassification_Error_to_CC )
plog ( logfile, '\n\n\n df_cell_misclassification_error_to_CC: \n\n', df_cell_misclassification_error_to_CC )

plog ( logfile, '\n\n\n distribution of Adjusted Rand Index: each individual clustering to consensus clusters: \n\n', df_Adjusted_Rand_Index_to_CC.describe ( percentiles = pctl_list ).transpose() )
plog ( logfile, '\n\n\n distribution of Misclassification Error: each individual clustering to consensus clusters: \n\n', df_Misclassification_Error_to_CC.describe ( percentiles = pctl_list ).transpose() )  
plog ( logfile, '\n\n\n distribution of CELL Misclassification Error: each individual clustering to consensus clusters: \n\n', df_cell_misclassification_error_to_CC.describe ( percentiles = pctl_list ).transpose() )  

pdline( logfile )   
 



plog ( logfile,  '\n\n\ncompare consensus clusterings created from 2-way split of all spectral clusterings' )   
result_list =  ARI_and_Misclassification_Error_for_subset_consensus_clusterings ( df_consensus_clusterings_H )
  
df_Adjusted_Rand_Index = result_list[0]
df_Misclassification_Error = result_list[1]  
 
plog ( logfile, '\n\n\n df_Adjusted_Rand_Index: \n\n', df_Adjusted_Rand_Index )  
plog ( logfile, '\n\n\n df_Misclassification_Error: \n\n', df_Misclassification_Error ) 
  
pdline( logfile )  
  
  
 
 
plog ( logfile,  '\n\n\n compare consensus clusterings created from 4-way split of all spectral clusterings' )   
result_list =  ARI_and_Misclassification_Error_for_subset_consensus_clusterings ( df_consensus_clusterings_Q )
  
df_Adjusted_Rand_Index = result_list[0]
df_Misclassification_Error = result_list[1]  
 
plog ( logfile, '\n\n df_Adjusted_Rand_Index: \n\n', df_Adjusted_Rand_Index )  
plog ( logfile, '\n\n\n df_Misclassification_Error: \n\n', df_Misclassification_Error ) 

plog ( logfile, '\n\n\n distribution of Adjusted Rand Index: \n\n', df_Adjusted_Rand_Index.describe ( percentiles = pctl_list ).transpose() )
plog ( logfile, '\n\n\n distribution of _Misclassification Error: \n\n', df_Misclassification_Error.describe ( percentiles = pctl_list ).transpose() )  
  
pdline( logfile )  



###############################################################################################################  

n_genes = 195

num_for_analysis =  1 + df_individual_clusterings['individual_clustering'].max().astype(int)
n_trees = 2000
n_computed_clusterings = df_Adjusted_Rand_Index_individual_clusterings.shape[1]


title0 = 'scRNA-seq data:  stability analysis for cell clusters  -- spectral clusters derived from random forest proximities \n'
title1 = data_directory + ' data,  filtered with binomial deviance to select  ' + str ( n_genes ) + ' genes \n'
title2 = str( num_for_analysis ) +  ' synthetic data sets used for random forest classifiers, each using  ' + str(  n_trees ) + '  trees\n'
title3 = 'distribution of Adjusted Rand index for pairs of clusterings'   
titl = title0 + title1 + title2 + title3   
  
fig = plt.figure()
plt.title ( titl, fontsize=5.5 )

boxplot_data_list = []  
for clusters in range ( 2,  n_computed_clusterings + 2 ):
  ARI = df_Adjusted_Rand_Index_individual_clusterings [ clusters ].values.tolist()
  boxplot_data_list.append ( ARI )
	
plt.violinplot( boxplot_data_list, showextrema = False,   showmedians = True )
 
plt.xlabel ( 'number of clusters', fontsize=6 )
plt.ylabel ( 'Adjusted Rand Index', fontsize=6 )

plt.xticks( fontsize= 5 )
plt.yticks( fontsize= 5 )

plt.xticks( list ( range ( 1,  n_computed_clusterings + 1) ), list ( range  ( 2,  n_computed_clusterings + 2 ) ) )
plt.hlines( y=np.arange(.1, 1, .1), xmin=0, xmax=n_computed_clusterings + 2, linestyle=':', color='grey', linewidth=.5 )

plt.ylim( 0, 1 )
  
pdf_pages.savefig(fig)  
 
  
  
  
title3 = 'distribution of Misclassification Error for pairs of clusterings'  
titl = title0 + title1 + title2 + title3   
  
fig = plt.figure()
plt.title ( titl, fontsize=5.5 )

boxplot_data_list = []  
for clusters in range ( 2,  n_computed_clusterings + 2 ):
  ME = df_Misclassification_Error_individual_clusterings [ clusters ].values.tolist()
  boxplot_data_list.append ( ME )
	
plt.violinplot( boxplot_data_list, showextrema = False,   showmedians = True )
 
plt.xlabel ( 'number of clusters', fontsize=6 )
plt.ylabel ( 'fraction misclassified', fontsize=6 )

plt.xticks( fontsize= 5 )
plt.yticks( fontsize= 5 )

plt.xticks( list ( range ( 1,  n_computed_clusterings + 1) ), list ( range  ( 2,  n_computed_clusterings + 2 ) ) )
plt.hlines( y=np.arange(.1, 1, .1), xmin=0, xmax=n_computed_clusterings + 2, linestyle=':', color='grey', linewidth=.5 )

plt.ylim( 0, 1 )
  
pdf_pages.savefig(fig)  
 
  
  


title3 = 'distribution of Misclassification Error: individual to consensus clusterings'  
titl = title0 + title1 + title2 + title3   
  
fig = plt.figure()
plt.title ( titl, fontsize=5.5 )

boxplot_data_list = []  
for clusters in range ( 2,  n_computed_clusterings + 2 ):
  ME = df_Misclassification_Error_to_CC [ clusters ].values.tolist()
  boxplot_data_list.append ( ME )
	
plt.violinplot( boxplot_data_list, showextrema = False,   showmedians = True )
 
plt.xlabel ( 'number of clusters', fontsize=6 )
plt.ylabel ( 'fraction misclassified', fontsize=6 )

plt.xticks( fontsize= 5 )
plt.yticks( fontsize= 5 )

plt.xticks( list ( range ( 1,  n_computed_clusterings + 1) ), list ( range  ( 2,  n_computed_clusterings + 2 ) ) )
plt.hlines( y=np.arange(.1, 1, .1), xmin=0, xmax=n_computed_clusterings + 2, linestyle=':', color='grey', linewidth=.5 )

plt.ylim( 0, 1 )
  
pdf_pages.savefig(fig)  
 
  
  
title3 = "distribution of indiviual cells' Misclassification Errors: all pairs of individual clusterings"
titl = title0 + title1 + title2 + title3   
  
fig = plt.figure()
plt.title ( titl, fontsize=5.5 )

boxplot_data_list = []  
for clusters in range ( 2,  n_computed_clusterings + 2 ):
  ME = df_cell_misclassification_error_individual_clusterings [ clusters ].values.tolist()
  boxplot_data_list.append ( ME )
	
plt.violinplot( boxplot_data_list, showextrema = False,   showmedians = True )
 
plt.xlabel ( 'number of clusters', fontsize=6 )
plt.ylabel ( 'fraction misclassified', fontsize=6 )

plt.xticks( fontsize= 5 )
plt.yticks( fontsize= 5 )

plt.xticks( list ( range ( 1,  n_computed_clusterings + 1) ), list ( range  ( 2,  n_computed_clusterings + 2 ) ) )
plt.hlines( y=np.arange(.1, 1, .1), xmin=0, xmax=n_computed_clusterings + 2, linestyle=':', color='grey', linewidth=.5 )

plt.ylim( 0, 1 )
  
pdf_pages.savefig(fig)    
  

 
title3 = "distribution of indiviual cells' Misclassification Errors: individual clusterings to consensus clusters"
titl = title0 + title1 + title2 + title3   
  
fig = plt.figure()
plt.title ( titl, fontsize=5.5 )

boxplot_data_list = []  
for clusters in range ( 2,  n_computed_clusterings + 2 ):
  ME = df_cell_misclassification_error_to_CC [ clusters ].values.tolist()
  boxplot_data_list.append ( ME )
	
plt.violinplot( boxplot_data_list, showextrema = False,   showmedians = True )
 
plt.xlabel ( 'number of clusters', fontsize=6 )
plt.ylabel ( 'fraction misclassified', fontsize=6 )

plt.xticks( fontsize= 5 )
plt.yticks( fontsize= 5 )

plt.xticks( list ( range ( 1,  n_computed_clusterings + 1) ), list ( range  ( 2,  n_computed_clusterings + 2 ) ) )
plt.hlines( y=np.arange(.1, 1, .1), xmin=0, xmax=n_computed_clusterings + 2, linestyle=':', color='grey', linewidth=.5 )

plt.ylim( 0, 1 )
  
pdf_pages.savefig(fig)  


  
   
###############################################################################################################  

df_Adjusted_Rand_Index_individual_clusterings.to_pickle ( ARI_dsn )  
df_Misclassification_Error_individual_clusterings.to_pickle ( ME_dsn )
df_cell_misclassification_error_individual_clusterings.to_pickle ( cell_ME_dsn )

df_Adjusted_Rand_Index_to_CC.to_pickle ( ARI_CC_dsn )
df_Misclassification_Error_to_CC.to_pickle ( ME_CC_dsn )
df_cell_misclassification_error_to_CC.to_pickle ( cell_ME_CC_dsn )




pdf_pages.close()

logfile.close()  



################################################################################################################################ 
##                                                                                                                            ##
##  end program:  analyze_cell_clustering_stability_drop_outliers.py                                                          ##                                                                                 
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################

