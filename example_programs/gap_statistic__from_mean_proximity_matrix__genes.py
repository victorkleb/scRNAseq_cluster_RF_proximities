
################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program: gap_statistic__from_mean_proximity_matrix__genes.py                                                        ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION  gap
#
################################################################################################################################



import pandas as pd
import numpy as np



 
from pathlib import Path
import os


import sys
sys.path.append("D:/scRNA-seq/Spec_clust_RFproximities_scRNAseq")


from  utilities import *
from FUNCTIONS_Spec_clust_RFproximities_scRNAseq import *



import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 
########################################################################################
   
consensus_clusterings_pkl =  "consensus_clusterings.pkl"
mean_proximity_pkl =  "mean proximity - genes.pkl"

logfile_txt = "gap_statistic_from_mean_proximity_matrix__genes.txt"
gap_pdf = "gap_statistic_from_mean_proximity_matrix__genes.pdf"
gap_table_pkl = "gap_statistic_from_mean_proximity_matrix__genes.pkl"


data_directory  = "Zhengmix4eq"
data_folder = r"D:/scRNA-seq/" + data_directory

trees_data_folder = data_folder + "/gene_data/1000 trees" 
trees_data_path = Path( trees_data_folder )

	

####  log output
logfile_dsn  =  trees_data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# PDF output  
plot_dsn = trees_data_path / gap_pdf
pdf_pages = PdfPages( plot_dsn )

# pickle output
gap_table_dsn = trees_data_path / gap_table_pkl


# pickle inputs 
consensus_clusterings_dsn = trees_data_path / consensus_clusterings_pkl
mean_proximity_dsn = trees_data_path / mean_proximity_pkl
 
########################################################################################

df_consensus_clusterings = pd.read_pickle ( consensus_clusterings_dsn )    
plog ( logfile, '\n df_consensus_clusterings : \n\n', df_consensus_clusterings )  
 
df_mean_proximity = pd.read_pickle ( mean_proximity_dsn )  
plog ( logfile, '\n df_mean_proximity : \n\n', df_mean_proximity )  
 
  
  
result_list = gap ( df_mean_proximity, df_consensus_clusterings )

df_gap = result_list[0]   
shuffles = result_list[1]   
plog ( logfile,'\n\n df_gap: \n\n',   df_gap )
 

 
n_genes = df_consensus_clusterings.shape[0]
n_synth_copies_clusterings = 100
n_trees = 1000


title0 = 'scRNA-seq data:  gap statistic for gene clusters -- spectral clusters and distances for gap derived from random forest proximities \n'
title1 = data_directory + ' data,  filtered with binomial deviance to select  ' + str ( n_genes ) + ' genes \n'
title2 = str( n_synth_copies_clusterings ) +  '  sets of synthetic data used for random forest classifiers, each using  ' + str(  n_trees ) + '  trees\n'
title3 = 'Gap calculated using  ' + str ( shuffles ) + ' shuffles'
 
 
titl = title0 + title1 + title2  + title3
  
fig = plt.figure()
plt.title ( titl, fontsize=5.5 )


plt.plot ( df_gap.index.values,  df_gap['gap'], color='k', label= 'Gap' )
plt.plot ( df_gap.index.values,  df_gap['gap - std'], color='k', linestyle=':',  label= 'Gap - std' )
plt.plot ( df_gap.index.values,  df_gap['gap + std'], color='k', linestyle=':',  label= 'Gap + std' )
plt.plot ( df_gap.index.values,  df_gap['gap - 2*std'], color='k', linestyle='--',  label= 'Gap - 2*std' )
plt.plot ( df_gap.index.values,  df_gap['gap + 2*std'], color='k', linestyle='--',  label= 'Gap + 2*std' )


plt.legend( loc='right' , prop={'size': 6}  )  
plt.xlabel ( 'number of gene clusters ', fontsize=7 )
plt.ylabel ( 'Gap statistic', fontsize=7 )
plt.xticks( fontsize= 7 )
plt.yticks( fontsize= 7 )

pdf_pages.savefig(fig) 



pdf_pages.close()   


df_gap.to_pickle ( gap_table_dsn ) 
 
 
logfile.close()  


################################################################################################################################ 
##                                                                                                                            ##
##  end program: gap_statistic__from_mean_proximity_matrix__genes.py                                                          ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################

