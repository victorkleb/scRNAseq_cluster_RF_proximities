
################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program:  Laplacian_scores_for_genes.py                                                                             ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION  Laplacian_scores
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
 
#########################################################################################    
  
pctl_list = [.01, .05, .1, .25, .5, .75, .9, .95, .99 ] 
      
########################################################################################

standardized_null_residuals_pkl = 'standardized_null_residuals.pkl'
cell_mean_proximities_pkl = "mean proximity - cells.pkl"
binomial_deviance_pkl = "binomial deviance - genuine data.pkl"


logfile_txt = "Laplacian_scores_for_genes.txt"
plot_pdf = "Laplacian_scores_vs_binomial_deviance.pdf"
Laplacian_scores_pkl = "Laplacian_scoress.pkl"


data_directory  = "Zhengmix4eq"
data_folder = r"D:/scRNA-seq/" + data_directory
data_path = Path( data_folder )

gene_trees_data_folder = data_folder + "/gene_data/1000 trees" 
gene_trees_data_path = Path( gene_trees_data_folder )

cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 
cell_trees_data_folder = cell_data_folder + "/2000 trees" 
cell_trees_data_path = Path( cell_trees_data_folder )


	

####  log output
logfile_dsn  =  cell_trees_data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# PDF output  
plot_dsn = cell_trees_data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn )

# pickle output
Laplacian_scores_dsn = cell_trees_data_path / Laplacian_scores_pkl


# pickle inputs
standardized_null_residuals_dsn = data_path / standardized_null_residuals_pkl	
cell_mean_proximities_dsn = cell_trees_data_path / cell_mean_proximities_pkl 

binomial_deviance_dsn = data_path / binomial_deviance_pkl

########################################################################################

df_proximity_array_mean = pd.read_pickle ( cell_mean_proximities_dsn )
plog ( logfile,  '\n\n df_proximity_array_mean: \n\n', df_proximity_array_mean )  
 
df_binomial_deviance = pd.read_pickle ( binomial_deviance_dsn ).rename ( columns={'genuine':'deviance'} )
plog ( logfile,  '\n\n df_binomial_deviance:\n' , df_binomial_deviance )

df_null_residuals_standardized = pd.read_pickle  ( standardized_null_residuals_dsn )
plog ( logfile,  '\n\n df_null_residuals_standardized:\n' , df_null_residuals_standardized )

pdline( logfile )



df_LS = Laplacian_scores ( df_proximity_array_mean, df_null_residuals_standardized )
 

df_stats =  df_binomial_deviance.merge ( df_LS, how='inner', left_index=True, right_index=True )
df_stats['log10 deviance'] = np.log10 ( df_stats['deviance']  )

plog ( logfile, '\n distribution of Laplacian scores: \n\n', df_stats[['Laplacian_score']].describe( percentiles = pctl_list ) )
plog ( logfile, '\n\n\n Spearman correlations of binomial deviance with Laplacian scores: \n\n', df_stats[['deviance', 'Laplacian_score']].corr ( method='spearman' ) )


n_genes = df_null_residuals_standardized.shape[0]

title1 = 'scRNA-seq data: ' + data_directory + ' data\n  plot binomial deviance vs Laplacian score \n'
title2 = 'calculated with random forest proximities \n ' + str ( n_genes ) + ' genes with largest binomial deviance -- grouped in 5 clusters \n'
 
fig, ax1 = plt.subplots()  

titl = title1 + title2
plt.title( titl, fontsize=5.5 )
  
ax1.scatter (  df_stats [ 'Laplacian_score' ],	df_stats [ 'log10 deviance' ],  c='k', s=1  )
ax1.set_ylabel ( 'log10 ( binomial deviance )', fontsize=7 )	
ax1.set_xlabel ( 'Laplacian score' , fontsize=7 )	
ax1.tick_params(labelsize=6)     
      
pdf_pages.savefig(fig, transparent=True )


pdf_pages.close()  
 
 
 

df_LS.to_pickle ( Laplacian_scores_dsn )
 
logfile.close()        
    


################################################################################################################################ 
##                                                                                                                            ##
##  end program:  Laplacian_scores_for_genes.py                                                                               ##                                                     
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################

