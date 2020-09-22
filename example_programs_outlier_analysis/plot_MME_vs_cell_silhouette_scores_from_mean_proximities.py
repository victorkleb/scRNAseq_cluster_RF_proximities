


#######################################################################################################
#######################################################################################################
#                                                                                                     #   
#  start program  plot_MME_vs_cell_silhouette_scores_from_mean_proximities.py                         #
#                                                                                                     #  
#######################################################################################################



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
  
def pv_table (  dsin0, col1,   col2 ):
  # print( 'col1: ', col1, '    col2: ', col2 )
  dsin = dsin0.copy()
  dsin ['count'] = 1
  pt = pd.pivot_table( dsin, values='count',  index=[ col1 ], columns=[ col2 ], aggfunc=np.sum, margins=True )
  ptfna = pt.fillna(0)
  pti = ptfna.astype(int)
  print ( '\n\n\n ',pti ,  file = logfile )
  g, p, dof, expctd = chi2_contingency( pti, lambda_="log-likelihood")
  print ( '\ng-statistic: ', format(g, '.2f') , '     p-value: ', format(p, '.6f'), '\n', file = logfile )
  
  return pti

  
  
pctl_list = [.01,.05, .10, .25, .5, .75, .90, .95, .99 ]

########################################################################################
   
    
plot_pdf = "plot_MME_vs_cell_silhouette_scores_from_mean_proximities.pdf"   
logfile_txt = "plot_MME_vs_cell_silhouette_scores_from_mean_proximities.txt"

cell_silhouette_scores_pkl = "cell_silhouette_scores_from_mean_proximities.pkl"
cell_ME_pkl = "cell_misclassification_error_individual_clusterings.pkl"



data_directory  = "Zhengmix8eq"
data_folder = r"D:/scRNA-seq/" + data_directory



cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 

trees_data_folder =  cell_data_folder +  "/8000 trees" 
trees_data_path = Path ( trees_data_folder )
	

####  log output
logfile_dsn  =  trees_data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')

# PDF output  
plot_dsn = trees_data_path / plot_pdf
pdf_pages = PdfPages( plot_dsn )
 

# pickle inputs
cell_silhouette_scores_dsn = trees_data_path / cell_silhouette_scores_pkl
cell_ME_dsn = trees_data_path / cell_ME_pkl 
 
########################################################################################


df_cell_silhouette_scores_all = pd.read_pickle ( cell_silhouette_scores_dsn )
plog ( logfile, '\n\n df_cell_silhouette_scores_all: \n\n', df_cell_silhouette_scores_all )
 
df_cell_misclassification_error_individual_clusterings_all = pd.read_pickle ( cell_ME_dsn )   
plog ( logfile, '\n\n\n df_cell_misclassification_error_individual_clusterings_all : \n\n', df_cell_misclassification_error_individual_clusterings_all ) 
 

for clusters in  df_cell_silhouette_scores_all.columns.values.tolist(): 

  df_cell_silhouette_scores = df_cell_silhouette_scores_all[[ clusters ]].rename ( columns={ clusters: 'silhouette_score'} )
  df_cell_misclassification_error_individual_clusterings = df_cell_misclassification_error_individual_clusterings_all[[ clusters ]] 
 
  df_cell_silhouette_scores ['ss_pct'] = df_cell_silhouette_scores ['silhouette_score'].rank ( pct=True ) 
  df_cell_silhouette_scores['silhouette_score_quartile'] = np.floor ( 3.999999* df_cell_silhouette_scores ['ss_pct'] ).astype(int)
  # plog ( logfile,  '\n\n df_cell_silhouette_scores:\n' , df_cell_silhouette_scores ) 
 

  df_combo = df_cell_misclassification_error_individual_clusterings[[ clusters ]].merge ( df_cell_silhouette_scores, how='inner', left_index=True, right_index=True )
  # plog ( logfile,  '\n\n df_combo:\n' , df_combo )

  arr_quartile_upper_limits = df_combo [[ 'silhouette_score', 'silhouette_score_quartile' ]].groupby ( ['silhouette_score_quartile'] ).max() [[ 'silhouette_score' ]].values


  df_corr_sp = df_combo.corr ( method='spearman' )
  plog ( logfile,  '\n\n df_corr_sp:\n' , df_corr_sp )
 
 
  df_combo [ 'ME_bin' ] = pd.cut ( df_combo [ clusters ], [0, 0.01, .1, .2, 1], right=False ) 
  pti = pv_table ( df_combo, 'ME_bin', 'silhouette_score_quartile' )   
  
  

  n_cells = df_combo.shape [0]
 
  title1 = 'Zhengmix8eq scRNA-seq data - ' + str ( n_cells ) + ' cells \n' 
  title2 = 'Mean Misclassification Error for ' + str ( clusters ) + ' individual clusters (8000 trees) vs Silhouette score'
 
  fig, ax1 = plt.subplots()

  titl = title1 + title2  
  plt.title( titl, fontsize=5. )

  ax1.scatter ( df_combo [ 'silhouette_score' ],	df_combo [ clusters ],c='k', s=0.05)
  ax1.set_ylabel ( 'cell Mean Misclassification Error \n 4950 pairs of individual clusterings', fontsize=5. )	
  ax1.set_xlabel ( 'Silhouette score \n blue vertical lines separate quartiles' , fontsize=5. )	
  ax1.tick_params(labelsize=5. )

  for i in range ( 3 ):
    ax1.axvline( x= arr_quartile_upper_limits[i], linewidth = 0.1 )
	
  pdf_pages.savefig(fig, transparent=True )
  pdline ( logfile )




pdf_pages.close()  


 
logfile.close()  


#######################################################################################################
#                                                                                                     #   
#  end   program  plot_MME_vs_cell_silhouette_scores_from_mean_proximities.py                         #
#                                                                                                     #  
#######################################################################################################
#######################################################################################################



