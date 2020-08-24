################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program: compare_cell_clusters_with_ground_truth.py                                                                 ##
##                                                                                                                            ## 
################################################################################################################################


import pandas as pd
import numpy  as np


from pathlib import Path

import sys
sys.path.append("D:/scRNA-seq/Spec_clust_RFproximities_scRNAseq")


from  utilities import *


from scipy.stats import chi2_contingency


from scipy.optimize import linear_sum_assignment




pd.options.display.width = 180
pd.set_option('display.max_columns', 30)
 

####################################################

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





def pv_table_noprint (  df, col1,   col2 ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ col1 ], columns=[ col2 ], fill_value=0, aggfunc=np.sum )
  pti = pt.astype(int)  
  return pti

  
	
def pdline ():
  print ('\n\n\n--------------------------------------------------------------------------\n\n', file=logfile)
 

########################################################################################

consensus_clusterings_pkl =  "consensus_clusterings.pkl"

logfile_txt = "compare cell clusters with ground truth.txt"


data_directory  = "Zhengmix4eq"
data_folder = r"D:/scRNA-seq/" + data_directory
data_path = Path ( data_folder )

cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 

trees_data_folder = cell_data_folder + "/2000 trees" 
trees_data_path = Path( trees_data_folder )




####  log output
logfile_dsn  =  trees_data_path /  logfile_txt
logfile = open ( logfile_dsn,'w')


# pickle input
consensus_clusterings_dsn = trees_data_path / consensus_clusterings_pkl

# csv input
cell_class_dsn = data_path / "sce_full_Zhengmix4eq_cell_class.csv"  
 
########################################################################################

df_consensus_clusterings = pd.read_pickle ( consensus_clusterings_dsn )    
plog ( logfile, '\n df_consensus_clusterings : \n\n', df_consensus_clusterings )  
 
df_cell_class = pd.read_csv ( cell_class_dsn, usecols=[1] ).rename ( columns={'x':'cell_class'} )
plog ( logfile, '\n\n df_cell_class: \n\n', df_cell_class )
 
 
 
df_compare = df_consensus_clusterings[[4]]
df_compare ['cell_class'] = df_cell_class['cell_class'].values
df_compare.rename ( columns={ 4:'consensus_clusters'}, inplace=True )
plog ( logfile, '\n\n df_compare: \n\n', df_compare )
 
pti = pv_table ( df_compare, 'consensus_clusters', 'cell_class' ) 
df_xtab = pv_table_noprint ( df_compare , 'consensus_clusters', 'cell_class' ) 
arr_xtab = df_xtab.values    

row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )
classified = arr_xtab [row_ind, col_ind].sum()
plog ( logfile,  'classifications agree: ', classified )
  
  
 
df_compare = df_consensus_clusterings[[5]]
df_compare ['cell_class'] = df_cell_class['cell_class'].values
df_compare.rename ( columns={ 5:'consensus_clusters'}, inplace=True )
plog ( logfile, '\n\n df_compare: \n\n', df_compare )
 
pti = pv_table ( df_compare, 'consensus_clusters', 'cell_class' ) 
df_xtab = pv_table_noprint ( df_compare , 'consensus_clusters', 'cell_class' ) 
arr_xtab = df_xtab.values    

row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )
classified = arr_xtab [row_ind, col_ind].sum() 
plog ( logfile,  'classifications agree: ', classified )


 
df_compare = df_consensus_clusterings[[4,5]]
df_compare.rename ( columns={ 4:'4'}, inplace=True )
df_compare.rename ( columns={ 5:'5'}, inplace=True )
plog ( logfile, '\n\n df_compare: \n\n', df_compare )
 
pti = pv_table ( df_compare, '4', '5' ) 
df_xtab = pv_table_noprint ( df_compare ,  '4', '5' ) 
arr_xtab = df_xtab.values    

row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )
classified = arr_xtab [row_ind, col_ind].sum() 
plog ( logfile,  'classifications agree: ', classified )

 
logfile.close()

 
################################################################################################################################ 
##                                                                                                                            ##
##  end program: compare_cell_clusters_with_ground_truth.py                                                                   ##
##                                                                                                                            ## 
################################################################################################################################ 
################################################################################################################################
