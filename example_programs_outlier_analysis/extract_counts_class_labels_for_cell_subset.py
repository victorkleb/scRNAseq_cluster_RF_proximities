

################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program: extract_counts_class_labels_for_cell_subset.py                                                             ##
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
#########################################################################################  

logfile_txt = "extract_counts_class_labels_for_cell_subset.txt"
cell_outliers_pkl = "cell_cluster_8_outliers_MME_0p1.pkl"


data_directory  = "Zhengmix8eq"
data_folder = r"D:/scRNA-seq/" + data_directory
data_path = Path ( data_folder)

output_data_directory  = "Zhengmix8eq_2877"
output_data_folder = r"D:/scRNA-seq/" + output_data_directory
output_data_path = Path ( output_data_folder)


cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 
cell_data_path = Path( cell_data_folder )

trees_data_folder =  cell_data_folder +  "/2000 trees" 
trees_data_path = Path ( trees_data_folder )

 
 
####  log output
logfile_dsn  =  data_path / logfile_txt 
logfile = open ( logfile_dsn,'w')


#### csv outputs
UMI_counts_out_dsn= output_data_path / "sce_full_Zhengmix8eq_counts.csv"
cell_class_out_dsn= output_data_path / "sce_full_Zhengmix8eq_cell_class.csv"  


### CSV inputs
UMI_counts_in_dsn= data_path / "sce_full_Zhengmix8eq_counts.csv"
cell_class_in_dsn= data_path / "sce_full_Zhengmix8eq_cell_class.csv"  

#### pickle input
cell_outliers_dsn = trees_data_path / cell_outliers_pkl

######################################################################################################################################

df_UMI_counts = pd.read_csv ( UMI_counts_in_dsn,  index_col = 0  )
plog ( logfile, '\n\n df_UMI_counts:\n\n' , df_UMI_counts )

df_cell_class = pd.read_csv ( cell_class_in_dsn, usecols=[1] ).rename ( columns={'x':'cell_class'} )
plog ( logfile, '\n\n df_cell_class: \n\n', df_cell_class )


df_cell_select_outliers = pd.read_pickle ( cell_outliers_dsn )
plog ( logfile, '\n\n\n df_cell_select_outliers : \n\n', df_cell_select_outliers )

df_cell_select_outliers['outliers'] = True
plog ( logfile, '\n\n\n df_cell_select_outliers : \n\n', df_cell_select_outliers )
pdline ( logfile )


df_cell = pd.DataFrame ( data = df_UMI_counts.columns.values, columns = ['cell'] )
df_cell_append_outliers = df_cell.merge ( df_cell_select_outliers, how='left', left_on = 'cell', right_index=True )
df_cell_append_outliers ['outliers'] = df_cell_append_outliers ['outliers'].fillna( False )
 
df_cell_drop_outliers = df_cell_append_outliers.loc [ ~ df_cell_append_outliers ['outliers'] ]. drop ( columns = ['outliers'] ).sort_index()
plog ( logfile, '\n\n\n df_cell_drop_outliers : \n\n', df_cell_drop_outliers )

list_cell_drop_outliers = df_cell_drop_outliers['cell'].values.tolist()

df_UMI_counts_out = df_UMI_counts [ list_cell_drop_outliers ]
plog ( logfile, '\n\n df_UMI_counts_out:\n\n' , df_UMI_counts_out )


df_cell_class_append_outliers = df_cell_class.merge ( df_cell_append_outliers, how='inner', left_index=True, right_index=True )
df_cell_class_drop_outliers = df_cell_class_append_outliers.loc [ ~ df_cell_class_append_outliers ['outliers'] ]. drop ( columns = ['outliers'] ).sort_index()
df_cell_class_drop_outliers['index_new'] = range ( df_cell_class_drop_outliers.shape[0] )
df_cell_class_out = df_cell_class_drop_outliers.set_index( 'index_new' ).rename ( columns={'cell_class':'x'} )
plog ( logfile, '\n\n\n df_cell_class_out : \n\n', df_cell_class_out )



df_UMI_counts_out.to_csv ( UMI_counts_out_dsn )
df_cell_class_out.to_csv ( cell_class_out_dsn )
  



logfile.close()


################################################################################################################################ 
##                                                                                                                            ##
##  end program: extract_counts_class_labels_for_cell_subset.py                                                               ##
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################