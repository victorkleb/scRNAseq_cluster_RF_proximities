################################################################################################################################ 
##                                                                                                                            ##
##  start program: compute_binomial_deviance.py                                                                               ##
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION compute_binomial_deviance__genuine_and_randomized_counts
#
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

genuine_data_pkl = 'binomial deviance - genuine data.pkl'
synthetic_data_pkl = 'binomial deviance - randomized data.pkl'
summary_pkl = 'binomial deviance summary table.pkl' 

data_folder = Path ( r"D:/scRNA-seq/Zhengmix4eq" )
 
 
####  log output
logfile_dsn  =  data_folder /  "compute_binomial_deviance.txt"
logfile = open ( logfile_dsn,'w')


#### pickle outputs - pandas data frames
genuine_out_dsn  = data_folder / genuine_data_pkl 
synthetic_out_dsn  = data_folder / synthetic_data_pkl 
summary_out_dsn = data_folder / summary_pkl 



# input
UMI_counts_dsn = data_folder / "sce_full_Zhengmix4eq_counts.csv"

######################################################################################################################################

df_UMI_counts = pd.read_csv ( UMI_counts_dsn,  index_col = 0  )
plog ( logfile, '\n\n df_UMI_counts:\n\n' , df_UMI_counts )
pdline( logfile )



result_list = compute_binomial_deviance__genuine_and_randomized_counts ( df_UMI_counts )

df_binomial_deviance__genuine_data = result_list[0]
df_binomial_deviance__randomized_data = result_list[1]
df_binomial_deviance_summary_table = result_list[2]

plog ( logfile, '\n\n\n df_binomial_deviance__genuine_data: \n', df_binomial_deviance__genuine_data )
plog ( logfile, '\n\n\n df_binomial_deviance__randomized_data: \n', df_binomial_deviance__randomized_data )
plog ( logfile, '\n\n\n df_binomial_deviance_summary_table: \n', df_binomial_deviance_summary_table )



df_binomial_deviance__genuine_data.to_pickle ( genuine_out_dsn )
df_binomial_deviance__randomized_data.to_pickle ( synthetic_out_dsn )
df_binomial_deviance_summary_table.to_pickle ( summary_out_dsn )


logfile.close()
 
################################################################################################################################ 
##                                                                                                                            ##
##  end program: compute_binomial_deviance.py                                                                                 ##
##                                                                                                                            ## 
################################################################################################################################ 