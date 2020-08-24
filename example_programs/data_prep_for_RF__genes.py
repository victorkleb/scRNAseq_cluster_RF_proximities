################################################################################################################################ 
##                                                                                                                            ##
##  start program: data_prep_for_RF__genes.py                                                                                 ##
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION data_prep_for_RF_genes
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

pctl_list = [.01,.05, .10, .25, .5, .75, .90, .95, .99 ]
 
########################################################################################

binomial_deviance_pkl = 'binomial deviance - genuine data.pkl'
standardized_null_residuals_pkl = 'standardized_null_residuals.pkl'
UMI_counts_csv =  "sce_full_Zhengmix4eq_counts.csv"


data_folder = r"D:/scRNA-seq/Zhengmix4eq"
data_path = Path ( data_folder )

gene_data_folder = data_folder + "/gene_data" 
gene_data_path = Path( gene_data_folder )

# Create folder if necessary	
if not os.path.exists(gene_data_path):
    os.mkdir(gene_data_path)

	

####  log output
logfile_dsn  =  data_path /  "data_prep_for_RF__genes.txt"
logfile = open ( logfile_dsn,'w')

#### pickle output - pandas data frame
standardized_null_residuals_dsn  = data_path / standardized_null_residuals_pkl 


# input - counts
UMI_counts_dsn = data_path / UMI_counts_csv

#### pickle input - pandas data frame
binomial_deviance_dsn  = data_path / binomial_deviance_pkl 


######################################################################################################################################

df_UMI_counts = pd.read_csv ( UMI_counts_dsn,  index_col = 0  )
plog ( logfile, '\n\n df_UMI_counts:\n\n' , df_UMI_counts )

 
df_binomial_deviance = pd.read_pickle ( binomial_deviance_dsn )
plog ( logfile, '\n\n df_binomial_deviance:\n' , df_binomial_deviance )

pdline( logfile )



result_list = data_prep_for_RF_genes  ( df_UMI_counts, df_binomial_deviance, 164, gene_data_path )
df_null_residuals_standardized = result_list[0]
plog ( logfile, '\n\n df_null_residuals_standardized:\n' , df_null_residuals_standardized )



df_null_residuals_standardized.to_pickle ( standardized_null_residuals_dsn )




logfile.close()



 
################################################################################################################################ 
##                                                                                                                            ##
##  end program: data_prep_for_RF__genes.py                                                                                   ##
##                                                                                                                            ## 
################################################################################################################################