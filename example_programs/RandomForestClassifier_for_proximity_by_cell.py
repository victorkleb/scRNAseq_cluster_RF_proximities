
################################################################################################################################
################################################################################################################################ 
##                                                                                                                            ##
##  start program: RandomForestClassifier_for_proximity_by_cell.py                                                            ##                                                                                
##                                                                                                                            ## 
################################################################################################################################
#
#  invoke FUNCTION  RandomForestClassifier_for_proximity_by_cell
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

data_directory  = "Zhengmix4eq"
data_folder = r"D:/scRNA-seq/" + data_directory
data_path = Path ( data_folder )

cell_data_folder = data_folder + "/cell_data - 5 gene clusters" 
cell_data_path = Path( cell_data_folder )


trees_data_folder =  cell_data_folder +  "/2000 trees" 
trees_data_path = Path ( trees_data_folder )

# Create folder if necessary	 
if not os.path.exists( trees_data_path ):
    os.mkdir( trees_data_path )

#############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

 
 
result_list = RandomForestClassifier_for_proximity_by_cell  ( 2000, cell_data_path, trees_data_path )

print ( '\n\n result_list: number of proximity arrays calculated: ', result_list )
	
	
################################################################################################################################ 
##                                                                                                                            ##
##  end program: RandomForestClassifier_for_proximity_by_cell.py                                                              ##                                                                                 
##                                                                                                                            ## 
################################################################################################################################
################################################################################################################################

