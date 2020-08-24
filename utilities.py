

import pandas as pd
import numpy  as np


def desc_list ( data_list, pctl_list = [.01,.05, .10, .25, .5, .75, .90, .95, .99 ] ):
  df_data = pd.DataFrame( data_list )
  pd_desc = df_data.describe( percentiles=  pctl_list )
  return pd_desc


  
  
def plog ( logfile, *args ): 
  for object in args: 
    print ( object, file = logfile, end='' )  
  
  
  
def pdline ( logfile ):
  plog ( logfile, '\n\n\n--------------------------------------------------------------------------\n\n' )




def pv_table_noprint (  df, col1,   col2 ):
  df_copy = df.copy()
  df_copy ['count'] = 1
  pt = pd.pivot_table( df_copy, values='count',  index=[ col1 ], columns=[ col2 ], fill_value=0, aggfunc=np.sum )
  pti = pt.astype(int)  
  return pti
    
