
import pandas as pd
import numpy  as np


import sys
sys.path.append("D:/scRNA-seq/Spec_clust_RFproximities_scRNAseq")

from  utilities import *



import os


from sklearn.ensemble import RandomForestClassifier

from sklearn.cluster import KMeans
from sklearn import preprocessing



from sklearn.metrics.cluster import adjusted_rand_score

from scipy.stats import chi2_contingency

from scipy.optimize import linear_sum_assignment


from random import shuffle



from numpy import eye
from numpy import linalg as LA


import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap 
   


import seaborn as sns

#########################################################################################
#########################################################################################


def  compute_binomial_deviance__genuine_and_randomized_counts ( df_UMI_counts, n_randomizations=11, pctl_list =  [.25,.5,.75,.9,.95, .96, .97, .98, .99,.995,.999 ] ):
  
  df_genuine = compute_binomial_deviance ( df_UMI_counts, 'genuine' )

  UMI_counts_index = df_UMI_counts.index
  df_deviance_randomized_list = []

  for rep in range ( n_randomizations ):
    print  ( 'randomization: ', rep )
  
    df_UMI_counts_randomized_list= []

    for col in df_UMI_counts.columns.values.tolist() :
      df_col_list =  df_UMI_counts [[ col ]].values.tolist()
      np.random.shuffle(df_col_list)
      df_col_randomized = pd.DataFrame ( index=UMI_counts_index, data = df_col_list, columns = [ col ] )
  
      df_UMI_counts_randomized_list.append ( df_col_randomized )

    df_UMI_counts_randomized = pd.concat ( df_UMI_counts_randomized_list, axis=1, sort=True )

    df_deviance_randomized = compute_binomial_deviance ( df_UMI_counts_randomized, rep )  
  
    df_deviance_randomized_list.append ( df_deviance_randomized )
    df_deviance_randomized_multis = pd.concat ( df_deviance_randomized_list,  axis=1, sort=True )


  
  df_percentiles_randomized = df_deviance_randomized_multis.describe( percentiles=pctl_list )   
  sel_index_list = df_percentiles_randomized.transpose().drop ( columns = [ 'count', 'mean', 'std', 'min' ] ).columns.values.tolist()  
  
  ser_median_percentiles = df_percentiles_randomized.median( axis = 1 )  
  
  
  df_tuple_list = []  
  for percentile in sel_index_list :
    threshold = ser_median_percentiles.loc [ percentile ]
    sel_gene_count = np.sum ( ( df_genuine ['genuine' ] >=  threshold )  )
    df_tuple_list.append ( ( percentile, threshold, sel_gene_count ) )
  
  df_summary_table = pd.DataFrame ( data = df_tuple_list, columns= ['percentile rank', 'median percentile - randomized data', '# genes selected' ] ).sort_index ( ascending=False ).set_index ( ['percentile rank'] )
  
  return  [ df_genuine, df_deviance_randomized_multis, df_summary_table ]
    
  
  
  
def compute_binomial_deviance ( df_UMI_counts, column_name ):  

  arr_gen_counts = df_UMI_counts.values 
  mat_gen_counts = np.matrix ( arr_gen_counts)

  n_vector = mat_gen_counts.sum(axis=0)
  p_vector = mat_gen_counts.sum(axis=1)

  pi_hat_vector = p_vector/ p_vector.sum()
  n_pi_hat = np.array ( pi_hat_vector * n_vector )
  term1 = arr_gen_counts * np.log ( 1e-20 + np.nan_to_num( arr_gen_counts / n_pi_hat ) )

  arr_ni_m_gen_counts = np.array( n_vector ) - arr_gen_counts
  n_1_m_pi_hat =  np.array( n_vector ) - n_pi_hat
  term2 = arr_ni_m_gen_counts * np.log ( 1e-20 +  arr_ni_m_gen_counts / n_1_m_pi_hat )

  rhs = term1 + term2
  rhs_values = np.ravel ( rhs )

  df_deviance = pd.DataFrame( index = df_UMI_counts.index, data = rhs.sum( axis=1 ), columns=[column_name] ) 

  return df_deviance
   
#########################################################################################  

  
def data_prep_for_RF_genes  ( df_UMI_counts, df_binomial_deviance, number_selected_genes, output_data_path, number_of_synthetic_data_sets=100, NR_ds = 'null_residuals_std', NR_synth_ds_base='null_residuals_std_synth_' ):

  NR_dsn = NR_ds + '.csv'
  NR_path_dsn = output_data_path / NR_dsn  

  df_binomial_deviance['rank'] = df_binomial_deviance['genuine'].rank ( ascending=False )
  df_bd_select = df_binomial_deviance.loc [ df_binomial_deviance['rank'] <= number_selected_genes ]
  df_filtered_counts = df_UMI_counts.loc [ df_bd_select.index ] 

  df_null_residuals = compute_null_residuals ( df_filtered_counts )
  df_NR_mean_0 =  df_null_residuals.subtract ( df_null_residuals.mean(axis=1), axis=0 )
  df_null_residuals_standardized = df_NR_mean_0.divide ( df_NR_mean_0.std(axis=1), axis=0 )


  list_col_names = list( df_null_residuals_standardized.columns.values )


  for rep in range( number_of_synthetic_data_sets ):
    print ( 'creating synthetic data set: ', rep )
    df_NR_copy = df_null_residuals_standardized.copy()
	
    for col in list_col_names:
      NR_col_list= df_NR_copy [ col ].values.tolist()
      np.random.shuffle( NR_col_list )
      df_NR_copy [ col ] = NR_col_list
     
    NR_synth_dsn = NR_synth_ds_base + str ( rep) + ".csv"  
    NR_synth_path_dsn	= output_data_path / NR_synth_dsn
    df_NR_copy.to_csv ( NR_synth_path_dsn ) 
  
 
 
  df_null_residuals_standardized.to_csv ( NR_path_dsn )
  
  return  [ df_null_residuals_standardized ]

    
  
  

def compute_null_residuals ( df_UMI_counts):

  arr_gen_counts = df_UMI_counts.values 
  mat_gen_counts = np.matrix ( arr_gen_counts )

  n_vector = mat_gen_counts.sum(axis=0)
  p_vector = mat_gen_counts.sum(axis=1) 
  pi_hat_vector = p_vector/ p_vector.sum()
  n_pi_hat = np.array ( pi_hat_vector * n_vector )
  n_1_m_pi_hat =  np.array( n_vector ) - n_pi_hat

  term1 = arr_gen_counts * np.log ( 1e-20 + np.nan_to_num( arr_gen_counts / n_pi_hat ) )

  arr_ni_m_gen_counts = np.array( n_vector ) - arr_gen_counts
  term2 = arr_ni_m_gen_counts * np.log ( 1e-20 +  arr_ni_m_gen_counts / n_1_m_pi_hat )

  rhs0 = 2* ( term1 + term2 )
  sign_arg = arr_gen_counts - n_pi_hat 
  sign_rhs = np.sign ( sign_arg )
  arr_null_residual = sign_rhs * np.sqrt ( rhs0 )

  df_null_residuals = pd.DataFrame( index = df_UMI_counts.index, data = arr_null_residual, columns = df_UMI_counts.columns ) 

  return df_null_residuals

    
 
#########################################################################################   
 
# https://www.manongdao.com/article-550147.html
# for creating proximity array 
 
# https://community.intel.com/t5/Intel-Distribution-for-Python/parallel-random-forest-scikit-learn/td-p/1092793
# n_jobs= -1 
 
def RandomForestClassifier_for_proximity   ( n_trees, input_data_path, output_data_path, NR_ds, NR_synth_ds_base, prox_ds_base, first_synth_copy ): 
 
  NR_dsn = NR_ds + '.csv'
  GM_path_dsn = input_data_path / NR_dsn


  df_genuine = pd.read_csv ( GM_path_dsn, index_col=0 )
  n_rows = df_genuine.shape[0]
  class_list = [1]*n_rows + [0]* n_rows

  synth_copy = first_synth_copy    
  GM_synth_exists = True

  while ( GM_synth_exists ):  
    GM_synth_dsn = NR_synth_ds_base + str ( synth_copy ) + ".csv"  
    GM_synth_path_dsn	= input_data_path / GM_synth_dsn
  
  
    if   os.path.exists( GM_synth_path_dsn ):  	
      print ( 'random forest classification for synthetic copy: ', synth_copy )	

      prox_dsn = prox_ds_base + str ( synth_copy ) + ".csv"  
      prox_path_dsn	= output_data_path / prox_dsn
  
      df_synth = pd.read_csv ( GM_synth_path_dsn, index_col=0 )
      X =  np.concatenate( ( df_genuine.values, df_synth.values ) )
  
      model = RandomForestClassifier(  n_estimators=n_trees,  min_samples_leaf=2, random_state=0, n_jobs=-1 ) 
      model.fit( X, class_list)
  
      terminals = model.apply(X)
      del (X)  
      nTrees = terminals.shape[1]

      a = terminals[:n_rows,0]
      proximity_array = 1*np.equal.outer(a, a)

      for i in range(1, nTrees):
        a = terminals[:n_rows,i]
        proximity_array += 1*np.equal.outer(a, a)

      df_proximity = pd.DataFrame ( index = df_genuine.index, data = proximity_array, columns = df_genuine.index.values ) 
      df_proximity.to_csv ( prox_path_dsn )
	  
      synth_copy += 1
	  
    else:
      GM_synth_exists = False	 	

  return [ synth_copy ]	
	

	
	

def RandomForestClassifier_for_proximity_by_gene   ( n_trees, input_data_path, output_data_path, NR_ds = 'null_residuals_std', NR_synth_ds_base='null_residuals_std_synth_', prox_ds_base = 'gene_proximities_', first_synth_copy=0 ):

  result_list = RandomForestClassifier_for_proximity  ( n_trees, input_data_path, output_data_path, NR_ds, NR_synth_ds_base, prox_ds_base, first_synth_copy ) 

  return result_list  

  
  
  


def RandomForestClassifier_for_proximity_by_cell   ( n_trees, input_data_path, output_data_path, GM_ds = 'gene_means_tr', GM_synth_ds_base='gene_means_tr_synth_', prox_ds_base = 'cell_proximities_',  first_synth_copy=0 ):

  result_list = RandomForestClassifier_for_proximity  ( n_trees, input_data_path, output_data_path, GM_ds, GM_synth_ds_base, prox_ds_base, first_synth_copy ) 

  return result_list  
  
  
#########################################################################################   	
 
def create_cluster_array ( df_clusters ): 
  arr_cvalues = np.unique ( df_clusters.values )
  
  n_obs =  df_clusters.shape[0] 
  arr_cluster = np.zeros ( ( n_obs, n_obs ) ) 
  
  for cvalue in arr_cvalues:    
    df_xt = df_clusters.loc [ df_clusters['clusters'] == cvalue ]
    ind_list = df_xt.index.values.tolist()	   
 
 
    for i in  range ( 1, len(ind_list) ):
       for j in range ( i ) :
          row_ = ind_list [i]
          col_ =  ind_list [j]
 
          arr_cluster [ row_ ] [ col_ ] = 1
          arr_cluster [ col_ ] [ row_ ] = 1

  return arr_cluster


  
  
def cons_clusters_for_n_clusters_and_range_list ( df_all_clusterings,  n_clusters, range_list ) :
# consensus clusters from individual clusterings - specified by range list
  
 
  for individual_clustering  in range_list[0:1]:
    df_clustering = df_all_clusterings[[ n_clusters ]].loc [ df_all_clusterings['individual_clustering'] == individual_clustering ].rename ( columns={ n_clusters:'clusters'} ).reset_index ().drop ( columns=['index'] )
    consensus_array = create_cluster_array ( df_clustering )
    
  for individual_clustering  in range_list[1:]:
    df_clustering = df_all_clusterings[[ n_clusters ]].loc [ df_all_clusterings['individual_clustering'] == individual_clustering ].rename ( columns={ n_clusters:'clusters'} ).reset_index ().drop ( columns=['index'] )
    c_array = create_cluster_array ( df_clustering )
    consensus_array = consensus_array + c_array
	
	
  consensus_array =  consensus_array / len ( range_list )
  Lm = Lsym ( consensus_array )
  eig_vals, eig_vecs = np.linalg.eigh( Lm ) 
  del ( Lm ) 
  
  ev_first_n  = (np.fliplr (eig_vecs)) [:, : n_clusters ]
  del ( eig_vecs )  
   
  cluster_list = spectral_clusters_from_eigenvectors( n_clusters, ev_first_n )  	
  df_cons_clusters  = pd.DataFrame (  data=cluster_list, columns=[n_clusters] )    
  return df_cons_clusters , consensus_array 
  
  
  
  
def cons_clusters_for_n_clusters_and_cons_array ( n_clusters, consensus_array ) :  
#  consensus clusters from consensus array

  Lm = Lsym ( consensus_array )

  eig_vals, eig_vecs = np.linalg.eigh( Lm ) 
  del ( Lm ) 
  
  ev_first_n  = (np.fliplr (eig_vecs)) [:, :(1 + n_clusters) ]
  del ( eig_vecs )   
   
  cluster_list = spectral_clusters_from_eigenvectors( n_clusters, ev_first_n )  	
  df_cons_clusters  = pd.DataFrame (  data=cluster_list, columns=[n_clusters] )  
  return df_cons_clusters 
  
  

 


def Lsym ( a_parm ):
  a = np.copy( a_parm )
  np.fill_diagonal( a, 0 )
  
  col_sum = np.sum(a, axis=0)
  col_sum_pwr = np.sqrt ( col_sum )
  col_sum_factor = np.reciprocal ( col_sum_pwr)
  col_sum_diag = np.diag( col_sum_factor )
 
  L = np.matmul ( col_sum_diag, np.matmul ( a, col_sum_diag ) ) 
  del ( col_sum_diag, a )
        
  return L

  
  

  
def  spectral_clusters_from_eigenvectors( n_clusters, ev_first_n ) : 
  ev_sel = ev_first_n [:, : n_clusters]
  ev_row_norm  = np.asmatrix  ( preprocessing.normalize( ev_sel, norm='l2')  ) 
  del (   ev_sel )
  
  kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit( ev_row_norm )
  cluster = kmeans.predict ( ev_row_norm )  
  
  cluster_list = cluster.tolist()  
  return   cluster_list 
  

  

def   spectral_clusterings ( max_clusters, input_data_path, prox_ds_base = 'gene_proximities_'  ):

  print ( 'calculate clusterings - one per synthetic data set:' ) 
 
  df_all_clusterings_list = []

  
  individual_clustering = 0  
  proximity_exists = True

  while ( proximity_exists ):  
  
    proximities_csv = prox_ds_base + str ( individual_clustering) + ".csv"  
    proximities_csv_dsn = input_data_path / proximities_csv  
  
    if   os.path.exists( proximities_csv_dsn ):  	
      print ( 'clustering synthetic data set: ', individual_clustering )	
  
      df_proximities = pd.read_csv ( proximities_csv_dsn, index_col=0 )
      index_values = df_proximities.index.values  
      index_name = df_proximities.index.name 
      df_proximities.index.name = 'index'
  
      proximity_array = df_proximities.values
      n_trees = proximity_array[0,0]
      proximity_array = proximity_array / n_trees
	  
      if ( individual_clustering == 0 ):
        print ( 'for reference - number of trees in random forest: ', n_trees )	  
      	  
 
 
      clusters_synth_set = [] 
  
      Lm = Lsym ( proximity_array )
      del ( proximity_array )

      eig_vals, eig_vecs = np.linalg.eigh( Lm )
      del ( Lm )
 
    
      ev_first_n  = (np.fliplr (eig_vecs)) [:, :(1 + max_clusters) ]
      del ( eig_vecs )  
  
  
      for clusters in range ( 2, 1 + max_clusters ):  
        cluster_list = spectral_clusters_from_eigenvectors( clusters, ev_first_n )  	
        clusters_synth_set.append ( cluster_list )  
      cluster_dict = dict ( zip ( list(range( 2,1 + max_clusters ) ) ,  clusters_synth_set )  )
  
      df_clusters = pd.DataFrame.from_dict (  cluster_dict )
      df_clusters['index'] = index_values
      df_clusters.set_index ( ['index'], inplace=True ) 
      df_clusters ['individual_clustering']  = individual_clustering
    
 
      df_all_clusterings_list.append ( df_clusters )
      individual_clustering = individual_clustering + 1

    else:
      proximity_exists = False	  
	  
	  
  df_all_clusterings = pd.concat ( df_all_clusterings_list,sort=False )
  df_all_clusterings.index.names = [ index_name ] 
  
  return df_all_clusterings 
 

 
 

 
 
 
 
def  consensus_clusterings  ( df_individual_clusterings ): 

  print  ( '\n\n calculate consensus clusterings' ) 
 
  n_synth_copies_clusterings = 1 + df_individual_clusterings['individual_clustering'].max().astype(int)

  column_list = df_individual_clusterings.columns.values.tolist()
  column_list.remove ( 'individual_clustering' )
  max_clusters = max ( column_list )
  

  cluster_index_values = df_individual_clusterings.loc [ df_individual_clusterings['individual_clustering'] == 0 ] .index.values
  index_name = df_individual_clusterings.index.name 
  df_individual_clusterings.index.name = 'index'  
  
  print ( '\n\n maximum number of clusters =', max_clusters )
  
  df_consensus_clusterings_list = []
  df_consensus_clusterings_H0_list = []
  df_consensus_clusterings_H1_list = []
  df_consensus_clusterings_Q0_list = []
  df_consensus_clusterings_Q1_list = []
  df_consensus_clusterings_Q2_list = []
  df_consensus_clusterings_Q3_list = []

  n_synth_copies_clusterings_Q0 = np.floor ( .6 + .25 * n_synth_copies_clusterings ).astype(int)
  n_synth_copies_clusterings_Q1 = np.floor ( .6 + .5  * n_synth_copies_clusterings ).astype(int)
  n_synth_copies_clusterings_Q2 = np.floor ( .6 + .75 * n_synth_copies_clusterings ).astype(int)


  Q0_range_list = list ( range ( n_synth_copies_clusterings_Q0 ) ) 
  Q1_range_list = list ( range ( n_synth_copies_clusterings_Q0, n_synth_copies_clusterings_Q1 ) ) 
  Q2_range_list = list ( range ( n_synth_copies_clusterings_Q1, n_synth_copies_clusterings_Q2 ) ) 
  Q3_range_list = list ( range ( n_synth_copies_clusterings_Q2, n_synth_copies_clusterings  ) ) 
  print ( '\n clusterings used in Q0 consensus clusters: \n', Q0_range_list )
  print ( '\n clusterings used in Q1 consensus clusters: \n', Q1_range_list )
  print ( '\n clusterings used in Q2 consensus clusters: \n', Q2_range_list )
  print ( '\n clusterings used in Q3 consensus clusters: \n', Q3_range_list )  
  print ( '\n\n' )  
  
  
  len_H0_range_list  = n_synth_copies_clusterings_Q1
  len_H1_range_list = n_synth_copies_clusterings  - n_synth_copies_clusterings_Q1
 
  
  for n_clusters in  range ( 2, 1 + max_clusters ):  
    print ('calculating consensus clusters: '  , n_clusters)
 
    df_cons_clusters_Q0, consensus_array_Q0 = cons_clusters_for_n_clusters_and_range_list ( df_individual_clusterings, n_clusters, Q0_range_list ) 
    df_consensus_clusterings_Q0_list.append ( df_cons_clusters_Q0 ) 
  
    df_cons_clusters_Q1, consensus_array_Q1 = cons_clusters_for_n_clusters_and_range_list ( df_individual_clusterings, n_clusters, Q1_range_list ) 
    df_consensus_clusterings_Q1_list.append ( df_cons_clusters_Q1 ) 
 
    consensus_array_H0 = (  ( len ( Q0_range_list ) * consensus_array_Q0 ) + ( len ( Q1_range_list ) * consensus_array_Q1 ) ) / ( len ( Q0_range_list ) + len ( Q1_range_list ) )
    df_cons_clusters_H0 = cons_clusters_for_n_clusters_and_cons_array ( n_clusters, consensus_array_H0 )  
    df_consensus_clusterings_H0_list.append ( df_cons_clusters_H0 )  
    del ( consensus_array_Q0, consensus_array_Q1 )
 
 
    df_cons_clusters_Q2, consensus_array_Q2 = cons_clusters_for_n_clusters_and_range_list ( df_individual_clusterings, n_clusters, Q2_range_list ) 
    df_consensus_clusterings_Q2_list.append ( df_cons_clusters_Q2 ) 
  
    df_cons_clusters_Q3, consensus_array_Q3 = cons_clusters_for_n_clusters_and_range_list ( df_individual_clusterings, n_clusters, Q3_range_list ) 
    df_consensus_clusterings_Q3_list.append ( df_cons_clusters_Q3 ) 
 
    consensus_array_H1 = (  ( len ( Q2_range_list ) * consensus_array_Q2 ) + ( len ( Q3_range_list ) * consensus_array_Q3 ) ) / ( len ( Q2_range_list ) + len ( Q3_range_list ) )
    df_cons_clusters_H1 = cons_clusters_for_n_clusters_and_cons_array ( n_clusters, consensus_array_H1 )  
    df_consensus_clusterings_H1_list.append ( df_cons_clusters_H1 )  
    del ( consensus_array_Q2, consensus_array_Q3 )	
	

    consensus_array = ( len_H0_range_list * consensus_array_H0  + len_H1_range_list * consensus_array_H1 ) / n_synth_copies_clusterings
    df_cons_clusters = cons_clusters_for_n_clusters_and_cons_array ( n_clusters, consensus_array )  
    df_consensus_clusterings_list.append ( df_cons_clusters )  
    del ( consensus_array_H0, consensus_array_H1 )  
    
  
  df_consensus_clusterings_Q0 = pd.concat ( df_consensus_clusterings_Q0_list, axis=1 )
  df_consensus_clusterings_Q0['cluster_set'] = 0
  df_consensus_clusterings_Q1 = pd.concat ( df_consensus_clusterings_Q1_list, axis=1 )
  df_consensus_clusterings_Q1['cluster_set'] = 1
  df_consensus_clusterings_Q2 = pd.concat ( df_consensus_clusterings_Q2_list, axis=1 )
  df_consensus_clusterings_Q2['cluster_set'] = 2 
  df_consensus_clusterings_Q3 = pd.concat ( df_consensus_clusterings_Q3_list, axis=1 )
  df_consensus_clusterings_Q3['cluster_set'] = 3
    
  df_consensus_clusterings_H0 = pd.concat ( df_consensus_clusterings_H0_list, axis=1 )
  df_consensus_clusterings_H0['cluster_set'] = 0  
  df_consensus_clusterings_H1 = pd.concat ( df_consensus_clusterings_H1_list, axis=1 )
  df_consensus_clusterings_H1['cluster_set'] = 1
   
  df_consensus_clusterings = pd.concat ( df_consensus_clusterings_list, axis=1 )
 
 
 
  df_list = [ df_consensus_clusterings_Q0, df_consensus_clusterings_Q1, df_consensus_clusterings_Q2, df_consensus_clusterings_Q3, df_consensus_clusterings_H0, df_consensus_clusterings_H1, df_consensus_clusterings ]
  for df in df_list:
    df['index'] = cluster_index_values
    df.set_index ( ['index'], inplace=True )
    df.index.names = [ index_name ]
 
 
  df_consensus_clusterings_Q = pd.concat ( [ df_consensus_clusterings_Q0, df_consensus_clusterings_Q1, df_consensus_clusterings_Q2, df_consensus_clusterings_Q3 ] )
  df_consensus_clusterings_H = pd.concat ( [ df_consensus_clusterings_H0, df_consensus_clusterings_H1 ] )

  return [ df_consensus_clusterings, df_consensus_clusterings_H, df_consensus_clusterings_Q ]

 
 
 
def   clusters__individual_and_consensus  ( max_clusters, input_data_path, prox_ds_base  ):

  df_individual_clusterings =  spectral_clusterings ( max_clusters, input_data_path, prox_ds_base  )
  index_name = df_individual_clusterings.index.name  

  cc_result_list = consensus_clusterings  ( df_individual_clusterings )
  df_individual_clusterings.index.name = index_name # changed in function due to call by reference
  
  return  [ df_individual_clusterings ] + cc_result_list
 
 
 
 
def   gene_clusters__individual_and_consensus  ( max_clusters, input_data_path, prox_ds_base = 'gene_proximities_'  ):

  result_list = clusters__individual_and_consensus  ( max_clusters, input_data_path, prox_ds_base )
  
  return   result_list
  
 
 
 
def   cell_clusters__individual_and_consensus  ( max_clusters, input_data_path, prox_ds_base = 'cell_proximities_'  ):

  result_list = clusters__individual_and_consensus  ( max_clusters, input_data_path, prox_ds_base )
  
  return   result_list
     
  

#########################################################################################  
  
def ARI_and_Misclassification_Error_for_individual_clusterings ( df_individual_clusterings ):

# kluge to get number of rows as divisor in Misclassification Error calculation
  df_synth_0 = df_individual_clusterings.loc [ df_individual_clusterings['individual_clustering'] == 0 ] 
  n_rows =   df_synth_0.shape[0]
  item_index = df_synth_0.index  

 
  n_synth_copies_clusterings = 1 + np.int ( df_individual_clusterings['individual_clustering'].max() )
  
  column_list = df_individual_clusterings.columns.values.tolist()
  column_list.remove ( 'individual_clustering' )
  clust_maxp1 = 1 +  max ( column_list )
  cluster_range_list = list ( range ( 2, clust_maxp1 ) )
  
  
  ARI_tuple_list = []  
  ME_tuple_list = [] 
  
 
  
  print ( '\n\n calculating comparisons for ', n_synth_copies_clusterings, ' individual clusterings' )


  df_items_misclassified_total = pd.DataFrame ( index = item_index, data=0, columns=column_list )

  for synth_1 in range ( 1, n_synth_copies_clusterings ):
    df_synth_1 = df_individual_clusterings.loc [ df_individual_clusterings['individual_clustering'] == synth_1 ]
  
    for synth_2 in range ( 0, synth_1 ):
      print ( synth_1, '  and  ', synth_2 ) 	
      df_synth_2 = df_individual_clusterings.loc [ df_individual_clusterings['individual_clustering'] == synth_2 ]  

      ARI_row_list = [synth_1, synth_2]
      ME_row_list = [synth_1, synth_2]
      df_items_misclassified = 0 * df_items_misclassified_total
	  
      for n_clusters in cluster_range_list:
        df_1 = 	df_synth_1[[n_clusters]].rename ( columns={ n_clusters:'individual_clustering_1'} )
        df_2 = 	df_synth_2[[n_clusters]].rename ( columns={ n_clusters:'individual_clustering_2'} )	  
        ARI= adjusted_rand_score( df_1['individual_clustering_1'].values, df_2['individual_clustering_2'].values )
        ARI_row_list.append ( ARI )

		
        df_1_2 = pd.concat ( [ df_1, df_2 ], axis=1, sort=False )
        df_1_2['error'] = 1		
        df_xtab = pv_table_noprint ( df_1_2 , 'individual_clustering_1', 'individual_clustering_2' )
        arr_xtab = df_xtab.values   
 
        row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )
        total = arr_xtab [row_ind, col_ind].sum() 
		
        indices = tuple ( zip ( row_ind, col_ind ) )		
        for row, column in indices:
          df_1_2['error'].loc [ ( df_1_2['individual_clustering_1']==row ) & ( df_1_2['individual_clustering_2']==column ) ] = 0		
		  		  
        ME = 1 - total / n_rows 	  
        ME_row_list.append ( ME )	
        df_items_misclassified [ n_clusters ] = df_1_2['error']		
		
      ARI_tuple_list.append ( tuple ( ARI_row_list ) )	  		
      ME_tuple_list.append ( tuple ( ME_row_list ) )
      df_items_misclassified_total = df_items_misclassified_total + df_items_misclassified 

 
	  
  df_Adjusted_Rand_Index = pd.DataFrame( data = ARI_tuple_list, columns=['individual_clustering_1', 'individual_clustering_2']+cluster_range_list ).set_index ( ['individual_clustering_1', 'individual_clustering_2'] )	  
  df_Misclassification_Error = pd.DataFrame( data = ME_tuple_list, columns=['individual_clustering_1', 'individual_clustering_2']+cluster_range_list ).set_index ( ['individual_clustering_1', 'individual_clustering_2'] )	  
  df_item_misclassification_error = df_items_misclassified_total /  df_Misclassification_Error.shape[0]  

  return [ df_Adjusted_Rand_Index, df_Misclassification_Error, df_item_misclassification_error ]  


  

  
  

def ARI_and_Misclassification_Error_for_individual_clusterings_wrt_consensus_clusters ( df_individual_clusterings, df_consensus_clusterings ):



  n_rows =   df_consensus_clusterings.shape[0]

  print ( '\n\n calculating comparisons with consensus clusters for all individual clusterings:' )
 
  n_synth_copies_clusterings = 1 + np.int ( df_individual_clusterings['individual_clustering'].max() )
  
  column_list = df_consensus_clusterings.columns.values.tolist()
  clust_maxp1 = 1 +  max ( column_list )
  cluster_range_list = list ( range ( 2, clust_maxp1 ) )    
  
  item_index = df_consensus_clusterings.index  
  

 
  ARI_tuple_list = []  
  ME_tuple_list = []
  df_items_misclassified_total = pd.DataFrame ( index = item_index, data=0, columns=column_list )
  
  
  for synth_1 in range  ( n_synth_copies_clusterings ):
    print ( '   ', synth_1 )   
    df_synth_1 = df_individual_clusterings.loc [ df_individual_clusterings['individual_clustering'] == synth_1 ]
  
    ARI_row_list = [synth_1]
    ME_row_list = [synth_1]
    df_items_misclassified = 0 * df_items_misclassified_total	
	
    for n_clusters in cluster_range_list:
      df_1 = 	df_synth_1[[n_clusters]].rename ( columns={ n_clusters:'individual_clustering_1'} )
	
      df_CC = df_consensus_clusterings[[n_clusters]].rename ( columns={ n_clusters:'CC'} )	  
      ARI= adjusted_rand_score( df_1['individual_clustering_1'].values, df_CC['CC'].values )
      ARI_row_list.append ( ARI )
	  	  
      df_1_CC = pd.concat ( [ df_1, df_CC ], axis=1, sort=False )
      df_1_CC['error'] = 1	  
	  
      df_xtab = pv_table_noprint ( df_1_CC , 'individual_clustering_1', 'CC' )
      arr_xtab = df_xtab.values   
	  
      row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )
      total = arr_xtab [row_ind, col_ind].sum() 	  
	  
      indices = tuple ( zip ( row_ind, col_ind ) )		
      for row, column in indices:
        df_1_CC['error'].loc [ ( df_1_CC['individual_clustering_1']==row ) & ( df_1_CC['CC']==column ) ] = 0	  
 		  		  
      ME = 1 - total / n_rows 	  
      ME_row_list.append ( ME )	
      df_items_misclassified [ n_clusters ] = df_1_CC['error']		
	  
		
    ARI_tuple_list.append ( tuple ( ARI_row_list ) )	  		
    ME_tuple_list.append ( tuple ( ME_row_list ) )	 
    df_items_misclassified_total = df_items_misclassified_total + df_items_misclassified 	  	  	
	  
  df_Adjusted_Rand_Index_to_CC = pd.DataFrame( data = ARI_tuple_list, columns=['individual_clustering_1']+cluster_range_list ).set_index ( ['individual_clustering_1'] )	  
  df_Misclassification_Error_to_CC = pd.DataFrame( data = ME_tuple_list, columns=['individual_clustering_1']+cluster_range_list ).set_index ( ['individual_clustering_1'] )	 
  df_item_misclassification_error = df_items_misclassified_total /  df_Misclassification_Error_to_CC.shape[0]  
  
  return [ df_Adjusted_Rand_Index_to_CC, df_Misclassification_Error_to_CC, df_item_misclassification_error ]  

 
 
 

 
 
def ARI_and_Misclassification_Error_for_subset_consensus_clusterings ( df_consensus_clusterings_S ):


# kluge to get number of rows as divisor in Misclassification Error calculation
  df_split_0 = df_consensus_clusterings_S.loc [ df_consensus_clusterings_S['cluster_set'] == 0 ] 
  n_rows =   df_split_0.shape[0]
  
 
  n_splits = 1 + np.int ( df_consensus_clusterings_S['cluster_set'].max() )
  
  column_list = df_consensus_clusterings_S.columns.values.tolist()
  column_list.remove ( 'cluster_set')
  clust_maxp1 = 1 +  max ( column_list )
  cluster_range_list = list ( range ( 2, clust_maxp1 ) )
  
  
  ARI_tuple_list = []  
  ME_tuple_list = []

  
  
  print ( '\n\n calculating comparisons for ', n_splits, ' sets of consensus clusterings derived from disjoint sets of clusterings:' )

  
  for split_1 in range ( 1, n_splits ):
    df_split_1 = df_consensus_clusterings_S.loc [ df_consensus_clusterings_S['cluster_set'] == split_1 ]
  
    for split_2 in range ( 0, split_1 ):
      print ( split_1, '  and  ', split_2 ) 	
      df_split_2 = df_consensus_clusterings_S.loc [ df_consensus_clusterings_S['cluster_set'] == split_2 ]  

      ARI_row_list = [split_1, split_2]
      ME_row_list = [split_1, split_2]
	
      for n_clusters in cluster_range_list:
        df_1 = 	df_split_1[[n_clusters]].rename ( columns={ n_clusters:'split_1'} )
        df_2 = 	df_split_2[[n_clusters]].rename ( columns={ n_clusters:'split_2'} )	  
        ARI= adjusted_rand_score( df_1['split_1'].values, df_2['split_2'].values )
        ARI_row_list.append ( ARI )
	  	  
        df_1_2 = pd.concat ( [ df_1, df_2 ], axis=1, sort=False )
        df_xtab = pv_table_noprint ( df_1_2 , 'split_1', 'split_2' )
        arr_xtab = df_xtab.values     	
 
        row_ind, col_ind = linear_sum_assignment( arr_xtab, maximize=True )
        total = arr_xtab [row_ind, col_ind].sum() 
			  		  
        ME = 1 - total / n_rows 	  
        ME_row_list.append ( ME )				
		
      ARI_tuple_list.append ( tuple ( ARI_row_list ) )	  		
      ME_tuple_list.append ( tuple ( ME_row_list ) )	 
	 	  	
	  
  df_Adjusted_Rand_Index = pd.DataFrame( data = ARI_tuple_list, columns=['split_1', 'split_2']+cluster_range_list ).set_index ( ['split_1', 'split_2'] )	  
  df_Misclassification_Error = pd.DataFrame( data = ME_tuple_list, columns=['split_1', 'split_2']+cluster_range_list ).set_index ( ['split_1', 'split_2'] )	  

  return [ df_Adjusted_Rand_Index, df_Misclassification_Error ]  
  
  
  
  
  
 
#########################################################################################  
  
def  genes_mean_proximities    ( input_data_path, prox_ds_base = 'gene_proximities_'  ):
  
  result_list = mean_proximity ( input_data_path, prox_ds_base )
  
  return result_list  

  
  
def  cells_mean_proximities    ( input_data_path, prox_ds_base = 'cell_proximities_'  ):
  
  result_list = mean_proximity ( input_data_path, prox_ds_base )
  
  return result_list  

  
 
  
 

def   mean_proximity ( input_data_path, prox_ds_base ):

  print ( 'calculate mean of RF proximity matrices:' )  
  
  individual_clustering = 0  
  proximity_exists = True

  while ( proximity_exists ):  
  
    proximities_csv = prox_ds_base + str ( individual_clustering) + ".csv"  
    proximities_csv_dsn = input_data_path / proximities_csv  
  
    if   os.path.exists( proximities_csv_dsn ):  	
      print ( 'synthetic data set: ', individual_clustering )	 
      df_proximities = pd.read_csv ( proximities_csv_dsn, index_col=0 )
	  
      if ( individual_clustering == 0 ):
        df_mean_proximity = 0 * df_proximities

      df_mean_proximity = df_mean_proximity + df_proximities		     	  
 
      individual_clustering = individual_clustering + 1

    else:
      proximity_exists = False	  

  index_values = df_mean_proximity.index.values
  value_0_0 = df_mean_proximity[ index_values[0] ].loc [ index_values[0] ]
	  
  df_mean_proximity = df_mean_proximity / value_0_0
  
  return [ df_mean_proximity ]
 


#########################################################################################
   

def  within_distance_sum ( array_distance, cluster_array ):

  clusters = 1 + np.max ( cluster_array )

  D_list = [] 
 
  for i in range ( 0, clusters ):
    select_array  =   ( cluster_array == i )
    n = np.sum ( select_array )	
    cluster_distance_subarry = array_distance [ select_array ]	[:, select_array ]
    D = np.sum ( cluster_distance_subarry ) 
    D_list.append ( 0.5 * D / n )	
 
  return ( sum ( D_list )  )
 
 
  

def gap ( df_mean_proximity, df_consensus_clusterings, shuffles=100 ):

  print ( 'computing gap statistic' ) 

  proximity_array_mean = df_mean_proximity.values
  array_distance = ( 1.0 - 0.999999 * proximity_array_mean )
  del ( proximity_array_mean )
 

  cons_cluster_columns_array = df_consensus_clusterings.columns.values

  log_W0_list = []

  for clusters in  cons_cluster_columns_array.tolist() :
    print ('  number of clusters: ',   clusters  )
    cluster_array = np.array ( df_consensus_clusterings [ clusters ] )
    W = within_distance_sum ( array_distance, cluster_array )
    log_W0_list.append ( np.log ( W ) )
  
  df_log_W0 = pd.DataFrame.from_dict ( dict ( zip ( cons_cluster_columns_array, log_W0_list ) ), orient='index'  )


              
  print ( '\n\n shuffle for null distribution' )

  df_log_W_list = []
  cluster_row_list = df_consensus_clusterings.index.values.tolist()


  std_factor = ( 1+ (1/shuffles) )**.5

  for shuffle_count in range ( shuffles ): 
    print ('  shuffle_count: ',   shuffle_count  )

    shuffle ( cluster_row_list )
    df_shuffled_clusters = df_consensus_clusterings.loc [ cluster_row_list ]

    log_W_list = []

    for clusters in cons_cluster_columns_array.tolist() :
      cluster_array = np.array ( df_shuffled_clusters [ clusters ] )
      W = within_distance_sum ( array_distance, cluster_array )
      log_W_list.append ( np.log ( W ) )
    
    df_log_W = pd.DataFrame.from_dict ( dict ( zip ( df_consensus_clusterings.columns.values, log_W_list ) ), orient='index'  ).transpose()  
    df_log_W_list.append ( df_log_W )  

  df_log_W = pd.concat ( df_log_W_list )  
  df_log_W_mean = df_log_W.mean ( axis=0 ).to_frame( name = 'mean' )
  df_log_W_std = std_factor * df_log_W.std ( axis=0 ).to_frame( name = 'std' )


  df_gap = pd.concat ( [ df_log_W_mean, df_log_W_std, df_log_W0 ],  axis=1 )
  df_gap['gap'] = df_gap ['mean'] - df_gap [0]
  df_gap['gap - std'] = df_gap ['gap'] - df_gap ['std']
  df_gap['gap + std'] = df_gap ['gap'] + df_gap ['std']
  df_gap['gap - 2*std'] = df_gap ['gap'] - 2*df_gap ['std']
  df_gap['gap + 2*std'] = df_gap ['gap'] + 2*df_gap ['std']
  df_gap['lag gap - std'] = df_gap['gap - std'] .shift(-1)
  df_gap['delta'] = df_gap['gap'] - df_gap['lag gap - std']

  return [ df_gap, shuffles ]
  

  
  
#########################################################################################  
  
def data_prep_for_RF_cells  ( df_null_residuals_standardized, df_consensus_clusterings, number_selected_gene_clusters, output_data_path, number_of_synthetic_data_sets=100, GM_ds = 'gene_means_tr', GM_synth_ds_base='gene_means_tr_synth_' ):

  GM_dsn = GM_ds + '.csv'
  GM_path_dsn = output_data_path / GM_dsn  


  df_selected_clustering = df_consensus_clusterings[[ number_selected_gene_clusters ]].rename ( columns = { number_selected_gene_clusters:'cluster' } )

  df_null_residuals_ri = df_null_residuals_standardized.reset_index().rename ( columns = { 'index':'gene' } )
  df_null_residuals_melt = df_null_residuals_ri.melt ( id_vars=['gene'] ) 
  df_null_residuals_melt_add_cluster = df_null_residuals_melt.set_index ( ['gene']).merge ( df_selected_clustering, how='inner', left_index=True, right_index=True  )
  df_null_residuals_melt_add_cluster_ri = df_null_residuals_melt_add_cluster.reset_index()

  pt = pd.pivot_table( df_null_residuals_melt_add_cluster_ri, values='value',  index=[ 'cluster' ], columns=[ 'variable' ], aggfunc=np.mean )
  df_gene_cluster_means = pt[ df_null_residuals_standardized.columns.values.tolist() ]
  df_gene_cluster_means_tr = df_gene_cluster_means.transpose()
  df_gene_cluster_means_tr.index.name = 'cell'

  list_col_names = list( df_gene_cluster_means_tr.columns.values )


  for rep in  range( number_of_synthetic_data_sets ):
    print ( ' create synthetic data set: ', rep )
    df_gcm_tr_copy = df_gene_cluster_means_tr.copy()  
  
    for col in list_col_names:
      gcm_col_list= df_gcm_tr_copy [ col ].values.tolist()
      np.random.shuffle( gcm_col_list )
      df_gcm_tr_copy [ col ] = gcm_col_list 

    GM_synth_dsn = GM_synth_ds_base + str ( rep) + ".csv"  
    GM_synth_path_dsn	= output_data_path / GM_synth_dsn
    df_gcm_tr_copy.to_csv ( GM_synth_path_dsn ) 	

  df_gene_cluster_means_tr.to_csv ( GM_path_dsn  )	
	
  return  [ df_gene_cluster_means ]


#########################################################################################  
  
def proximity_array_means ( df_mean_proximities, df_clusterings, number_of_clusters ):
## diagonal [should be all=1] is excluded from numerator and denominator

  proximity_array = df_mean_proximities.values

  df_number_of_clusters = df_clusterings[ [ number_of_clusters ] ].rename ( columns={ number_of_clusters: 'cluster' } )
  df_number_of_clusters['indx'] = range ( proximity_array.shape[0] )

  pa_cluster_count = np.zeros(( number_of_clusters, number_of_clusters  ))
  pa_cluster_totals = np.zeros(( number_of_clusters, number_of_clusters  ))

  for   row in range ( number_of_clusters ):
    for column in range ( number_of_clusters ):
  
      iter_row = df_number_of_clusters['indx'].loc [ df_number_of_clusters[ 'cluster' ] == row ].values
      iter_column = df_number_of_clusters['indx'].loc [ df_number_of_clusters[ 'cluster' ] == column ].values	
      pa_number_of_clusters_row_col = proximity_array [ iter_row] [:, iter_column ]
      pa_cluster_totals [row, column ] =  np.sum ( pa_number_of_clusters_row_col )
      pa_cluster_count [row, column ] =  pa_number_of_clusters_row_col.shape[0] * pa_number_of_clusters_row_col.shape[1]
      if ( row == column ):
        block_size = 	( df_number_of_clusters[ 'cluster' ] == row  ).sum()
        pa_cluster_totals [row, column ] =  pa_cluster_totals [row, column ] - block_size
        pa_cluster_count [row, column ]  =   pa_cluster_count [row, column ] - block_size

  pa_cluster_mean = pa_cluster_totals / pa_cluster_count

  return pa_cluster_mean  

 
   
   

def reorder ( x_array ):

  col_sum = np.sum(x_array, axis=0)
  d_inv  = np.asmatrix ( np.diag ( 1./col_sum )  )

  s_mat = np.asmatrix ( x_array )
  I_DInvS = eye ( len(col_sum) ) - d_inv*s_mat
  w, v = LA.eig ( I_DInvS )
 

  arg_sort = np.argsort ( w )
  evec_p = np.asarray ( np.asarray(v)[:,arg_sort[1] ] )	
  cluster_permutation = np.argsort ( evec_p )	
  
  return cluster_permutation  

  
  

def plot_cluster_proximity_array ( plot_data, title2, plot_dsn ):

  fig = plt.figure()

  plt.title( 'NJW Consensus Cluster Proximity Array\n' + title2 )
  plt.imshow ( plot_data , cmap=plt.cm.Greys )
  ticks = list ( range ( plot_data.shape[0] )  )
  plt.xticks ( ticks )
  plt.yticks ( ticks )  
  plt.colorbar()

  plt.savefig( plot_dsn, dpi=300 )

   




# https://stackoverflow.com/questions/3279560/reverse-colormap-in-matplotlib

cdict1 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.1),
                   (1.0, 1.0, 1.0)),

         'blue': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'green':  ((0.0, 0.0, 1.0),
                   (0.5, 0.1, 0.0),
                   (1.0, 0.0, 0.0))
         }
		 
green_red1 = LinearSegmentedColormap('GreenRed1', cdict1)
  
   
   

 
# to remove ticks
# https://stackoverflow.com/questions/55506917/how-do-i-remove-the-axis-tick-marks-on-seaborn-heatmap

 
def data_prep_and_heat_maps  ( data_descriptor, number_of_gene_clusters, number_of_cell_clusters, df_null_residuals_standardized, df_gene_consensus_clusterings, \
df_cell_consensus_clusterings, df_gene_mean_proximities, df_cell_mean_proximities, heatmap_data_path ):
 


  # def HM_plot  ( df_heat, titl, plot_dsn  ):  #  before upgrade 2020 08 13
  
    # fig = plt.figure()  
    # sns.set(font_scale=0.50) 

    # g = sns.clustermap ( df_heat , xticklabels=blank_x_ticks, yticklabels=blank_y_ticks, cmap=green_red1 , row_cluster=False, col_cluster=False, 
           # row_colors=rbar_colors , col_colors=cbar_colors )
	   
    # plt.suptitle ( titl, fontsize=12 )
    # plt.savefig( plot_dsn, dpi=300 )
    # return


  def HM_plot  ( df_heat, titl, plot_dsn  ):
  
    fig = plt.figure()  
    sns.set(font_scale=0.50) 

    g = sns.clustermap ( df_heat , xticklabels=blank_x_ticks, yticklabels=blank_y_ticks, cmap=green_red1 , row_cluster=False, col_cluster=False, 
           row_colors=rbar_colors , col_colors=cbar_colors )
		   
    g.ax_heatmap.tick_params(tick2On=False, labelsize=False, length=0)

		   
    plt.suptitle ( titl, fontsize=12 )
    plt.savefig( plot_dsn, dpi=300 )
    return
     
 
 
  gene_mean_proximities_png = 'gene mean proximities.png'  
  gene_mean_proximities_reordered_png = 'gene mean proximities reordered.png'
  cell_mean_proximities_png = 'cell mean proximities.png'
  cell_mean_proximities_reordered_png = 'cell mean proximities reordered.png'
  heat_map_clip_1_png = 'heat map - clip-1.png' 
  heat_map_percentiles_png = 'heat map - percentiles.png'  

  gene_cluster_mean_proximities_dsn = heatmap_data_path / gene_mean_proximities_png  
  gene_cluster_mean_proximities_reordered_dsn = heatmap_data_path / gene_mean_proximities_reordered_png 
  cell_cluster_mean_proximities_dsn = heatmap_data_path / cell_mean_proximities_png 
  cell_cluster_mean_proximities_reordered_dsn = heatmap_data_path / cell_mean_proximities_reordered_png  
  heat_map_clip_1_dsn = heatmap_data_path / heat_map_clip_1_png
  heat_map_percentiles_dsn = heatmap_data_path / heat_map_percentiles_png   
 
 
  n_genes = df_null_residuals_standardized.shape[0]
  n_cells = df_null_residuals_standardized.shape[1]

  pa_cluster_genes_mean =  proximity_array_means ( df_gene_mean_proximities, df_gene_consensus_clusterings, number_of_gene_clusters )
  pa_cluster_cells_mean =  proximity_array_means ( df_cell_mean_proximities, df_cell_consensus_clusterings, number_of_cell_clusters )
   

  gene_permutation = reorder ( pa_cluster_genes_mean )
  g_reordered = np.take ( np.take ( pa_cluster_genes_mean, gene_permutation, axis=0 ), gene_permutation, axis=1 )
  plot_cluster_proximity_array ( pa_cluster_genes_mean, "gene clusters", gene_cluster_mean_proximities_dsn )
  plot_cluster_proximity_array ( g_reordered, "gene clusters - reordered with Ding's algorithm",  gene_cluster_mean_proximities_reordered_dsn )

  cell_permutation = reorder ( pa_cluster_cells_mean )
  g_reordered = np.take ( np.take ( pa_cluster_cells_mean, cell_permutation, axis=0 ), cell_permutation, axis=1 )
  plot_cluster_proximity_array ( pa_cluster_cells_mean, "cell clusters", cell_cluster_mean_proximities_dsn )
  plot_cluster_proximity_array ( g_reordered, "cell clusters - reordered with Ding's algorithm",  cell_cluster_mean_proximities_reordered_dsn )

  ##################
 
  df_cell_cluster_sel = df_cell_consensus_clusterings[ [ number_of_cell_clusters ] ].rename ( columns = { number_of_cell_clusters: 'cell cluster'} )
  df_gene_cluster_sel = df_gene_consensus_clusterings[ [ number_of_gene_clusters ] ].rename ( columns = { number_of_gene_clusters: 'gene cluster'} )

  df_null_residuals_standardized ['RF_gene_seq'] = range ( n_genes ) 
  df_null_residuals_add_gene_clusters = df_null_residuals_standardized.merge ( df_gene_cluster_sel, how='inner', left_index = True, right_index=True )

  df_gene_permutation = pd.DataFrame ( index = list ( gene_permutation ), data = list ( range ( number_of_gene_clusters ) ), columns = ['permuted gene cluster' ] )
  df_null_residuals_add_gene_permutation = df_null_residuals_add_gene_clusters.merge ( df_gene_permutation, how='inner', left_on='gene cluster',  right_index=True ).sort_values ( ['permuted gene cluster'] )
  df_null_residuals_add_gene_permutation [ 'permuted gene seq' ] = range ( n_genes ) 

  df_gene_cluster_order = df_null_residuals_add_gene_permutation[  [ 'RF_gene_seq', 'gene cluster', 'permuted gene cluster', 'permuted gene seq' ]  ] 
  df_null_residuals_gene_ordered = df_null_residuals_add_gene_permutation.drop ( columns= [ 'RF_gene_seq', 'gene cluster', 'permuted gene cluster', 'permuted gene seq' ] ) 

  df_cell_cluster_sel['RF_cell_seq'] = range ( n_cells )
  df_cell_permutation = pd.DataFrame ( index = list ( cell_permutation ), data = list ( range ( len ( cell_permutation ) ) ), columns = ['permuted cell cluster' ] )
  df_cell_cluster_order = df_cell_cluster_sel.merge ( df_cell_permutation, how='inner', left_on='cell cluster',  right_index=True ).sort_values (['permuted cell cluster'])  
  df_cell_cluster_order [ 'permuted cell seq' ] = range ( n_cells ) 


  df_null_residuals_gene_cell_ordered = df_null_residuals_gene_ordered [ df_cell_cluster_order.index.values.tolist() ]

  ##################

  df_null_residuals_gene_cell_ordered.index.name = None
  df_cells = df_cell_cluster_order[['permuted cell cluster','permuted cell seq']]

  dfr_ = df_gene_cluster_order  [ ['permuted gene cluster' ] ]
  gene_axis = list ( dfr_[ 'permuted gene cluster' ] ) 
 


  blank_y_ticks = len ( gene_axis ) * [' '] 

  arr_gene_clust_numbers = dfr_ ['permuted gene cluster' ].unique()
  df_cluster_ = pd.DataFrame ( index= arr_gene_clust_numbers, data =  arr_gene_clust_numbers % 2 )
  df_color_default = pd.DataFrame ( index=[0,1], data = [ 'white', 'grey' ] )
  df_gene_cluster_colors = df_cluster_.merge( df_color_default, how='inner', left_on=[0],right_index=True )
  df_gene_cluster_colors.drop ( [ 'key_0', '0_x' ],axis=1, inplace=True )
  df_gene_cluster_colors.rename ( columns = {'0_y':'color'}, inplace=True )

  dfr_colors = dfr_.merge ( df_gene_cluster_colors, how= 'left', left_on ='permuted gene cluster', right_index=True )
  rbar_colors = dfr_colors ['color'].rename("") 
  rbar_colors.index.name = None
 
  cell_axis = df_cells.index.tolist() 
  arr_cell_clust_numbers = df_cells ['permuted cell cluster' ].unique()
  df_cluster_ = pd.DataFrame ( index= arr_cell_clust_numbers, data =  arr_cell_clust_numbers % 2 )
  df_color_default = pd.DataFrame ( index=[0,1], data = [ 'white', 'grey' ] )
  df_cell_cluster_colors = df_cluster_.merge( df_color_default, how='inner', left_on=[0],right_index=True )
  df_cell_cluster_colors.drop ( [ 'key_0', '0_x' ],axis=1, inplace=True )
  df_cell_cluster_colors.rename ( columns = {'0_y':'color'}, inplace=True )

  dfc_colors0 = df_cells.merge ( df_cell_cluster_colors, how= 'left', left_on ='permuted cell cluster', right_index=True )
  cbar_colors = dfc_colors0 ['color'].rename("") 
  cbar_colors.index.name = None
    
  blank_x_ticks = len ( cell_axis ) * [' '] 


  title0 = 'scRNA-seq data:  heat map of clustered binomial deviance null residuals \n spectral clusters derived from random forest proximities \n'
  title1 = data_descriptor + ' data,  filtered with binomial deviance to select  ' + str ( n_genes ) + ' genes  -- grouped in ' + str (number_of_gene_clusters ) + ' clusters \n'
  title2 = str ( n_cells ) + ' cells  -- grouped in ' + str (number_of_cell_clusters) + ' clusters \n'
  title3 = 'gene-standardized null residuals'

  titl = title0 + title1 + title2 + title3
 
  HM_plot ( df_null_residuals_gene_cell_ordered.clip( -1.0 , 1.0 ), titl, heat_map_clip_1_dsn )


  title3 = 'gene percentiles'
  titl = title0 + title1 + title2 + title3

  df_heat_row_pct = df_null_residuals_gene_cell_ordered .rank ( axis=1, pct=True )
  HM_plot ( df_heat_row_pct, titl, heat_map_percentiles_dsn )
 

  pa_cluster_mean_list =  [ pa_cluster_genes_mean, pa_cluster_cells_mean ]  
  permutation_list = [ gene_permutation, cell_permutation ]
  dataframe_list = [ df_gene_cluster_order, df_cell_cluster_order, df_null_residuals_gene_cell_ordered ]
  
  return [ pa_cluster_mean_list, permutation_list,  dataframe_list ]  
 
   
  

#########################################################################################  
     
def Laplacian_scores ( df_proximity_array_mean, df_null_residuals_standardized ):

  proximities = df_proximity_array_mean.values 
  null_residuals_mat = np.asmatrix ( df_null_residuals_standardized.values ) 

  col_sum = np.sum( proximities, axis=0)
  D_mat = np.matrix ( np.diag( col_sum ) )
  L_mat = D_mat - np.asmatrix ( proximities )
  one_vector = np.ones ( ( L_mat.shape[0], 1 ) )

  f_T_D_1 = null_residuals_mat * D_mat * one_vector
  one_T_D_1 = np.transpose ( one_vector ) * D_mat * one_vector
  to_subtract = ( one_vector * np.transpose ( f_T_D_1 ) ) / one_T_D_1
  f_tilde =  np.transpose( null_residuals_mat ) - to_subtract

  LS_num = np.transpose ( f_tilde ) * L_mat * f_tilde 
  LS_denom = np.transpose ( f_tilde ) * D_mat * f_tilde 
  LS = np.diagonal ( LS_num ) / np.diagonal ( LS_denom )

  df_LS = pd.DataFrame ( index =  df_null_residuals_standardized.index, data=LS, columns=[ 'Laplacian_score' ] )
  
  return df_LS     
  
    
  
  
  
  