# scRNAseq_cluster_RF_proximities

This set of programs clusters single cell RNA-Seq data, represented by UMI counts, with no metric assumptions and no logarithmic transformation to treat zero counts.  



The approach uses
- binomial deviance to rank genes for filtering
- random forest classification to produce proximities to cluster genes and cells
- spectral consensus clustering using the random forest proximities
- distributions of adjusted Rand index or Misclassification Error distance to determine an appropriate number of clusters
- Laplacian scores computed with random forest cell proximities to confirm the gene ranks and provide additional - validation of the methodology.



For details, please refer to 
- overview.doc 
- example_programs: folder containing a set of programs to analyze the Zhengmix4eq data set
- functions that perform the analysis, and documentation:
  - FUNCTIONS_Spec_clust_RFproximities_scRNAseq.py
  - function documentation.docx
