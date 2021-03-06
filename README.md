## scRNAseq_cluster_RF_proximities

This set of programs clusters single cell RNA-Seq data, represented by UMI counts, with no metric assumptions and no logarithmic transformation to treat zero counts.  
<br />
<br />
The approach uses
- binomial deviance to rank genes for filtering
- random forest classification to produce proximities to cluster genes and cells
- spectral consensus clustering using the random forest proximities
- distributions of adjusted Rand index or Misclassification Error distance to determine an appropriate number of clusters
- Laplacian scores computed with random forest cell proximities to confirm the gene ranks and provide additional - validation of the methodology.
<br />

For details, please refer to 
- example_programs.docx
- the folder example_programs: a 17-program stream to analyze the Zhengmix4eq data set
- functions that perform the analysis, and documentation:
  - FUNCTIONS_Spec_clust_RFproximities_scRNAseq.py
  - utilities.py
  - function documentation.docx
<br />

Update September 22, 2020
<br />
Programs to analyze outliers have been added, along with a document including analytic results and program descriptions:
- folder: example_programs_outlier_analysis
- outliers.docx