################################################################################################################################ 
##                                                                                                                            ##
##  start program: extract_data_from_package_DuoClustering2018.r                                                              ##
##                                                                                                                            ## 
################################################################################################################################


library(DuoClustering2018)


csv_counts     <-  "D:/scRNA-seq/Zhengmix4eq/sce_full_Zhengmix4eq_counts.csv"
csv_cell_class <-  "D:/scRNA-seq/Zhengmix4eq/sce_full_Zhengmix4eq_cell_class.csv"

sce<-sce_full_Zhengmix4eq()

gcounts <- counts(sce)
cell_class <- sce$phenoid

write.csv( gcounts, file = csv_counts )
write.csv( cell_class, file = csv_cell_class )

################################################################################################################################ 
##                                                                                                                            ##
##  end program: extract data from package DuoClustering2018.r                                                                ##
##                                                                                                                            ## 
################################################################################################################################