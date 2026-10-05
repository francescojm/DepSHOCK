library(data.table)
library(edgeR)
library(limma)

source('DepSHOCK_DE_analysis_library.R')

DE<-DS_run_experiment_DE('M14-BRAF',
                             path_location = '/Volumes/iorio/raw_data_inHouse/DepSHOCK/bulkRNAseq/',
                             path_file = '../DepSHOCK_data/raw_data_paths.xlsx',
                             output_path = '../DepSHOCK_data/DiffExp/',overwrite = TRUE)




DS_plot_DE_PCA(DE,padj_thr = 0.05)
