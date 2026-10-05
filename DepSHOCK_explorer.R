
library(dplyr)
library(ggplot2)

source('DepSHOCK_library.R')

DS_data<-DS_loading_and_formatting()
X<-DS_create_analysis_matrix(DS_data$lfc,DS_data$padj)

PCA<-DS_pca(X)


par(mfrow=c(2,2))
DS_plot_pca(PCA,meta = DS_data$metadata,pc1 = 1,pc2 = 2,color_by = 'oncogene_cellline')
DS_plot_pca(PCA,meta = DS_data$metadata,pc1 = 2,pc2 = 3,color_by = 'oncogene_cellline')
DS_plot_pca(PCA,meta = DS_data$metadata,pc1 = 1,pc2 = 3,color_by = 'oncogene_cellline')
DS_plot_pca(legend_only = TRUE,meta = meta,color_by = 'oncogene_cellline')


par(mfrow=c(2,2))
DS_plot_pca(PCA,meta = DS_data$metadata,pc1 = 1,pc2 = 2,color_by = 'time')
DS_plot_pca(PCA,meta = DS_data$metadata,pc1 = 2,pc2 = 3,color_by = 'time')
DS_plot_pca(PCA,meta = DS_data$metadata,pc1 = 1,pc2 = 3,color_by = 'time')
DS_plot_pca(legend_only = TRUE,meta = meta,color_by = 'time')

par(mfrow=c(2,2))
DS_plot_pca(PCA,meta = DS_data$metadata,pc1 = 1,pc2 = 2)
DS_plot_pca(PCA,meta = DS_data$metadata,pc1 = 2,pc2 = 3)
DS_plot_pca(PCA,meta = DS_data$metadata,pc1 = 1,pc2 = 3)
DS_plot_pca(legend_only = TRUE)

DS_plot_pca_by_time(
  pca = pca,
  meta = meta,
  color_by = "oncogene_cellline"
)

DS_plot_pca_by_time(
  pca = pca,
  meta = meta,
  color_by = "perturbation"
)

DS_plot_pca_by_cellline_perturbation(
  pca = pca,
  meta = meta,
  cell_lines = "M14"
)

DSB_obj <- DS_compute_DSB_component(
  lfc,
  meta
)

DSB_mat <- DSB_obj$DSB
DSB_meta <- DSB_obj$metadata



DSB_similarity <- DS_compare_DSB_similarity(
  lfc,
  meta,
  DSB_obj
)

aggregate(
  cbind(
    raw_AAVS1_cosine,
    DSB_cosine
  ) ~ time_h,
  data = DSB_similarity,
  FUN = mean
)

common_obj <- DS_compute_common_trace(
  lfc,
  meta,
  exclude_targets = "NT",
  consensus = "median"
)

common_NT_similarity <-
  DS_validate_common_trace_with_NT(
    common_obj,
    lfc,
    meta
  )

aggregate(
  cosine_common_NT ~ time_h,
  data = common_NT_similarity,
  FUN = mean
)

common_similarity <-
  DS_compare_common_trace(
    common_obj
  )

aggregate(
  cosine ~ time_h,
  data = common_similarity,
  FUN = mean
)

