
DS_loading_and_formatting <- function(
    path = '../DepSHOCK_data/depshock_dataset.rds',
    max_na_lfc = 1,
    max_na_padj = 1
) {
  
  DepSHOCK_obj <- readRDS(path)
  
  lfc  <- as.matrix(DepSHOCK_obj$logFC)
  padj <- as.matrix(DepSHOCK_obj$padj)
  meta <- DepSHOCK_obj$metadata
  
  storage.mode(lfc)  <- "numeric"
  storage.mode(padj) <- "numeric"
  
  ## genes -> rows, experiments -> columns
  lfc  <- t(lfc)
  padj <- t(padj)
  
  
  ## ========================================================
  ## Convert Ensembl IDs to gene symbols
  ## ========================================================
  
  ensembl_ids <- sub(
    "\\..*$",
    "",
    rownames(lfc)
  )
  
  gene_symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  gene_symbols <- unname(gene_symbols)
  
  
  ## Remove unmapped genes
  keep <- !is.na(gene_symbols) &
    gene_symbols != ""
  
  lfc  <- lfc[keep, , drop = FALSE]
  padj <- padj[keep, , drop = FALSE]
  gene_symbols <- gene_symbols[keep]
  
  
  ## Remove duplicated symbols
  unique_gene <- !duplicated(gene_symbols)
  
  lfc  <- lfc[unique_gene, , drop = FALSE]
  padj <- padj[unique_gene, , drop = FALSE]
  gene_symbols <- gene_symbols[unique_gene]
  
  rownames(lfc)  <- gene_symbols
  rownames(padj) <- gene_symbols
  
  
  ## ========================================================
  ## Filter genes according to NA fraction
  ## ========================================================
  
  na_fraction_lfc <- rowMeans(is.na(lfc))
  na_fraction_padj <- rowMeans(is.na(padj))
  
  keep_na <- (
    na_fraction_lfc <= max_na_lfc &
      na_fraction_padj <= max_na_padj
  )
  
  cat(
    "Genes before NA filtering:",
    nrow(lfc),
    "\n"
  )
  
  cat(
    "Genes removed because of NA threshold:",
    sum(!keep_na),
    "\n"
  )
  
  lfc  <- lfc[keep_na, , drop = FALSE]
  padj <- padj[keep_na, , drop = FALSE]
  
  cat(
    "Genes retained:",
    nrow(lfc),
    "\n"
  )
  
  ## Final NA percentages
  na_pct_lfc  <- 100 * mean(is.na(lfc))
  na_pct_padj <- 100 * mean(is.na(padj))
  
  cat(
    "Final NA percentage in lfc:",
    round(na_pct_lfc, 3),
    "%\n"
  )
  
  cat(
    "Final NA percentage in padj:",
    round(na_pct_padj, 3),
    "%\n"
  )
  
  
  ## ========================================================
  ## Align metadata
  ## ========================================================
  
  meta <- meta[
    match(colnames(lfc), meta$guide_id),
    ,
    drop = FALSE
  ]
  
  
  ## ========================================================
  ## Sanity checks
  ## ========================================================
  
  stopifnot(
    identical(meta$guide_id, colnames(lfc)),
    identical(rownames(lfc), rownames(padj)),
    identical(colnames(lfc), colnames(padj))
  )
  
  
  return(
    list(
      lfc = lfc,
      padj = padj,
      metadata = meta
    )
  )
}

DS_load_gex_dep_data_matrices<-function(cancerType,centering=TRUE){
  
  load('/Volumes/Macintosh HD/Users/francesco.iorio/Dropbox/CODING/vus/VUS 2025/data/Robj/basa_exp_FPKM.RData')

  CMP_annot <- read_csv('/Volumes/Macintosh HD/Users/francesco.iorio/Dropbox/CODING/vus/VUS 2025/data/raw/model_list_20241120.csv')
  
  duplicatedGenes<-rownames(basal_exp)[which(duplicated(rownames(basal_exp)))]
  
  for (i in duplicatedGenes){
    ids<-which(rownames(basal_exp)==i)
    basal_exp<-basal_exp[-ids[order(rowSums(is.na(basal_exp[ids,])),decreasing=TRUE)],]
  }
  
  cls<-intersect(CMP_annot$model_name[which(CMP_annot$cancer_type==cancerType)],colnames(basal_exp))
  basal_exp<-basal_exp[,cls]
  
  scaled_depFC<-read.csv('/Volumes/Macintosh HD/Users/francesco.iorio/Dropbox/CODING/vus/VUS 2025/data/raw/CRISPRGeneEffect.csv',row.names = 1)
  colnames(scaled_depFC)<-gsub("\\..*","",colnames(scaled_depFC))
  
  toremove<-which(is.na(CMP_annot$model_id[match(rownames(scaled_depFC),CMP_annot$BROAD_ID)]))
  scaled_depFC<-scaled_depFC[-toremove,]
  rownames(scaled_depFC)<-CMP_annot$model_name[match(rownames(scaled_depFC),CMP_annot$BROAD_ID)]
  
  scaled_depFC<-t(scaled_depFC)
  
  duplicatedGenes<-rownames(scaled_depFC)[which(duplicated(rownames(scaled_depFC)))]
  
  while(length(duplicatedGenes)>0){
      ids<-which(rownames(scaled_depFC)==duplicatedGenes[1])
      scaled_depFC<-scaled_depFC[-ids[order(rowSums(is.na(scaled_depFC[ids,])),decreasing=TRUE)],]
      duplicatedGenes<-rownames(scaled_depFC)[which(duplicated(rownames(scaled_depFC)))]
    }

  
  
  cls<-intersect(cls,colnames(scaled_depFC))
  
  basal_exp<-basal_exp[,cls]
  scaled_depFC<-scaled_depFC[,cls]
  
  fpkm<-basal_exp
  # fpkm: genes in rows, melanoma cell lines in columns
  fpkm <- as.matrix(fpkm)
  

  # remove cell lines with NAs
  
  fpkm<-fpkm[which(rowSums(is.na(fpkm))==0),]

  # Check input; resolve missing values or duplicate identifiers explicitly
  stopifnot(
    is.numeric(fpkm),
    ncol(fpkm) >= 2,
    !is.null(rownames(fpkm)),
    !is.null(colnames(fpkm)),
    !anyDuplicated(rownames(fpkm)),
    !anyDuplicated(colnames(fpkm)),
    all(is.finite(fpkm)),
    all(fpkm >= 0)
  )
  
  
  scaled_depFC[is.na(scaled_depFC)]<-0
  
  stopifnot(
    is.numeric(scaled_depFC),
    ncol(scaled_depFC) >= 2,
    !is.null(rownames(scaled_depFC)),
    !is.null(colnames(scaled_depFC)),
    !anyDuplicated(rownames(scaled_depFC)),
    !anyDuplicated(colnames(scaled_depFC)),
    all(is.finite(scaled_depFC))
  )
  
  # 1. Keep genes expressed at >= 1 FPKM in >= 10% of lines.
  # This is a pragmatic starting threshold, not a biological boundary.
  keep <- rowMeans(fpkm >= 1) >= 0.10
  expr_log <- log2(fpkm[keep, , drop = FALSE] + 1)
  
  # 2. Remove genes with no variation across cell lines
  gene_sd <- apply(expr_log, 1, sd)
  expr_log <- expr_log[gene_sd > 0, , drop = FALSE]
  
  if(centering){
    gene_mean <- rowMeans(expr_log)
    expr_centered <- sweep(expr_log, 1, gene_mean, FUN = "-")
    expr_log<-expr_centered
  }
  # 3. Centre each gene cell lines
  
  return(list(EXP=expr_log,DEP=scaled_depFC))
  
}

DS_create_analysis_matrix <- function(lfc, padj) {
  
  ## Significant in at least one experiment
  ## NA padj values are treated as non-significant
  sig_gene <- colSums(
    !is.na(padj) & padj < 0.05
  ) >= 1
  
  ## Keep only genes with finite logFC values in every experiment
  complete_gene <- colSums(
    !is.finite(lfc)
  ) == 0
  
  ## Keep only variable genes
  variable_gene <- apply(
    lfc,
    2,
    sd
  ) > 0
  
  ## Combine filters
  genes_use <- sig_gene & complete_gene & variable_gene
  
  X <- lfc[, genes_use, drop = FALSE]
  
  return(X)
}
DS_pca<-function(X){
  pca <- prcomp(X,center = TRUE,scale. = FALSE)
  var_explained <- 100 *pca$sdev^2 / sum(pca$sdev^2)
  plot(var_explained[1:20],type = "b",pch = 16,xlab = "Principal component",ylab = "Variance explained (%)",bty = "l")
  return(pca)
}

DS_plot_pca <- function(
    pca = NULL,
    meta = NULL,
    pc1 = 1,
    pc2 = 2,
    color_by = c("perturbation", "oncogene_cellline", "time"),
    legend_only = FALSE,
    xlim = NULL,
    ylim = NULL
) {
  
  color_by <- match.arg(color_by)
  
  control_targets <- c("RPL15", "AAVS1", "NT")
  
  
  ## ========================================================
  ## Define groups and colours
  ## ========================================================
  
  if (color_by == "perturbation") {
    
    cols <- c(
      "Oncogene"      = "#2166AC",
      "RPL15"         = "#B2182B",
      "AAVS1"         = "#4D9221",
      "Non-targeting" = "grey55"
    )
    
    if (!legend_only) {
      
      group <- ifelse(
        meta$target == "RPL15", "RPL15",
        ifelse(
          meta$target == "AAVS1", "AAVS1",
          ifelse(
            meta$target == "NT", "Non-targeting",
            "Oncogene"
          )
        )
      )
    }
    
    
  } else if (color_by == "oncogene_cellline") {
    
    if (is.null(meta))
      stop("meta must be provided when color_by = 'oncogene_cellline'")
    
    tmp <- meta[
      !(meta$target %in% control_targets),
      c("cell_line", "target"),
      drop = FALSE
    ]
    
    tmp <- unique(tmp)
    
    oncogene_by_cell <- setNames(
      tmp$target,
      tmp$cell_line
    )
    
    group <- paste0(
      meta$cell_line,
      " — ",
      oncogene_by_cell[meta$cell_line]
    )
    
    group_levels <- unique(group)
    
    cols <- setNames(
      hcl.colors(
        length(group_levels),
        palette = "Dark 3"
      ),
      group_levels
    )
    
    
  } else if (color_by == "time") {
    
    if (is.null(meta))
      stop("meta must be provided when color_by = 'time'")
    
    time_days <- meta$time_h / 24
    
    group_levels <- paste0(
      "Day ",
      sort(unique(time_days))
    )
    
    group <- factor(
      paste0("Day ", time_days),
      levels = group_levels
    )
    
    cols <- setNames(
      rev(
        hcl.colors(
          length(group_levels),
          palette = "Viridis"
        )
      ),
      group_levels
    )
  }
  
  
  ## ========================================================
  ## Legend only
  ## ========================================================
  
  if (legend_only) {
    
    plot.new()
    
    legend(
      "center",
      legend = names(cols),
      col = cols,
      pch = 16,
      bty = "n",
      cex = 1.2
    )
    
    return(invisible(NULL))
  }
  
  
  ## ========================================================
  ## PCA + metadata required from here
  ## ========================================================
  
  if (is.null(pca))
    stop("pca must be provided when legend_only = FALSE")
  
  if (is.null(meta))
    stop("meta must be provided when legend_only = FALSE")
  
  
  ## Align metadata
  if (!is.null(rownames(pca$x)) && "guide_id" %in% colnames(meta)) {
    
    meta <- meta[
      match(rownames(pca$x), meta$guide_id),
      ,
      drop = FALSE
    ]
  }
  
  
  ## Variance explained
  var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  
  
  ## Point colours
  point_cols <- cols[as.character(group)]
  
  
  ## If limits were not supplied, determine automatically
  if (is.null(xlim))
    xlim <- range(pca$x[, pc1], na.rm = TRUE)
  
  if (is.null(ylim))
    ylim <- range(pca$x[, pc2], na.rm = TRUE)
  
  
  ## ========================================================
  ## PCA plot
  ## ========================================================
  
  plot(
    pca$x[, pc1],
    pca$x[, pc2],
    pch = 16,
    cex = 1.2,
    col = point_cols,
    xlim = xlim,
    ylim = ylim,
    xlab = paste0(
      "PC", pc1,
      " (", round(var_exp[pc1], 1), "%)"
    ),
    ylab = paste0(
      "PC", pc2,
      " (", round(var_exp[pc2], 1), "%)"
    ),
    bty = "l"
  )
}

DS_plot_pca_by_time <- function(
    pca,
    meta,
    color_by = c("oncogene_cellline", "perturbation"),
    pc1 = 1,
    pc2 = 2,
    times = NULL,
    ncol = NULL,
    add_legend = TRUE,
    padding = 0.05
) {
  
  color_by <- match.arg(color_by)
  
  
  ## ========================================================
  ## Align metadata
  ## ========================================================
  
  meta <- meta[
    match(rownames(pca$x), meta$guide_id),
    ,
    drop = FALSE
  ]
  
  
  ## ========================================================
  ## Time points
  ## ========================================================
  
  if (is.null(times)) {
    times <- sort(unique(meta$time_h))
  }
  
  
  ## ========================================================
  ## Common PCA limits
  ## ========================================================
  
  idx_all <- meta$time_h %in% times
  
  x_range <- range(
    pca$x[idx_all, pc1],
    na.rm = TRUE
  )
  
  y_range <- range(
    pca$x[idx_all, pc2],
    na.rm = TRUE
  )
  
  ## Add a little whitespace around extremes
  x_pad <- diff(x_range) * padding
  y_pad <- diff(y_range) * padding
  
  xlim <- x_range + c(-x_pad, x_pad)
  ylim <- y_range + c(-y_pad, y_pad)
  
  
  ## ========================================================
  ## Panel layout
  ## ========================================================
  
  n_panels <- length(times) + ifelse(add_legend, 1, 0)
  
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_panels))
  }
  
  nrow <- ceiling(n_panels / ncol)
  
  oldpar <- par(
    mfrow = c(nrow, ncol),
    mar = c(4, 4, 3, 1)
  )
  
  on.exit(par(oldpar))
  
  
  ## ========================================================
  ## Plot time points
  ## ========================================================
  
  for (tt in times) {
    
    idx <- which(meta$time_h == tt)
    
    pca_t <- pca
    
    pca_t$x <- pca$x[
      idx,
      ,
      drop = FALSE
    ]
    
    meta_t <- meta[
      idx,
      ,
      drop = FALSE
    ]
    
    DS_plot_pca(
      pca = pca_t,
      meta = meta_t,
      pc1 = pc1,
      pc2 = pc2,
      color_by = color_by,
      legend_only = FALSE,
      xlim = xlim,
      ylim = ylim
    )
    
    title(
      main = paste0(
        "Day ",
        tt / 24
      )
    )
  }
  
  
  ## ========================================================
  ## Shared legend
  ## ========================================================
  
  if (add_legend) {
    
    DS_plot_pca(
      meta = meta,
      color_by = color_by,
      legend_only = TRUE
    )
  }
  
  
  invisible(NULL)
}

DS_plot_pca_by_cellline_perturbation <- function(
    pca,
    meta,
    pc1 = 1,
    pc2 = 2,
    cell_lines = NULL,
    targets = NULL,
    ncol = NULL,
    add_legend = TRUE,
    padding = 0.05
) {
  
  ## ========================================================
  ## Align metadata
  ## ========================================================
  
  meta <- meta[
    match(rownames(pca$x), meta$guide_id),
    ,
    drop = FALSE
  ]
  
  
  ## ========================================================
  ## Optional filtering
  ## ========================================================
  
  if (!is.null(cell_lines)) {
    meta_keep <- meta$cell_line %in% cell_lines
  } else {
    meta_keep <- rep(TRUE, nrow(meta))
  }
  
  if (!is.null(targets)) {
    meta_keep <- meta_keep & meta$target %in% targets
  }
  
  meta_sub <- meta[
    meta_keep,
    ,
    drop = FALSE
  ]
  
  
  ## ========================================================
  ## Define cell line × perturbation groups
  ## ========================================================
  
  groups <- unique(
    meta_sub[, c("cell_line", "target")]
  )
  
  groups <- groups[
    order(groups$cell_line, groups$target),
    ,
    drop = FALSE
  ]
  
  
  ## ========================================================
  ## Common PCA limits
  ## ========================================================
  
  idx_all <- rownames(pca$x) %in% meta_sub$guide_id
  
  x_range <- range(
    pca$x[idx_all, pc1],
    na.rm = TRUE
  )
  
  y_range <- range(
    pca$x[idx_all, pc2],
    na.rm = TRUE
  )
  
  x_pad <- diff(x_range) * padding
  y_pad <- diff(y_range) * padding
  
  xlim <- x_range + c(-x_pad, x_pad)
  ylim <- y_range + c(-y_pad, y_pad)
  
  
  ## ========================================================
  ## Panel layout
  ## ========================================================
  
  n_panels <- nrow(groups) + ifelse(add_legend, 1, 0)
  
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_panels))
  }
  
  nrow <- ceiling(n_panels / ncol)
  
  oldpar <- par(
    mfrow = c(nrow, ncol),
    mar = c(4, 4, 3, 1)
  )
  
  on.exit(par(oldpar))
  
  
  ## ========================================================
  ## One panel per cell line × perturbation
  ## ========================================================
  
  for (i in seq_len(nrow(groups))) {
    
    cl <- groups$cell_line[i]
    tg <- groups$target[i]
    
    idx <- which(
      meta$cell_line == cl &
        meta$target == tg
    )
    
    pca_tmp <- pca
    
    pca_tmp$x <- pca$x[
      idx,
      ,
      drop = FALSE
    ]
    
    meta_tmp <- meta[
      idx,
      ,
      drop = FALSE
    ]
    
    
    DS_plot_pca(
      pca = pca_tmp,
      meta = meta_tmp,
      pc1 = pc1,
      pc2 = pc2,
      color_by = "time",
      legend_only = FALSE,
      xlim = xlim,
      ylim = ylim
    )
    
    
    title(
      main = paste0(
        cl,
        " — ",
        tg
      )
    )
  }
  
  
  ## ========================================================
  ## Shared time legend
  ## ========================================================
  
  if (add_legend) {
    
    DS_plot_pca(
      meta = meta_sub,
      color_by = "time",
      legend_only = TRUE
    )
  }
  
  
  invisible(NULL)
}


DS_compute_DSB_component <- function(lfc, meta) {
  
  cell_lines <- unique(meta$cell_line)
  times <- sort(unique(meta$time_h))
  
  out <- list()
  out_meta <- list()
  
  k <- 1
  
  for (cl in cell_lines) {
    
    for (tt in times) {
      
      m <- meta[
        meta$cell_line == cl &
          meta$time_h == tt,
        ,
        drop = FALSE
      ]
      
      A_ids <- m$guide_id[m$target == "AAVS1"]
      N_ids <- m$guide_id[m$target == "NT"]
      
      if (length(A_ids) == 0 || length(N_ids) == 0)
        next
      
      A <- colMeans(
        lfc[A_ids, , drop = FALSE],
        na.rm = TRUE
      )
      
      N <- colMeans(
        lfc[N_ids, , drop = FALSE],
        na.rm = TRUE
      )
      
      E <- A - N
      
      id <- paste(cl, tt, sep = "_")
      
      out[[id]] <- E
      
      out_meta[[k]] <- data.frame(
        id = id,
        cell_line = cl,
        time_h = tt
      )
      
      k <- k + 1
    }
  }
  
  E_mat <- do.call(rbind, out)
  
  E_meta <- do.call(
    rbind,
    out_meta
  )
  
  rownames(E_meta) <- E_meta$id
  
  return(
    list(
      DSB = E_mat,
      metadata = E_meta
    )
  )
}

DS_cosine <- function(x, y) {
  
  ok <- is.finite(x) & is.finite(y)
  
  sum(x[ok] * y[ok]) /
    (
      sqrt(sum(x[ok]^2)) *
        sqrt(sum(y[ok]^2))
    )
}

DS_compare_DSB_similarity <- function(lfc, meta, DSB_obj) {
  
  times <- sort(unique(meta$time_h))
  
  res <- data.frame()
  
  for (tt in times) {
    
    cls <- unique(
      meta$cell_line[
        meta$time_h == tt
      ]
    )
    
    if (length(cls) < 2)
      next
    
    pairs <- combn(cls, 2)
    
    for (j in seq_len(ncol(pairs))) {
      
      cl1 <- pairs[1, j]
      cl2 <- pairs[2, j]
      
      m1 <- meta[
        meta$cell_line == cl1 &
          meta$time_h == tt,
      ]
      
      m2 <- meta[
        meta$cell_line == cl2 &
          meta$time_h == tt,
      ]
      
      A1 <- colMeans(
        lfc[
          m1$guide_id[m1$target == "AAVS1"],
          ,
          drop = FALSE
        ],
        na.rm = TRUE
      )
      
      A2 <- colMeans(
        lfc[
          m2$guide_id[m2$target == "AAVS1"],
          ,
          drop = FALSE
        ],
        na.rm = TRUE
      )
      
      E1 <- DSB_obj$DSB[
        paste(cl1, tt, sep = "_"),
      ]
      
      E2 <- DSB_obj$DSB[
        paste(cl2, tt, sep = "_"),
      ]
      
      res <- rbind(
        res,
        data.frame(
          time_h = tt,
          cell_line_1 = cl1,
          cell_line_2 = cl2,
          raw_AAVS1_cosine = DS_cosine(A1, A2),
          DSB_cosine = DS_cosine(E1, E2)
        )
      )
    }
  }
  
  return(res)
}

DS_cosine <- function(x, y) {
  
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) == 0)
    return(NA_real_)
  
  nx <- sqrt(sum(x[ok]^2))
  ny <- sqrt(sum(y[ok]^2))
  
  if (nx == 0 || ny == 0)
    return(NA_real_)
  
  sum(x[ok] * y[ok]) / (nx * ny)
}

DS_compute_common_trace <- function(
    lfc,
    meta,
    exclude_targets = "NT",
    consensus = c("median", "mean")
) {
  
  consensus <- match.arg(consensus)
  
  cell_lines <- unique(meta$cell_line)
  times <- sort(unique(meta$time_h))
  
  traces <- list()
  trace_meta <- list()
  
  k <- 1
  
  for (cl in cell_lines) {
    
    for (tt in times) {
      
      m <- meta[
        meta$cell_line == cl &
          meta$time_h == tt,
        ,
        drop = FALSE
      ]
      
      ## Exclude NT from estimation
      m <- m[
        !(m$target %in% exclude_targets),
        ,
        drop = FALSE
      ]
      
      if (nrow(m) == 0)
        next
      
      
      ## -----------------------------------------------
      ## First average guides within each target
      ## -----------------------------------------------
      
      targets <- unique(m$target)
      
      target_profiles <- matrix(
        NA_real_,
        nrow = length(targets),
        ncol = ncol(lfc)
      )
      
      rownames(target_profiles) <- targets
      colnames(target_profiles) <- colnames(lfc)
      
      for (i in seq_along(targets)) {
        
        tg <- targets[i]
        
        ids <- m$guide_id[
          m$target == tg
        ]
        
        target_profiles[i, ] <- colMeans(
          lfc[ids, , drop = FALSE],
          na.rm = TRUE
        )
      }
      
      
      ## -----------------------------------------------
      ## Consensus across perturbations
      ## -----------------------------------------------
      
      if (consensus == "median") {
        
        common <- apply(
          target_profiles,
          2,
          median,
          na.rm = TRUE
        )
        
      } else {
        
        common <- colMeans(
          target_profiles,
          na.rm = TRUE
        )
      }
      
      
      id <- paste(cl, tt, sep = "_")
      
      traces[[id]] <- common
      
      trace_meta[[k]] <- data.frame(
        id = id,
        cell_line = cl,
        time_h = tt,
        n_targets = length(targets)
      )
      
      k <- k + 1
    }
  }
  
  
  trace_mat <- do.call(
    rbind,
    traces
  )
  
  trace_meta <- do.call(
    rbind,
    trace_meta
  )
  
  rownames(trace_meta) <- trace_meta$id
  
  
  return(
    list(
      trace = trace_mat,
      metadata = trace_meta
    )
  )
}

DS_validate_common_trace_with_NT <- function(
    common_obj,
    lfc,
    meta
) {
  
  res <- data.frame()
  
  for (i in seq_len(nrow(common_obj$metadata))) {
    
    cl <- common_obj$metadata$cell_line[i]
    tt <- common_obj$metadata$time_h[i]
    id <- common_obj$metadata$id[i]
    
    m <- meta[
      meta$cell_line == cl &
        meta$time_h == tt &
        meta$target == "NT",
      ,
      drop = FALSE
    ]
    
    if (nrow(m) == 0)
      next
    
    NT <- colMeans(
      lfc[m$guide_id, , drop = FALSE],
      na.rm = TRUE
    )
    
    common <- common_obj$trace[id, ]
    
    res <- rbind(
      res,
      data.frame(
        cell_line = cl,
        time_h = tt,
        cosine_common_NT = DS_cosine(
          common,
          NT
        )
      )
    )
  }
  
  return(res)
}


DS_compare_common_trace <- function(common_obj) {
  
  times <- sort(
    unique(common_obj$metadata$time_h)
  )
  
  res <- data.frame()
  
  for (tt in times) {
    
    d <- common_obj$metadata[
      common_obj$metadata$time_h == tt,
      ,
      drop = FALSE
    ]
    
    if (nrow(d) < 2)
      next
    
    pairs <- combn(
      d$cell_line,
      2
    )
    
    for (j in seq_len(ncol(pairs))) {
      
      cl1 <- pairs[1, j]
      cl2 <- pairs[2, j]
      
      id1 <- d$id[
        d$cell_line == cl1
      ]
      
      id2 <- d$id[
        d$cell_line == cl2
      ]
      
      x <- common_obj$trace[id1, ]
      y <- common_obj$trace[id2, ]
      
      res <- rbind(
        res,
        data.frame(
          time_h = tt,
          cell_line_1 = cl1,
          cell_line_2 = cl2,
          cosine = DS_cosine(x, y)
        )
      )
    }
  }
  
  return(res)
}



