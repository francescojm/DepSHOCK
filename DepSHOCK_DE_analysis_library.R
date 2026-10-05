DS_consensus_clusters <- function(
    clusterings,
    k = NULL,
    plot = TRUE
) {
  
  ## --------------------------------------------------------
  ## Basic checks
  ## --------------------------------------------------------
  
  if (!is.list(clusterings) || length(clusterings) < 2) {
    stop("clusterings must be a list containing at least two clusterings.")
  }
  
  for (i in seq_along(clusterings)) {
    
    if (is.null(names(clusterings[[i]]))) {
      stop("Each clustering vector must have gene names.")
    }
    
    if (anyDuplicated(names(clusterings[[i]]))) {
      stop("Duplicated gene names found in clustering ", i)
    }
  }
  
  
  ## --------------------------------------------------------
  ## Use genes represented in ALL clusterings
  ## --------------------------------------------------------
  
  common_genes <- Reduce(
    intersect,
    lapply(
      clusterings,
      names
    )
  )
  
  if (length(common_genes) < 2) {
    stop("Fewer than two genes are shared across all clusterings.")
  }
  
  cat(
    "Clusterings:",
    length(clusterings),
    "\n"
  )
  
  cat(
    "Genes common to all clusterings:",
    length(common_genes),
    "\n"
  )
  
  
  ## --------------------------------------------------------
  ## Consensus co-clustering matrix
  ## --------------------------------------------------------
  
  consensus <- matrix(
    0,
    nrow = length(common_genes),
    ncol = length(common_genes),
    dimnames = list(
      common_genes,
      common_genes
    )
  )
  
  
  for (cl in clusterings) {
    
    labels <- cl[
      common_genes
    ]
    
    same_cluster <- outer(
      labels,
      labels,
      FUN = "=="
    )
    
    consensus <- consensus +
      same_cluster
  }
  
  
  consensus <- consensus /
    length(clusterings)
  
  
  diag(consensus) <- 1
  
  
  ## --------------------------------------------------------
  ## Consensus distance
  ## --------------------------------------------------------
  
  consensus_distance <- as.dist(
    1 - consensus
  )
  
  
  ## --------------------------------------------------------
  ## Hierarchical clustering of consensus matrix
  ## --------------------------------------------------------
  
  hc <- hclust(
    consensus_distance,
    method = "average"
  )
  
  
  ## --------------------------------------------------------
  ## Optional consensus cluster assignment
  ## --------------------------------------------------------
  
  consensus_clusters <- NULL
  cluster_strength <- NULL
  gene_stability <- NULL
  
  
  if (!is.null(k)) {
    
    consensus_clusters <- cutree(
      hc,
      k = k
    )
    
    
    ## ------------------------------------------------------
    ## Mean within-cluster consensus
    ## ------------------------------------------------------
    
    cluster_strength <- sapply(
      seq_len(k),
      function(cl) {
        
        genes <- names(
          consensus_clusters[
            consensus_clusters == cl
          ]
        )
        
        if (length(genes) == 1) {
          return(NA_real_)
        }
        
        M <- consensus[
          genes,
          genes,
          drop = FALSE
        ]
        
        mean(
          M[
            upper.tri(M)
          ],
          na.rm = TRUE
        )
      }
    )
    
    names(cluster_strength) <- paste0(
      "C",
      seq_len(k)
    )
    
    
    ## ------------------------------------------------------
    ## Per-gene stability
    ##
    ## How consistently does each gene co-cluster with
    ## the other members of its consensus cluster?
    ## ------------------------------------------------------
    
    gene_stability <- sapply(
      common_genes,
      function(g) {
        
        cl <- consensus_clusters[g]
        
        mates <- names(
          consensus_clusters[
            consensus_clusters == cl
          ]
        )
        
        mates <- setdiff(
          mates,
          g
        )
        
        if (length(mates) == 0) {
          return(NA_real_)
        }
        
        mean(
          consensus[
            g,
            mates
          ],
          na.rm = TRUE
        )
      }
    )
  }
  
  
  ## --------------------------------------------------------
  ## Plot consensus matrix
  ## --------------------------------------------------------
  
  if (plot) {
    
    dend <- as.dendrogram(
      hc
    )
    
    heatmap(
      consensus,
      Rowv = dend,
      Colv = dend,
      scale = "none",
      symm = TRUE,
      labRow = NA,
      labCol = NA,
      margins = c(2, 2),
      xlab = "",
      ylab = ""
    )
  }
  
  
  ## --------------------------------------------------------
  ## Return
  ## --------------------------------------------------------
  
  return(
    list(
      consensus = consensus,
      distance = consensus_distance,
      hclust = hc,
      clusters = consensus_clusters,
      cluster_strength = cluster_strength,
      gene_stability = gene_stability,
      common_genes = common_genes
    )
  )
}

DS_cluster_temporal_profiles <- function(
    subFC,
    k = 6,
    distance = c("euclidean", "correlation"),
    linkage = "average",
    time_h = c(24, 48, 72, 96, 120),
    plot = TRUE
) {
  
  distance <- match.arg(distance)
  
  X <- as.matrix(subFC)
  storage.mode(X) <- "numeric"
  
  ## --------------------------------------------------------
  ## Remove problematic genes
  ## --------------------------------------------------------
  
  keep <- apply(
    X,
    1,
    function(x) {
      all(is.finite(x)) &&
        sd(x) > 0
    }
  )
  
  X <- X[
    keep,
    ,
    drop = FALSE
  ]
  
  cat(
    "Genes retained:",
    nrow(X),
    "\n"
  )
  
  
  ## --------------------------------------------------------
  ## Distance between temporal profiles
  ## NO scaling: actual logFC values are retained
  ## --------------------------------------------------------
  
  if (distance == "correlation") {
    
    D <- as.dist(
      1 - cor(
        t(X),
        method = "pearson"
      )
    )
    
  } else {
    
    D <- dist(
      X,
      method = "euclidean"
    )
  }
  
  
  ## --------------------------------------------------------
  ## Hierarchical clustering
  ## --------------------------------------------------------
  
  hc <- hclust(
    D,
    method = linkage
  )
  
  clusters <- cutree(
    hc,
    k = k
  )
  
  
  ## --------------------------------------------------------
  ## Cluster colours
  ## --------------------------------------------------------
  
  cluster_cols <- hcl.colors(
    k,
    palette = "Dark 3"
  )
  
  names(cluster_cols) <- paste0(
    "C",
    seq_len(k)
  )
  
  
  ## --------------------------------------------------------
  ## Cluster sizes
  ## --------------------------------------------------------
  
  cluster_sizes <- as.numeric(
    table(
      factor(
        clusters,
        levels = seq_len(k)
      )
    )
  )
  
  names(cluster_sizes) <- paste0(
    "C",
    seq_len(k)
  )
  
  
  ## --------------------------------------------------------
  ## Mean RAW logFC profile for each cluster
  ## --------------------------------------------------------
  
  cluster_profiles <- matrix(
    NA_real_,
    nrow = k,
    ncol = ncol(X),
    dimnames = list(
      paste0(
        "Cluster_",
        seq_len(k)
      ),
      colnames(X)
    )
  )
  
  for (cl in seq_len(k)) {
    
    genes_cl <- names(
      clusters[
        clusters == cl
      ]
    )
    
    cluster_profiles[
      cl,
    ] <- colMeans(
      X[
        genes_cl,
        ,
        drop = FALSE
      ]
    )
  }
  
  
  ## --------------------------------------------------------
  ## Plot
  ## --------------------------------------------------------
  
  if (plot) {
    
    oldpar <- par(
      mfrow = c(1, 2)
    )
    
    on.exit(
      par(oldpar)
    )
    
    
    ## ======================================================
    ## Dendrogram
    ## ======================================================
    
    plot(
      hc,
      labels = FALSE,
      hang = -1,
      main = paste(
        "Clustering -",
        distance,
        "distance"
      ),
      xlab = "",
      sub = ""
    )
    
    rect.hclust(
      hc,
      k = k
    )
    
    
    ## ======================================================
    ## Temporal profiles
    ## ======================================================
    
    ## Add artificial t = 0 ONLY for visualization
    
    profiles_plot <- cbind(
      "0h" = 0,
      cluster_profiles
    )
    
    time_plot <- c(
      0,
      time_h
    )
    
    
    matplot(
      time_plot[1:nrow(t(profiles_plot))],
      t(profiles_plot),
      type = "o",
      pch = 16,
      lty = 1,
      lwd = 2,
      col = cluster_cols,
      xlab = "Time after induction (h)",
      ylab = "Mean logFC",
      main = "Cluster temporal profiles"
    )
    
    
    abline(
      h = 0,
      lty = 2,
      col = "grey60"
    )
    
    
    legend(
      "topright",
      legend = paste0(
        "C",
        seq_len(k),
        " (n=",
        cluster_sizes,
        ")"
      ),
      col = cluster_cols,
      lty = 1,
      pch = 16,
      lwd = 2,
      bty = "n"
    )
  }
  
  
  ## --------------------------------------------------------
  ## Return
  ## --------------------------------------------------------
  
  return(
    list(
      clusters = clusters,
      hclust = hc,
      fc = X,
      cluster_profiles = cluster_profiles,
      cluster_sizes = cluster_sizes,
      cluster_cols = cluster_cols
    )
  )
}


DS_split_DE_by_perturbation <- function(DE_obj) {
  
  meta <- DE_obj$contrast_metadata
  fc   <- DE_obj$fc
  adjp <- DE_obj$adjp
  
  stopifnot(
    identical(colnames(fc), colnames(adjp)),
    all(meta$experiment_id %in% colnames(fc))
  )
  
  guides <- unique(meta$guide)
  
  out <- setNames(
    vector("list", length(guides)),
    guides
  )
  
  for (g in guides) {
    
    m <- meta[
      meta$guide == g,
      ,
      drop = FALSE
    ]
    
    ## sort by increasing time
    m <- m[
      order(m$time_h),
      ,
      drop = FALSE
    ]
    
    ids <- m$experiment_id
    
    out[[g]] <- list(
      
      lfc = fc[
        ,
        ids,
        drop = FALSE
      ],
      
      adjp = adjp[
        ,
        ids,
        drop = FALSE
      ]
      
    )
    
  }
  
  return(out)
}


DS_BRAF_weighted_score <- function(
    DE_obj,
    gex,
    time_h,
    padj_thr = 0.05,
    min_logFC = 0
) {
  
  meta <- DE_obj$contrast_metadata
  fc   <- DE_obj$fc
  adjp <- DE_obj$adjp
  
  ## BRAF contrasts at selected time
  m <- meta[
    meta$target == "BRAF" &
      meta$time_h == time_h,
    ,
    drop = FALSE
  ]
  
  ids <- m$experiment_id
  
  fc_sub <- fc[, ids, drop = FALSE]
  p_sub  <- adjp[, ids, drop = FALSE]
  
  ## Require concordant downregulation across guides
  keep <- (
    rowSums(fc_sub < -min_logFC) == ncol(fc_sub) &
      rowSums(p_sub < padj_thr) == ncol(p_sub)
  )
  
  fc_sub <- fc_sub[keep, , drop = FALSE]
  
  ## Consensus KO effect
  mean_fc <- rowMeans(fc_sub)
  
  genes <- intersect(
    names(mean_fc),
    colnames(gex)
  )
  
  mean_fc <- mean_fc[genes]
  
  ## Standardise each gene across melanoma cell lines
  X <- gex[, genes, drop = FALSE]
  
  X <- scale(
    X,
    center = TRUE,
    scale = TRUE
  )
  
  ## stronger KO-down genes get larger positive weights
  w <- -mean_fc
  
  ## prevent total magnitude depending on signature size
  w <- w / sum(w)
  
  score <- as.numeric(
    X %*% w
  )
  
  names(score) <- rownames(X)
  
  return(
    list(
      score = score,
      genes = genes,
      weights = w
    )
  )
}

DS_get_KO_down_signature <- function(
    DE_obj,
    target = "BRAF",
    time_h,
    padj_thr = 0.05,
    min_logFC = 0,
    require_all_guides = TRUE
) {
  
  meta <- DE_obj$contrast_metadata
  fc   <- DE_obj$fc
  adjp <- DE_obj$adjp
  
  ## contrasts for target/time
  m <- meta[
    meta$target == target &
      meta$time_h == time_h,
    ,
    drop = FALSE
  ]
  
  if (nrow(m) == 0) {
    stop("No contrasts found for target/time.")
  }
  
  ids <- m$experiment_id
  
  fc_sub <- fc[, ids, drop = FALSE]
  p_sub  <- adjp[, ids, drop = FALSE]
  
  ## downregulated and significant
  down <- (
    fc_sub < -min_logFC &
      p_sub < padj_thr
  )
  
  if (require_all_guides) {
    
    keep <- rowSums(down) == ncol(down)
    
  } else {
    
    ## at least one guide significant,
    ## but direction must agree across all guides
    keep <- rowSums(down) >= 1 &
      rowSums(fc_sub < 0) == ncol(fc_sub)
  }
  
  genes <- rownames(fc_sub)[keep]
  
  cat("Target:", target, "\n")
  cat("Time:", time_h, "h\n")
  cat("Guides:", paste(m$guide, collapse = ", "), "\n")
  cat("Downregulated signature genes:", length(genes), "\n")
  
  return(genes)
}

# ============================================================
# Helper: find a column robustly
# ============================================================

DS_find_column <- function(x, candidates) {
  
  clean <- function(z) {
    tolower(gsub("[^a-z0-9]", "", z))
  }
  
  nms <- clean(colnames(x))
  cand <- clean(candidates)
  
  hit <- match(cand, nms)
  hit <- hit[!is.na(hit)]
  
  if (length(hit) == 0) {
    stop(
      "Could not identify required column. Tried: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  colnames(x)[hit[1]]
}


# ============================================================
# 1. Build library -> biological sample mapping
# ============================================================

DS_build_sample_map <- function(
    sgRNA_id_path,
    id_map_path,
    cell_line
) {
  
  sg <- fread(
    sgRNA_id_path,
    data.table = FALSE,
    check.names = FALSE
  )
  
  id <- fread(
    id_map_path,
    data.table = FALSE,
    check.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Identify relevant columns
  # ----------------------------------------------------------
  
  id_library_col <- DS_find_column(
    id,
    c(
      "Sample_ID",
      "Sample ID"
    )
  )
  
  id_user_col <- DS_find_column(
    id,
    c(
      "User_Sample_ID",
      "User Sample ID"
    )
  )
  
  sg_user_col <- DS_find_column(
    sg,
    c(
      "Sample name",
      "Sample_name",
      "Sample"
    )
  )
  
  sg_full_col <- DS_find_column(
    sg,
    c(
      "Full sample name",
      "Full_sample_name",
      "Full Sample Name"
    )
  )
  
  
  # ----------------------------------------------------------
  # Run-end report can contain repeated records
  # ----------------------------------------------------------
  
  run_map <- unique(
    data.frame(
      library_id = as.character(id[[id_library_col]]),
      user_id    = as.character(id[[id_user_col]]),
      stringsAsFactors = FALSE
    )
  )
  
  
  bio_map <- unique(
    data.frame(
      user_id = as.character(sg[[sg_user_col]]),
      sample_name = trimws(as.character(sg[[sg_full_col]])),
      stringsAsFactors = FALSE
    )
  )
  
  
  # ----------------------------------------------------------
  # Check one-to-one mapping
  # ----------------------------------------------------------
  
  if (anyDuplicated(run_map$library_id)) {
    stop("A sequencing library maps to multiple user IDs.")
  }
  
  if (anyDuplicated(bio_map$user_id)) {
    stop("A user ID maps to multiple biological sample names.")
  }
  
  
  map <- merge(
    run_map,
    bio_map,
    by = "user_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  
  if (anyNA(map$sample_name)) {
    stop(
      "Some sequencing libraries could not be mapped ",
      "to a biological sample."
    )
  }
  
  
  # ----------------------------------------------------------
  # Parse biological sample names
  #
  # Expected examples:
  #
  # sgAAVS1 DMSO R1
  # sgAAVS1 24h R2
  # sgBRAF #1 day3 R3
  # M14 parental R1
  # ----------------------------------------------------------
  
  map$cell_line <- cell_line
  
  
  # replicate
  map$replicate <- as.integer(
    sub(
      ".*\\bR([0-9]+)\\b.*",
      "\\1",
      map$sample_name,
      ignore.case = TRUE
    )
  )
  
  
  # condition
  map$condition <- NA_character_
  
  map$condition[
    grepl("\\bDMSO\\b", map$sample_name, ignore.case = TRUE)
  ] <- "DMSO"
  
  map$condition[
    grepl("\\b24h\\b", map$sample_name, ignore.case = TRUE)
  ] <- "24h"
  
  for (d in 2:5) {
    
    hit <- grepl(
      paste0("\\bday\\s*", d, "\\b"),
      map$sample_name,
      ignore.case = TRUE
    )
    
    map$condition[hit] <- paste0("day", d)
  }
  
  map$condition[
    grepl("\\bparental\\b", map$sample_name, ignore.case = TRUE)
  ] <- "parental"
  
  
  # time in hours
  map$time_h <- NA_real_
  
  map$time_h[map$condition == "24h"]  <- 24
  map$time_h[map$condition == "day2"] <- 48
  map$time_h[map$condition == "day3"] <- 72
  map$time_h[map$condition == "day4"] <- 96
  map$time_h[map$condition == "day5"] <- 120
  
  
  # ----------------------------------------------------------
  # Parse sgRNA
  # ----------------------------------------------------------
  
  is_sg <- grepl(
    "\\bsg",
    map$sample_name,
    ignore.case = TRUE
  )
  
  map$guide <- NA_character_
  
  
  tmp <- map$sample_name[is_sg]
  
  # remove replicate
  tmp <- gsub(
    "\\s+R[0-9]+\\b.*$",
    "",
    tmp,
    ignore.case = TRUE
  )
  
  # remove condition
  tmp <- gsub(
    "\\s+(DMSO|24h|day\\s*[2-5])\\b.*$",
    "",
    tmp,
    ignore.case = TRUE
  )
  
  tmp <- trimws(tmp)
  
  # normalize '# 1' / ' #1' etc.
  tmp <- gsub("\\s*#\\s*", "#", tmp)
  tmp <- gsub("\\s+", "", tmp)
  
  map$guide[is_sg] <- tmp
  
  
  # target from guide
  map$target <- NA_character_
  
  map$target[is_sg] <- sub(
    "^sg",
    "",
    map$guide[is_sg],
    ignore.case = TRUE
  )
  
  map$target[is_sg] <- sub(
    "#.*$",
    "",
    map$target[is_sg]
  )
  
  
  # ----------------------------------------------------------
  # QC
  # ----------------------------------------------------------
  
  de_samples <- map$condition != "parental"
  
  if (any(is.na(map$condition))) {
    
    warning(
      "Unparsed sample names:\n",
      paste(
        map$sample_name[is.na(map$condition)],
        collapse = "\n"
      )
    )
  }
  
  if (any(is.na(map$guide[de_samples]))) {
    
    warning(
      "Some non-parental samples have no parsed guide."
    )
  }
  
  
  map <- map[
    order(map$library_id),
    ,
    drop = FALSE
  ]
  
  return(map)
}


# ============================================================
# 2. Read Salmon gene count matrix
# ============================================================

DS_read_salmon_counts <- function(
    rawCount_path,
    sample_map
) {
  
  x <- fread(
    rawCount_path,
    data.table = FALSE,
    check.names = FALSE
  )
  
  
  if (!all(c("gene_id", "gene_name") %in% colnames(x))) {
    stop(
      "Count file must contain gene_id and gene_name columns."
    )
  }
  
  
  gene_annotation <- x[
    ,
    c("gene_id", "gene_name"),
    drop = FALSE
  ]
  
  
  sample_cols <- sample_map$library_id
  
  missing_samples <- setdiff(
    sample_cols,
    colnames(x)
  )
  
  if (length(missing_samples) > 0) {
    
    stop(
      "These libraries are missing from the count matrix:\n",
      paste(missing_samples, collapse = ", ")
    )
  }
  
  
  counts <- as.matrix(
    x[
      ,
      sample_cols,
      drop = FALSE
    ]
  )
  
  storage.mode(counts) <- "numeric"
  
  rownames(counts) <- gene_annotation$gene_id
  
  
  if (anyNA(counts)) {
    stop("Count matrix contains NA values.")
  }
  
  if (any(!is.finite(counts))) {
    stop("Count matrix contains non-finite values.")
  }
  
  if (any(counts < 0)) {
    stop("Count matrix contains negative values.")
  }
  
  
  return(
    list(
      counts = counts,
      genes = gene_annotation
    )
  )
}


# ============================================================
# 3. Differential expression
#
# Each guide/time is compared with the SAME guide under DMSO
# ============================================================


library(readxl)

DS_run_experiment_DE <- function(
    experiment,
    path_file,
    path_location,
    output_path,
    sheet = 1,
    min_replicates = 2,
    overwrite = FALSE
) {
  
  ## ========================================================
  ## Load experiment configuration
  ## ========================================================
  
  config <- readxl::read_excel(
    path_file,
    sheet = sheet
  )
  
  config <- as.data.frame(
    config,
    stringsAsFactors = FALSE
  )
  
  
  ## ========================================================
  ## Check required columns
  ## ========================================================
  
  required_cols <- c(
    "experiment",
    "cell_line",
    "sgRNA_id_path",
    "id_map_path",
    "rawCount_path"
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(config)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      "Missing required columns in path file: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  
  ## ========================================================
  ## Find requested experiment
  ## ========================================================
  
  hit <- which(
    config$experiment == experiment
  )
  
  if (length(hit) == 0) {
    
    stop(
      "Experiment '",
      experiment,
      "' not found in configuration file."
    )
  }
  
  if (length(hit) > 1) {
    
    stop(
      "Experiment '",
      experiment,
      "' occurs more than once in configuration file."
    )
  }
  
  cfg <- config[
    hit,
    ,
    drop = FALSE
  ]
  
  
  ## ========================================================
  ## Helper to prepend common path
  ## ========================================================
  
  prefix_path <- function(prefix, x) {
    
    x <- as.character(x)
    
    if (is.na(x) || x == "") {
      stop("Empty path found in configuration file.")
    }
    
    ## Remove leading / from config path so file.path()
    ## treats it as relative to path_location
    x <- sub("^/+", "", x)
    
    file.path(
      path.expand(prefix),
      x
    )
  }
  
  
  ## ========================================================
  ## Extract and construct paths
  ## ========================================================
  
  cell_line <- as.character(
    cfg$cell_line
  )
  
  sgRNA_id_path <- prefix_path(
    path_location,
    cfg$sgRNA_id_path
  )
  
  id_map_path <- prefix_path(
    path_location,
    cfg$id_map_path
  )
  
  rawCount_path <- prefix_path(
    path_location,
    cfg$rawCount_path
  )
  
  
  ## ========================================================
  ## Check input files
  ## ========================================================
  
  input_files <- c(
    sgRNA_id_path = sgRNA_id_path,
    id_map_path = id_map_path,
    rawCount_path = rawCount_path
  )
  
  missing_files <- input_files[
    !file.exists(input_files)
  ]
  
  if (length(missing_files) > 0) {
    
    stop(
      "The following input files do not exist:\n",
      paste(
        names(missing_files),
        "=",
        missing_files,
        collapse = "\n"
      )
    )
  }
  
  
  ## ========================================================
  ## Prepare output directory
  ## ========================================================
  
  if (!dir.exists(output_path)) {
    
    dir.create(
      output_path,
      recursive = TRUE
    )
  }
  
  
  output_file <- file.path(
    output_path,
    paste0(
      experiment,
      "_DE.Robj"
    )
  )
  
  
  if (
    file.exists(output_file) &&
    !overwrite
  ) {
    
    stop(
      "Output file already exists:\n",
      output_file,
      "\nUse overwrite = TRUE to replace it."
    )
  }
  
  
  ## ========================================================
  ## Run differential expression
  ## ========================================================
  
  cat(
    "\n========================================\n",
    "Running experiment: ", experiment, "\n",
    "Cell line: ", cell_line, "\n",
    "========================================\n\n",
    sep = ""
  )
  
  cat(
    "sgRNA annotation:\n",
    sgRNA_id_path,
    "\n\n"
  )
  
  cat(
    "Library ID mapping:\n",
    id_map_path,
    "\n\n"
  )
  
  cat(
    "Count matrix:\n",
    rawCount_path,
    "\n\n"
  )
  
  
  DE_obj <- DS_build_DE_matrices(
    sgRNA_id_path = sgRNA_id_path,
    id_map_path = id_map_path,
    rawCount_path = rawCount_path,
    cell_line = cell_line,
    min_replicates = min_replicates
  )
  
  
  ## ========================================================
  ## Store provenance
  ## ========================================================
  
  DE_obj$experiment <- experiment
  
  DE_obj$path_location <- path_location
  
  DE_obj$input_files <- list(
    sgRNA_id_path = sgRNA_id_path,
    id_map_path = id_map_path,
    rawCount_path = rawCount_path
  )
  
  
  ## ========================================================
  ## Save object
  ## ========================================================
  
  save(
    DE_obj,
    file = output_file
  )
  
  
  cat(
    "\nDE analysis completed.\n",
    "Saved to:\n",
    output_file,
    "\n"
  )
  
  
  invisible(DE_obj)
}

DS_run_differential_expression <- function(
    counts,
    gene_annotation,
    sample_map,
    min_replicates = 2
) {
  
  # ----------------------------------------------------------
  # Remove parental samples from differential expression
  # ----------------------------------------------------------
  
  keep_sample <- (
    !is.na(sample_map$guide) &
      !is.na(sample_map$condition) &
      sample_map$condition != "parental"
  )
  
  meta <- sample_map[
    keep_sample,
    ,
    drop = FALSE
  ]
  
  counts <- counts[
    ,
    meta$library_id,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Define guide-condition groups
  # ----------------------------------------------------------
  
  meta$group <- paste(
    meta$guide,
    meta$condition,
    sep = "__"
  )
  
  meta$group <- factor(meta$group)
  
  
  design <- model.matrix(
    ~ 0 + group,
    data = meta
  )
  
  group_levels <- levels(meta$group)
  
  safe_group_names <- make.names(
    group_levels,
    unique = TRUE
  )
  
  colnames(design) <- safe_group_names
  
  group_to_design <- setNames(
    safe_group_names,
    group_levels
  )
  
  
  # ----------------------------------------------------------
  # edgeR object + expression filtering
  # ----------------------------------------------------------
  
  y <- DGEList(
    counts = counts
  )
  
  
  keep_gene <- filterByExpr(
    y,
    design = design
  )
  
  cat(
    "Genes before expression filtering:",
    nrow(y),
    "\n"
  )
  
  cat(
    "Genes retained:",
    sum(keep_gene),
    "\n"
  )
  
  
  y <- y[
    keep_gene,
    ,
    keep.lib.sizes = FALSE
  ]
  
  
  # TMM normalization
  y <- calcNormFactors(y)
  
  
  # ----------------------------------------------------------
  # Voom
  # ----------------------------------------------------------
  
  v <- voom(
    y,
    design,
    plot = FALSE
  )
  
  fit <- lmFit(
    v,
    design
  )
  
  
  # ----------------------------------------------------------
  # Construct all guide/time vs same-guide DMSO contrasts
  # ----------------------------------------------------------
  
  guides <- unique(meta$guide)
  
  contrast_list <- list()
  contrast_meta <- list()
  
  k <- 1
  
  
  for (g in guides) {
    
    gmeta <- meta[
      meta$guide == g,
      ,
      drop = FALSE
    ]
    
    
    baseline_group <- paste(
      g,
      "DMSO",
      sep = "__"
    )
    
    
    if (!(baseline_group %in% group_levels)) {
      
      warning(
        "No DMSO baseline for guide ",
        g,
        "; skipping."
      )
      
      next
    }
    
    
    # require sufficient DMSO replicates
    if (
      sum(gmeta$condition == "DMSO") < min_replicates
    ) {
      
      warning(
        "Insufficient DMSO replicates for ",
        g,
        "; skipping."
      )
      
      next
    }
    
    
    conditions <- unique(
      gmeta$condition[
        gmeta$condition != "DMSO"
      ]
    )
    
    
    for (cond in conditions) {
      
      n_rep <- sum(
        gmeta$condition == cond
      )
      
      if (n_rep < min_replicates) {
        
        warning(
          "Skipping ",
          g,
          " ",
          cond,
          ": only ",
          n_rep,
          " replicate(s)."
        )
        
        next
      }
      
      
      treatment_group <- paste(
        g,
        cond,
        sep = "__"
      )
      
      
      contrast <- rep(
        0,
        ncol(design)
      )
      
      names(contrast) <- colnames(design)
      
      
      contrast[
        group_to_design[treatment_group]
      ] <- 1
      
      contrast[
        group_to_design[baseline_group]
      ] <- -1
      
      
      contrast_name <- paste(
        g,
        cond,
        sep = "__"
      )
      
      
      contrast_list[[contrast_name]] <- contrast
      
      
      time_h <- unique(
        gmeta$time_h[
          gmeta$condition == cond
        ]
      )
      
      target <- unique(
        gmeta$target
      )
      
      
      contrast_meta[[k]] <- data.frame(
        contrast = contrast_name,
        guide = g,
        target = target[1],
        condition = cond,
        time_h = time_h[1],
        stringsAsFactors = FALSE
      )
      
      k <- k + 1
    }
  }
  
  
  if (length(contrast_list) == 0) {
    stop("No valid differential-expression contrasts generated.")
  }
  
  
  contrast_matrix <- do.call(
    cbind,
    contrast_list
  )
  
  
  # ----------------------------------------------------------
  # Differential expression
  # ----------------------------------------------------------
  
  fit2 <- contrasts.fit(
    fit,
    contrast_matrix
  )
  
  fit2 <- eBayes(
    fit2,
    robust = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Build gene_id x contrast matrices
  # ----------------------------------------------------------
  
  gene_ids <- rownames(v$E)
  
  fc <- matrix(
    NA_real_,
    nrow = length(gene_ids),
    ncol = ncol(contrast_matrix),
    dimnames = list(
      gene_ids,
      colnames(contrast_matrix)
    )
  )
  
  adjp <- fc
  
  
  for (i in seq_len(ncol(contrast_matrix))) {
    
    tt <- topTable(
      fit2,
      coef = i,
      number = Inf,
      sort.by = "none",
      adjust.method = "BH"
    )
    
    tt <- tt[
      gene_ids,
      ,
      drop = FALSE
    ]
    
    fc[, i] <- tt$logFC
    adjp[, i] <- tt$adj.P.Val
  }
  
  
  # ----------------------------------------------------------
  # Gene annotations for retained genes
  # ----------------------------------------------------------
  
  ann <- gene_annotation[
    match(
      gene_ids,
      gene_annotation$gene_id
    ),
    ,
    drop = FALSE
  ]
  
  
  # Mean expression used to resolve duplicated symbols
  mean_expression <- rowMeans(
    v$E,
    na.rm = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Convert from gene_id to gene symbol
  # ----------------------------------------------------------
  
  symbol <- trimws(
    as.character(ann$gene_name)
  )
  
  valid_symbol <- (
    !is.na(symbol) &
      symbol != ""
  )
  
  
  fc <- fc[
    valid_symbol,
    ,
    drop = FALSE
  ]
  
  adjp <- adjp[
    valid_symbol,
    ,
    drop = FALSE
  ]
  
  ann <- ann[
    valid_symbol,
    ,
    drop = FALSE
  ]
  
  symbol <- symbol[
    valid_symbol
  ]
  
  mean_expression <- mean_expression[
    valid_symbol
  ]
  
  
  # ----------------------------------------------------------
  # Resolve duplicate gene symbols
  #
  # Keep the Ensembl/gene_id entry having the highest
  # mean expression across the experiment.
  # ----------------------------------------------------------
  
  ord <- order(
    symbol,
    -mean_expression
  )
  
  fc <- fc[ord, , drop = FALSE]
  adjp <- adjp[ord, , drop = FALSE]
  ann <- ann[ord, , drop = FALSE]
  symbol <- symbol[ord]
  mean_expression <- mean_expression[ord]
  
  
  keep_symbol <- !duplicated(symbol)
  
  fc <- fc[
    keep_symbol,
    ,
    drop = FALSE
  ]
  
  adjp <- adjp[
    keep_symbol,
    ,
    drop = FALSE
  ]
  
  ann <- ann[
    keep_symbol,
    ,
    drop = FALSE
  ]
  
  symbol <- symbol[
    keep_symbol
  ]
  
  
  rownames(fc) <- symbol
  rownames(adjp) <- symbol
  
  
  # ----------------------------------------------------------
  # Cleaner column names, compatible with old DepSHOCK object
  # ----------------------------------------------------------
  
  contrast_info <- do.call(
    rbind,
    contrast_meta
  )
  
  time_label <- ifelse(
    contrast_info$time_h == 24,
    "24h",
    paste0(
      "day",
      contrast_info$time_h / 24
    )
  )
  
  
  final_names <- paste0(
    contrast_info$guide,
    "_",
    time_label
  )
  
  
  colnames(fc) <- final_names
  colnames(adjp) <- final_names
  
  contrast_info$experiment_id <- final_names
  
  
  # ----------------------------------------------------------
  # Final QC
  # ----------------------------------------------------------
  
  stopifnot(
    identical(
      rownames(fc),
      rownames(adjp)
    ),
    identical(
      colnames(fc),
      colnames(adjp)
    )
  )
  
  
  cat(
    "Final dimensions:",
    nrow(fc),
    "genes x",
    ncol(fc),
    "contrasts\n"
  )
  
  cat(
    "NA in logFC matrix:",
    round(100 * mean(is.na(fc)), 4),
    "%\n"
  )
  
  cat(
    "NA in adjusted-P matrix:",
    round(100 * mean(is.na(adjp)), 4),
    "%\n"
  )
  
  
  return(
    list(
      fc = fc,
      adjp = adjp,
      contrast_metadata = contrast_info,
      sample_metadata = meta,
      gene_annotation = ann,
      voom = v,
      fit = fit2
    )
  )
}


# ============================================================
# 4. Complete wrapper
# ============================================================

DS_build_DE_matrices <- function(
    sgRNA_id_path,
    id_map_path,
    rawCount_path,
    cell_line,
    min_replicates = 2
) {
  
  cat("Building sample map...\n")
  
  sample_map <- DS_build_sample_map(
    sgRNA_id_path = sgRNA_id_path,
    id_map_path = id_map_path,
    cell_line = cell_line
  )
  
  
  cat(
    "Mapped libraries:",
    nrow(sample_map),
    "\n"
  )
  
  
  cat("Reading Salmon counts...\n")
  
  counts_obj <- DS_read_salmon_counts(
    rawCount_path = rawCount_path,
    sample_map = sample_map
  )
  
  
  cat(
    "Count matrix:",
    nrow(counts_obj$counts),
    "genes x",
    ncol(counts_obj$counts),
    "samples\n"
  )
  
  
  cat("Running differential expression...\n")
  
  res <- DS_run_differential_expression(
    counts = counts_obj$counts,
    gene_annotation = counts_obj$genes,
    sample_map = sample_map,
    min_replicates = min_replicates
  )
  
  
  return(res)
}


DS_create_analysis_matrix <- function(fc, adjp, padj_thr = 0.05) {
  
  stopifnot(
    identical(rownames(fc), rownames(adjp)),
    identical(colnames(fc), colnames(adjp))
  )
  
  ## Keep genes significantly DE in >= 1 experiment
  sig_gene <- rowSums(
    !is.na(adjp) & adjp < padj_thr
  ) >= 1
  
  X <- fc[sig_gene, , drop = FALSE]
  
  ## Remove genes containing NA / Inf or with zero variance
  good_gene <- apply(
    X,
    1,
    function(x) {
      all(is.finite(x)) &&
        sd(x) > 0
    }
  )
  
  X <- X[good_gene, , drop = FALSE]
  
  cat("Genes significant in >=1 contrast:", sum(sig_gene), "\n")
  cat("Genes retained for PCA:", nrow(X), "\n")
  
  ## PCA needs experiments on rows, genes on columns
  X <- t(X)
  
  return(X)
}

DS_plot_DE_PCA <- function(
    DE_obj,
    pc1 = 1,
    pc2 = 2,
    padj_thr = 0.05,
    plot_scree = TRUE,
    color_by = c("perturbation", "time"),
    cex = 1.3
) {
  
  ## ========================================================
  ## Extract matrices
  ## ========================================================
  
  fc   <- DE_obj$fc
  adjp <- DE_obj$adjp
  meta <- DE_obj$contrast_metadata
  
  stopifnot(
    identical(rownames(fc), rownames(adjp)),
    identical(colnames(fc), colnames(adjp))
  )
  
  
  ## ========================================================
  ## Select informative genes
  ## ========================================================
  
  sig_gene <- rowSums(
    !is.na(adjp) & adjp < padj_thr
  ) >= 1
  
  X <- fc[
    sig_gene,
    ,
    drop = FALSE
  ]
  
  
  ## Remove non-finite / invariant genes
  good_gene <- apply(
    X,
    1,
    function(x) {
      all(is.finite(x)) &&
        sd(x) > 0
    }
  )
  
  X <- X[
    good_gene,
    ,
    drop = FALSE
  ]
  
  
  cat(
    "Genes significant in >=1 contrast:",
    sum(sig_gene),
    "\n"
  )
  
  cat(
    "Genes retained for PCA:",
    nrow(X),
    "\n"
  )
  
  
  ## ========================================================
  ## Experiments on rows, genes on columns
  ## ========================================================
  
  X <- t(X)
  
  
  ## ========================================================
  ## PCA
  ## ========================================================
  
  pca <- prcomp(
    X,
    center = TRUE,
    scale. = FALSE
  )
  
  var_exp <- 100 *
    pca$sdev^2 /
    sum(pca$sdev^2)
  
  
  cat(
    "Variance explained:\n"
  )
  
  print(
    round(
      var_exp[1:min(10, length(var_exp))],
      2
    )
  )
  
  
  ## ========================================================
  ## Align metadata
  ## ========================================================
  
  if (!"experiment_id" %in% colnames(meta)) {
    stop(
      "contrast_metadata must contain an 'experiment_id' column."
    )
  }
  
  meta <- meta[
    match(
      rownames(pca$x),
      meta$experiment_id
    ),
    ,
    drop = FALSE
  ]
  
  if (anyNA(meta$experiment_id)) {
    stop(
      "Could not align PCA samples with contrast metadata."
    )
  }
  
  stopifnot(
    identical(
      rownames(pca$x),
      meta$experiment_id
    )
  )
  
  
  ## ========================================================
  ## Set layout
  ## ========================================================
  
  nplots <- length(color_by) +
    as.integer(plot_scree)
  
  oldpar <- par(
    mfrow = c(
      1,
      nplots
    )
  )
  
  on.exit(
    par(oldpar)
  )
  
  
  ## ========================================================
  ## Scree plot
  ## ========================================================
  
  if (plot_scree) {
    
    n_pc <- min(
      20,
      length(var_exp)
    )
    
    plot(
      seq_len(n_pc),
      var_exp[seq_len(n_pc)],
      type = "b",
      pch = 16,
      xlab = "Principal component",
      ylab = "Variance explained (%)",
      bty = "l"
    )
  }
  
  
  ## ========================================================
  ## PCA plots
  ## ========================================================
  
  for (cb in color_by) {
    
    if (cb == "perturbation") {
      
      group <- ifelse(
        meta$target == "RPL15",
        "RPL15",
        ifelse(
          meta$target == "AAVS1",
          "AAVS1",
          ifelse(
            meta$target == "NT",
            "Non-targeting",
            "Oncogene"
          )
        )
      )
      
      
      cols <- c(
        "Oncogene"      = "#2166AC",
        "RPL15"         = "#B2182B",
        "AAVS1"         = "#4D9221",
        "Non-targeting" = "grey55"
      )
      
      
      point_cols <- cols[
        group
      ]
      
      
      plot(
        pca$x[, pc1],
        pca$x[, pc2],
        pch = 16,
        cex = cex,
        col = point_cols,
        xlab = paste0(
          "PC",
          pc1,
          " (",
          round(var_exp[pc1], 1),
          "%)"
        ),
        ylab = paste0(
          "PC",
          pc2,
          " (",
          round(var_exp[pc2], 1),
          "%)"
        ),
        main = "Perturbation",
        bty = "l"
      )
      
      
      legend(
        "bottomright",
        legend = names(cols),
        col = cols,
        pch = 16,
        bty = "n",
        cex = 0.9
      )
    }
    
    
    if (cb == "time") {
      
      times <- sort(
        unique(
          meta$time_h
        )
      )
      
      time_cols <- setNames(
        rev(
          hcl.colors(
            length(times),
            palette = "Viridis"
          )
        ),
        times
      )
      
      
      point_cols <- time_cols[
        as.character(
          meta$time_h
        )
      ]
      
      
      plot(
        pca$x[, pc1],
        pca$x[, pc2],
        pch = 16,
        cex = cex,
        col = point_cols,
        xlab = paste0(
          "PC",
          pc1,
          " (",
          round(var_exp[pc1], 1),
          "%)"
        ),
        ylab = paste0(
          "PC",
          pc2,
          " (",
          round(var_exp[pc2], 1),
          "%)"
        ),
        main = "Time",
        bty = "l"
      )
      
      
      legend(
        "bottomright",
        legend = paste0(
          "Day ",
          times / 24
        ),
        col = time_cols,
        pch = 16,
        bty = "n",
        cex = 0.9
      )
    }
  }
  
  
  ## ========================================================
  ## Return useful objects
  ## ========================================================
  
  invisible(
    list(
      pca = pca,
      X = X,
      metadata = meta,
      variance_explained = var_exp,
      genes_used = rownames(
        fc
      )[sig_gene][good_gene]
    )
  )
}