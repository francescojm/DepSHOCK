
#DepSHOCK_obj<-readRDS('/Volumes/iorio/raw_data_inHouse/DepSHOCK/bulkRNAseq/depshock_dataset.rds')

## ============================================================
## DepSHOCK first-pass signature deconvolution
## ============================================================

## Functions
mean_signature <- function(ids, lfc) {
  
  if(length(ids) == 0)
    return(NULL)
  
  colMeans(
    lfc[ids, , drop = FALSE],
    na.rm = TRUE
  )
}


lfc  <- DepSHOCK_obj$logFC
padj  <- DepSHOCK_obj$padj

lfc  <- as.matrix(DepSHOCK_obj$logFC)
padj <- as.matrix(DepSHOCK_obj$padj)

library(AnnotationDbi)
library(org.Hs.eg.db)

## Remove Ensembl version suffixes if present
ensembl_ids <- sub("\\..*$", "", colnames(lfc))

## Map Ensembl IDs -> HGNC gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

gene_symbols <- unname(gene_symbols)

sum(!is.na(gene_symbols))
sum(is.na(gene_symbols))

head(
  data.frame(
    ensembl = ensembl_ids,
    symbol = gene_symbols
  ),
  20
)

gene_names <- gene_symbols

gene_names[
  is.na(gene_names) |
    gene_names == ""
] <- ensembl_ids[
  is.na(gene_names) |
    gene_names == ""
]

sum(duplicated(gene_names))

gene_names[duplicated(gene_names)]

gene_names <- make.unique(gene_names)

colnames(lfc)  <- gene_names
colnames(padj) <- gene_names

storage.mode(lfc)  <- "numeric"
storage.mode(padj) <- "numeric"

meta <- DepSHOCK_obj$metadata

## Make sure metadata are aligned
meta <- meta[match(rownames(lfc), meta$guide_id), ]

stopifnot(
  identical(rownames(lfc), meta$guide_id)
)

control_targets <- c("NT", "AAVS1", "RPL15")


deconv <- list()

summary_list <- list()

k <- 1

for(cl in unique(meta$cell_line)) {
  
  times <- sort(
    unique(meta$time_h[meta$cell_line == cl])
  )
  
  for(tt in times) {
    
    m <- meta[
      meta$cell_line == cl &
        meta$time_h == tt,
    ]
    
    ## ------------------------
    ## Control signatures
    ## ------------------------
    
    NT_ids <- m$guide_id[
      m$target == "NT"
    ]
    
    AAVS1_ids <- m$guide_id[
      m$target == "AAVS1"
    ]
    
    RPL15_ids <- m$guide_id[
      m$target == "RPL15"
    ]
    
    N <- mean_signature(NT_ids, lfc)
    A <- mean_signature(AAVS1_ids, lfc)
    R <- mean_signature(RPL15_ids, lfc)
    
    ## Skip if controls are missing
    if(
      is.null(N) |
      is.null(A) |
      is.null(R)
    ) next
    
    
    ## ------------------------
    ## Fundamental components
    ## ------------------------
    
    ## DSB/editing response
    E <- A - N
    
    ## RPL15-derived generic
    ## loss-of-fitness direction
    D <- R - A
    
    
    ## ------------------------
    ## Oncogene guides
    ## ------------------------
    
    onc_rows <- which(
      !(m$target %in% control_targets)
    )
    
    for(i in onc_rows) {
      
      guide_id <- m$guide_id[i]
      guide    <- m$guide[i]
      target   <- m$target[i]
      
      O <- lfc[guide_id, ]
      
      ## Oncogene response beyond AAVS1
      Q <- O - A
      
      
      ## ------------------------
      ## Projection onto generic
      ## loss-of-fitness axis
      ## ------------------------
      
      ok <- is.finite(Q) &
        is.finite(D)
      
      beta <- sum(
        Q[ok] * D[ok]
      ) /
        sum(
          D[ok]^2
        )
      
      
      ## Generic component actually
      ## present in oncogene response
      G <- beta * D
      
      
      ## Dependency-specific residual
      S <- Q - G
      
      
      ## ------------------------
      ## Cosine similarity
      ## ------------------------
      
      cosine <- sum(
        Q[ok] * D[ok]
      ) /
        (
          sqrt(sum(Q[ok]^2)) *
            sqrt(sum(D[ok]^2))
        )
      
      cosine <- max(
        -1,
        min(1, cosine)
      )
      
      angle <- acos(cosine) *
        180 / pi
      
      
      ## Fraction of squared response
      ## aligned with generic axis
      
      fraction_generic <- sum(
        G[ok]^2
      ) /
        sum(
          Q[ok]^2
        )
      
      
      ## ------------------------
      ## Store results
      ## ------------------------
      
      deconv[[guide_id]] <- list(
        
        cell_line = cl,
        target = target,
        guide = guide,
        time_h = tt,
        
        NT = N,
        
        editing = E,
        
        generic_axis = D,
        
        oncogene_excess = Q,
        
        generic_component = G,
        
        specific_component = S,
        
        beta = beta,
        
        cosine = cosine,
        
        angle = angle,
        
        fraction_generic = fraction_generic
      )
      
      
      summary_list[[k]] <- data.frame(
        
        guide_id = guide_id,
        cell_line = cl,
        target = target,
        guide = guide,
        time_h = tt,
        
        beta = beta,
        cosine = cosine,
        angle = angle,
        fraction_generic = fraction_generic
      )
      
      k <- k + 1
    }
  }
}

deconv_summary <- do.call(
  rbind,
  summary_list
)

rownames(deconv_summary) <- NULL

guide_benchmark <- data.frame()

for(cl in unique(deconv_summary$cell_line)) {
  
  targets <- unique(
    deconv_summary$target[
      deconv_summary$cell_line == cl
    ]
  )
  
  for(target in targets) {
    
    times <- sort(unique(
      deconv_summary$time_h[
        deconv_summary$cell_line == cl &
          deconv_summary$target == target
      ]
    ))
    
    for(tt in times) {
      
      d <- deconv_summary[
        deconv_summary$cell_line == cl &
          deconv_summary$target == target &
          deconv_summary$time_h == tt,
      ]
      
      if(nrow(d) != 2)
        next
      
      id1 <- d$guide_id[1]
      id2 <- d$guide_id[2]
      
      Q1 <- deconv[[id1]]$oncogene_excess
      Q2 <- deconv[[id2]]$oncogene_excess
      
      S1 <- deconv[[id1]]$specific_component
      S2 <- deconv[[id2]]$specific_component
      
      G1 <- deconv[[id1]]$generic_component
      G2 <- deconv[[id2]]$generic_component
      
      guide_benchmark <- rbind(
        guide_benchmark,
        data.frame(
          cell_line = cl,
          target = target,
          time_h = tt,
          
          cor_Q = cor(
            Q1, Q2,
            use = "pairwise.complete.obs",
            method = "pearson"
          ),
          
          cor_S = cor(
            S1, S2,
            use = "pairwise.complete.obs",
            method = "pearson"
          ),
          
          cor_G = cor(
            G1, G2,
            use = "pairwise.complete.obs",
            method = "pearson"
          )
        )
      )
    }
  }
}


cell_lines <- unique(deconv_summary$cell_line)

n <- length(cell_lines)
nc <- ceiling(sqrt(n))
nr <- ceiling(n / nc)

oldpar <- par(
  mfrow = c(nr, nc),
  mar = c(4, 4.5, 3, 1)
)

for(cl in cell_lines) {
  
  d <- deconv_summary[
    deconv_summary$cell_line == cl,
  ]
  
  oncogene <- paste(unique(d$target), collapse = "/")
  
  plot(
    NA,
    xlim = range(d$time_h),
    ylim = c(-1, 1),
    xlab = "Time after Cas9 activation (h)",
    ylab = "Cosine similarity with RPL15 axis",
    main = paste0(cl, " — ", oncogene),
    xaxt = "n",
    bty = "l"
  )
  
  axis(
    1,
    at = sort(unique(d$time_h))
  )
  
  abline(
    h = 0,
    lty = 2,
    col = "grey70"
  )
  
  for(g in unique(d$guide)) {
    
    dg <- d[d$guide == g, ]
    dg <- dg[order(dg$time_h), ]
    
    lines(
      dg$time_h,
      dg$cosine,
      type = "b",
      pch = 16,
      lwd = 2
    )
  }
}

par(oldpar)


oldpar <- par(
  mfrow = c(nr, nc),
  mar = c(4, 4.5, 3, 1)
)

for(cl in cell_lines) {
  
  d <- deconv_summary[
    deconv_summary$cell_line == cl,
  ]
  
  oncogene <- paste(unique(d$target), collapse = "/")
  
  plot(
    NA,
    xlim = range(d$time_h),
    ylim = range(deconv_summary$beta, na.rm = TRUE),
    xlab = "Time after Cas9 activation (h)",
    ylab = "RPL15-axis coefficient (beta)",
    main = paste0(cl, " — ", oncogene),
    xaxt = "n",
    bty = "l"
  )
  
  axis(
    1,
    at = sort(unique(d$time_h))
  )
  
  abline(
    h = 0,
    lty = 2,
    col = "grey70"
  )
  
  for(g in unique(d$guide)) {
    
    dg <- d[d$guide == g, ]
    dg <- dg[order(dg$time_h), ]
    
    lines(
      dg$time_h,
      dg$beta,
      type = "b",
      pch = 16,
      lwd = 2
    )
  }
}

par(oldpar)


braf_ids <- deconv_summary$guide_id[
  deconv_summary$target == "BRAF" &
    deconv_summary$time_h == 72
]

braf_specific <- do.call(
  rbind,
  lapply(
    braf_ids,
    function(id) deconv[[id]]$specific_component
  )
)

braf_specific_mean <- colMeans(
  braf_specific,
  na.rm = TRUE
)


braf_down <- sort(
  braf_specific_mean[braf_specific_mean < 0]
)

head(braf_down, 30)

