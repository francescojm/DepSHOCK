DepSHOCK_obj<-readRDS('/Volumes/iorio/raw_data_inHouse/DepSHOCK/bulkRNAseq/depshock_dataset.rds')

library(dplyr)
library(ggplot2)

## ============================================================
## 1. Extract the three components of the DepSHOCK object
## ============================================================

lfc  <- DepSHOCK_obj$logFC
padj <- DepSHOCK_obj$padj
meta <- DepSHOCK_obj$metadata


## ============================================================
## 2. Make sure metadata and matrices are aligned
## ============================================================

meta <- meta[
  match(rownames(lfc), meta$guide_id),
]

stopifnot(
  identical(rownames(lfc), rownames(padj)),
  identical(colnames(lfc), colnames(padj)),
  identical(rownames(lfc), meta$guide_id)
)


## ============================================================
## 3. Calculate basic transcriptional-response metrics
## ============================================================

## Number of DE genes, padj < 0.05
meta$n_DEG <- rowSums(
  padj < 0.05,
  na.rm = TRUE
)

## Number of significantly upregulated genes
meta$n_up <- rowSums(
  padj < 0.05 & lfc > 0,
  na.rm = TRUE
)

## Number of significantly downregulated genes
meta$n_down <- rowSums(
  padj < 0.05 & lfc < 0,
  na.rm = TRUE
)

## Overall magnitude of the transcriptional response:
## Euclidean norm of the logFC vector
meta$response_norm <- sqrt(
  rowSums(lfc^2, na.rm = TRUE)
)


## ============================================================
## 4. Define perturbation classes
## ============================================================

## Change "NT" here if the target has a different label
control_targets <- c("AAVS1", "RPL15", "NT")

meta$type <- ifelse(
  meta$target == "RPL15", "RPL15",
  ifelse(
    meta$target == "AAVS1", "AAVS1",
    ifelse(
      meta$target == "NT", "Non-targeting",
      "Oncogene"
    )
  )
)


## Quick check
table(meta$target, meta$type)

## ============================================================
## 5. Colours
## ============================================================

cols <- c(
  "Oncogene"      = "#2166AC",
  "RPL15"         = "#B2182B",
  "AAVS1"         = "#4D9221",
  "Non-targeting" = "grey55"
)


## ============================================================
## 6. Function for longitudinal plots
## ============================================================

plot_trajectories <- function(meta, variable, ylab) {

  cell_lines <- unique(meta$cell_line)

  n <- length(cell_lines)

  nc <- ceiling(sqrt(n))
  nr <- ceiling(n / nc)

  oldpar <- par(
    mfrow = c(nr, nc),
    mar = c(4, 4.5, 3, 1),
    oma = c(0, 0, 0, 0),
    bty = "n"
  )

  for (cl in cell_lines) {

    d <- meta[meta$cell_line == cl, ]

    oncogene <- unique(
      d$target[d$type == "Oncogene"]
    )

    oncogene <- paste(oncogene, collapse = "/")

    y <- d[[variable]]

    plot(
      NA,
      xlim = range(d$time_h, na.rm = TRUE),
      ylim = range(y, na.rm = TRUE),
      xlab = "Time after Cas9 activation (h)",
      ylab = ylab,
      main = paste0(cl, " — ", oncogene),
      xaxt = "n"
    )

    axis(
      1,
      at = sort(unique(d$time_h))
    )

    for (g in unique(d$guide)) {

      dg <- d[d$guide == g, ]
      dg <- dg[order(dg$time_h), ]

      type_g <- unique(dg$type)

      lines(
        dg$time_h,
        dg[[variable]],
        type = "b",
        pch = 16,
        lwd = 2,
        col = cols[type_g]
      )
    }

    box()
  }

  ## Shared legend
  par(xpd = NA)

  legend(
    "right",
    inset = c(-1, -0.05),
    legend = names(cols),
    col = cols,
    lwd = 2,
    pch = 16,
    horiz = FALSE,
    bty = "n",
    xpd = NA
  )

  par(oldpar)
}

plot_trajectories(
  meta,
  variable = "response_norm",
  ylab = "Transcriptional response magnitude"
)

plot_trajectories(
  meta,
  variable = "n_DEG",
  ylab = "Number of differentially expressed genes"
)

padj < 0.05

plot_trajectories(
  meta,
  variable = "n_up",
  ylab = "Number of upregulated genes"
)

plot_trajectories(
  meta,
  variable = "n_down",
  ylab = "Number of downregulated genes"
)
