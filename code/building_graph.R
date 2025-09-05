library(Matrix)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(pheatmap)
library(Rtsne)
library(flashier)
library(gbcd)
library(tidyverse)
library(dbplyr)

data(hnscc)
dim(hnscc$Y)
head(hnscc$info)

setwd("~/Documents/Prof Dat/Graph-MF")
# L = loadRDS(here("data","cancer_L.rds"))
L = readRDS('data/cancer_L.rds')
dim(L)


anno <- data.frame(sample = hnscc$info$sample, subtype = hnscc$info$subtype)
rownames(anno) <- rownames(L)
anno_colors <- list(sample = hnscc$sample_col, subtype = hnscc$subtype_col)
cols <- colorRampPalette(c("gray96", "red"))(50)
brks <- seq(0, 1, 0.02)

pheatmap(L[order(anno$sample), -c(1)], cluster_rows = FALSE,
         cluster_cols = FALSE, show_rownames = FALSE, annotation_row = anno,
         annotation_colors = anno_colors, annotation_names_row = FALSE,
         angle_col = 45, fontsize = 9, color = cols, breaks = brks,
         main = "")

V = (L>2e-05)*1
V_filter = V[,1+c(1,2,3,7,8,9,11,12,15,19,20)]
motifs = unique(V_filter)
counts <- as.data.frame(V_filter) %>%
     group_by(across(everything())) %>%
     summarise(count = n(), .groups = "drop") %>%
     mutate(freq = count / sum(count) * 100)
b = counts$freq

# keep the most frequent rows until their cumulative frequency explains 95% of the data
filtered_counts = counts %>%arrange(desc(freq)) %>%
     mutate(cum_freq = cumsum(freq)) %>%
     filter(cum_freq <= 95)

#####
df11 <- filtered_counts[, -( (ncol(filtered_counts)-2):ncol(filtered_counts) )]   # keep only the first 11 columns

# convert to matrix for speed
mat <- as.matrix(df11)

# number of columns
p <- ncol(mat)

# initialize result matrix
dmat <- matrix(0, nrow = p, ncol = p,
               dimnames = list(colnames(mat), colnames(mat)))

# compute d(i,j)
for (i in 1:p) {
     idx_i <- which(mat[, i] != 0)               # nonzero indices in column i
     denom <- length(idx_i)                      # denominator
     if (denom > 0) {
          for (j in 1:p) {
               idx_j <- which(mat[, j] != 0)
               numer <- length(intersect(idx_i, idx_j)) # common nonzero indices
               dmat[i, j] <- numer / denom
          }
     } else {
          dmat[i, ] <- NA  # or 0 if you prefer
     }
}
library(igraph)
g <- graph_from_adjacency_matrix(dmat,
                                 mode = "undirected", 
                                 weighted = TRUE,
                                 diag = FALSE)

## Function to generate semicircle coordinates
semi_circle <- function(k, radius = 1, translate = 0, top = TRUE) {
     # angles equally spaced from 0 to pi (top) or pi to 2pi (bottom)
     if (top) {
          angles <- seq(0, pi, length.out = k + 2)[-c(1, k + 2)] # exclude ends
     } else {
          angles <- seq(pi, 2*pi, length.out = k + 2)[-c(1, k + 2)]
     }
     x <- radius * cos(angles) 
     y <- radius * sin(angles) + translate
     cbind(x, y)
}
# c(1,2,3,7,8,9,11,12,15,19,20)
# GEP_top = c(4, 5, 6, 10, 13, 14, 16, 17, 18)
# GEP_mid = c(1, 7, 2, 3)
col_sums <- colSums(dmat)
# compute tertile cutoffs (33% and 67%)
qs <- quantile(col_sums, probs = c(0.25, 0.5, 0.75))

GEP_top = which(col_sums>qs[3])
# GEP_mid = which(col_sums >= qs[2] & col_sums< qs[3])
GEP_bot = setdiff(c(1:11), GEP_top)
# Generate coordinates
coords_top <- semi_circle(length(GEP_top), radius = 5, translate = 10, top = TRUE)
coords_mid <- semi_circle(length(GEP_mid), radius = 5, top = TRUE)
coords_bottom <- semi_circle(length(GEP_bot), radius = 5, translate = -5, top = FALSE)

layout_coords <- matrix(0, nrow = 11, ncol = 2)
layout_coords[GEP_top, ] = coords_top
# layout_coords[GEP_mid, ] = coords_mid
layout_coords[GEP_bot, ] = coords_bottom

plot(g,
     layout = layout_coords,
     vertex.size = 15,
     vertex.label.cex = 1.2,
     vertex.color = "skyblue",
     edge.color = "gray40",
     edge.width = E(g)$weight*5,   # scale edge thickness by weight
     layout = layout_with_fr)
