# remotes::install_github('stephenslab/gbcd')

library(Matrix)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(pheatmap)
library(Rtsne)
library(flashier)
library(gbcd)

data(hnscc)
dim(hnscc$Y)
head(hnscc$info)

## t-sne plot of the data
set.seed(1)
cols <- order(apply(hnscc$Y,2,sd),decreasing = TRUE)
cols <- cols[1:3000]
res  <- Rtsne(as.matrix(hnscc$Y[,cols]),normalize = TRUE)
colnames(res$Y) <- c("tsne1","tsne2")
pdat <- cbind(res$Y,hnscc$info)
ggplot(pdat,aes(x = tsne1,y = tsne2,color = sample,shape = subtype)) +
  geom_point(size = 1) + scale_color_manual(values = hnscc$sample_col) +
  scale_shape_manual(values = c(15,16,18)) +
  theme_cowplot(font_size = 10)


## gbcd analysis
res.gbcd <- fit_gbcd(Y = hnscc$Y, Kmax = 12, maxiter1 = 100,
                     maxiter2 = 50, maxiter3 = 50, 
                     prior = flash_ebnm(prior_family = "generalized_binary",
                                        scale = 0.04))

head(round(res.gbcd$L, 3))

## plot

anno <- data.frame(sample = hnscc$info$sample, subtype = hnscc$info$subtype)
rownames(anno) <- rownames(res.gbcd$L)
anno_colors <- list(sample = hnscc$sample_col, subtype = hnscc$subtype_col)
cols <- colorRampPalette(c("gray96", "red"))(50)
brks <- seq(0, 1, 0.02)

pheatmap(res.gbcd$L[order(anno$sample), -c(1)], cluster_rows = FALSE,
         cluster_cols = FALSE, show_rownames = FALSE, annotation_row = anno,
         annotation_colors = anno_colors, annotation_names_row = FALSE,
         angle_col = 45, fontsize = 9, color = cols, breaks = brks,
         main = "")



L = res.gbcd$L

setwd("~/Documents/single-cell-jamboree/visualization-graph")
saveRDS(L, file='cancer_L.rds')


dim(L)

#############
# Visualization
L = loadRDS('cancer_L.rds')
dim(L)
L = L[, c(2:21, 1)] ## put baseline factor on the bottom

L = res.gbcd$L[, c(2:21)]

library(igraph)
A = t(L) %*% L
diag(A) = 0
A = A / max(A) 
rownames(A) <- c(1:20)
colnames(A) <- c(1:20)

q_trunc = quantile(as.vector(A), prob=0.5)
A_trunc = A * (A > q_trunc) * 3

g <- graph_from_adjacency_matrix(A_trunc, mode = "undirected", weighted=TRUE)

# Function to generate semicircle coordinates
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

GEP_top = c(4, 5, 6, 10, 13, 14, 16, 17, 18)
GEP_mid = c(1, 11, 2, 3)
GEP_bot = c(7, 8, 9, 12, 15, 19, 20)
# Generate coordinates
coords_top <- semi_circle(length(GEP_top), radius = 5, translate = 10, top = TRUE)
coords_mid <- semi_circle(length(GEP_mid), radius = 5, top = TRUE)
coords_bottom <- semi_circle(length(GEP_bot), radius = 5, translate = -5, top = FALSE)

layout_coords <- matrix(0, nrow = 20, ncol = 2)
layout_coords[GEP_top, ] = coords_top
layout_coords[GEP_mid, ] = coords_mid
layout_coords[GEP_bot, ] = coords_bottom

plot(g,
     layout = layout_coords,
     vertex.size = 15,
     vertex.label.cex = 1.2,
     vertex.color = "skyblue",
     edge.color = "gray40",
     edge.width = E(g)$weight,   # scale edge thickness by weight
     layout = layout_with_fr)

## only mid and bottom
sub_indx = c(GEP_mid, GEP_bot) 
A_sub = A_trunc[sub_indx, sub_indx]
sub_layout = layout_coords[sub_indx, ]
g <- graph_from_adjacency_matrix(A_sub, mode = "undirected", weighted=TRUE)

plot(g,
     layout = sub_layout,
     vertex.size = 15,
     vertex.label.cex = 1.2,
     vertex.color = "skyblue",
     edge.color = "gray40",
     edge.width = E(g)$weight,   # scale edge thickness by weight
     layout = layout_with_fr)



