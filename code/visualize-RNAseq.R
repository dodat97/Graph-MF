library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

set.seed(1)

subsample_cell_types <- function (x, n = 1000) {
  cells <- NULL
  groups <- levels(x)
  for (g in groups) {
    i  <-  which(x == g)
    n0 <- min(n,length(i))
    i  <- sample(i,n0)
    cells <- c(cells,i)
  }
  return(sort(cells))
}

setwd("~/Documents/Graph-MF")
load("data/pancreas.RData")
load("data/pancreas_factors.RData")
timings0 <- timings
load("data/pancreas_factors2.RData")
timings <- c(timings0,timings)

cells <- subsample_cell_types(sample_info$celltype,n = 500)
L <- fl_nmf_ldf$L
k <- ncol(L)
colnames(L) <- paste0("k",1:k)
batch_factors    <- c(2:5,7:8,20)
celltype_factors <- c(6,11:19,21)
other_factors    <- c(1,9:10,22:23)
celltype <- sample_info$celltype
celltype <-
  factor(celltype,
         c("acinar","ductal","activated_stellate","quiescent_stellate",
           "endothelial","macrophage","mast","schwann","t_cell","alpha",
           "beta","delta","gamma","epsilon"))
p1 <- structure_plot(L,topics = batch_factors,grouping = sample_info$tech,
                     gap = 3,perplexity = 70,n = 500) +
  labs(y = "membership",title = "data-set factors",
       fill = "factor",color = "factor")
p2 <- structure_plot(L[cells,],topics = celltype_factors,
                     grouping = celltype[cells],gap = 3,
                     n = 500,perplexity = 70) +
  labs(y = "membership",title = "cell-type factors",
       fill = "factor",color = "factor")
p3 <- structure_plot(L[cells,],topics = other_factors,
                     grouping = celltype[cells],gap = 3,
                     n = 500,perplexity = 70) +
  labs(y = "membership",title = "other factors",
       fill = "factor",color = "factor")
plot_grid(p1,p2,p3,nrow = 3,ncol = 1)

structure_plot(L[cells,],topics = c(1:23),
               grouping = celltype[cells],gap = 3,
               n = 500,perplexity = 70) +
  labs(y = "membership",title = "cell-type factors",
       fill = "factor",color = "factor")


## plot the graph
L_norm = L / rowSums(L)
A = t(L_norm ) %*% L_norm 
diag(A) = 0
A = A / max(A) * 8
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted=TRUE)

colors <- brewer.pal(23, "Set1")   # qualitative, good for categories
colors

library(Polychrome)
data(glasbey)
# colors <- glasbey()[-1]
colors <- glasbey

class1 <- c(1:5, 7:10, 20, 22, 23)
class2 <- c(6, 11:19, 21)

# pos1 = 

# Plot it
plot(g,
     vertex.size = 30,
     vertex.label.cex = 1.2,
     vertex.color = colors,
     edge.color = "gray40",
     edge.width = E(g)$weight,   # scale edge thickness by weight
     layout = layout_with_fr)

