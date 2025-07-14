# UNSUPERVISED ANALYSIS OF RNA-seq BREAST CANCER DATA

# Install required packages
install.packages(c("readr", "caret", "ggplot2", "plotly", "vegan", "Rtsne", "FNN"))
# Load the libraries
library(readr)
library(caret)
library(ggplot2)
library(plotly)
library(stats)
library(vegan)
library(Rtsne)
library(FNN)

# Load data
gene_expression <- read_delim("data/gene_expression.csv", 
                              delim = ";", escape_double = FALSE, 
                              col_names = read_lines("data/column_names.txt"), trim_ws = TRUE)

classes <- read_delim("data/classes.csv", delim = ";", 
                      escape_double = FALSE, col_names = c("Sample", "Class"), trim_ws = TRUE)

# Remove columns with zero variance
sum_vector <- colSums(gene_expression)
zero_columns <- names(sum_vector[sum_vector == 0])
print(zero_columns)
data <- gene_expression[, !names(gene_expression) %in% zero_columns]
data$Class <- classes$Class

# --- PCA ---
data_pca <- data[,-ncol(data)]
pca_result <- prcomp(data_pca, scale. = TRUE)
print(summary(pca_result))
cumulative_proportion <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
print(cumulative_proportion)

data_2d <- data.frame(pca_result$x[,1:2], Class = data$Class)
data_3d <- data.frame(pca_result$x[,1:3], Class = data$Class)

ggplot(data_2d, aes(x = PC1, y = PC2, color = Class)) +
  geom_point() +
  ggtitle("PCA - 2 Principal Components")

plot_ly(data_3d, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Class) %>%
  add_markers() %>%
  layout(title = 'PCA - 3 Principal Components')

# --- MDS ---
data_mds <- data[,-ncol(data)]
dist_matrix <- dist(data_mds)
mds_2d <- cmdscale(dist_matrix, k = 2, eig = TRUE)
mds_3d <- cmdscale(dist_matrix, k = 3, eig = TRUE)

data_2d <- data.frame(mds_2d$points, Class = data$Class)
data_3d <- data.frame(mds_3d$points, Class = data$Class)

ggplot(data_2d, aes(x = V1, y = V2, color = Class)) +
  geom_point() +
  ggtitle("MDS - 2 Dimensions")

plot_ly(data_3d, x = ~V1, y = ~V2, z = ~V3, color = ~Class) %>%
  add_markers() %>%
  layout(title = 'MDS - 3 Dimensions')

stress_2d <- sqrt(sum((dist_matrix - dist(mds_2d$points))^2) / sum(dist_matrix^2))
stress_3d <- sqrt(sum((dist_matrix - dist(mds_3d$points))^2) / sum(dist_matrix^2))
print(paste("MDS Stress 2D:", stress_2d))
print(paste("MDS Stress 3D:", stress_3d))

# --- Isomap ---
data_iso <- data[,-ncol(data)]
dist_matrix <- dist(data_iso)
iso_2d <- isomap(dist_matrix, k = 10, ndim = 2)
iso_3d <- isomap(dist_matrix, k = 10, ndim = 3)

stress_2d <- sqrt(sum((dist_matrix - dist(iso_2d$points))^2) / sum(dist_matrix^2))
stress_3d <- sqrt(sum((dist_matrix - dist(iso_3d$points))^2) / sum(dist_matrix^2))
print(paste("Isomap Stress 2D:", stress_2d))
print(paste("Isomap Stress 3D:", stress_3d))

data_2d <- data.frame(iso_2d$points, Class = data$Class)
data_3d <- data.frame(iso_3d$points, Class = data$Class)

ggplot(data_2d, aes(x = Dim1, y = Dim2, color = Class)) +
  geom_point() +
  ggtitle("Isomap - 2 Dimensions")

plot_ly(data_3d, x = ~Dim1, y = ~Dim2, z = ~Dim3, color = ~Class) %>%
  add_markers() %>%
  layout(title = 'Isomap - 3 Dimensions')

# --- t-SNE ---
data_tsne <- data[,-ncol(data)]
set.seed(30)
tsne_2d <- Rtsne(data_tsne, dims = 2, perplexity = 30, verbose = FALSE)
tsne_3d <- Rtsne(data_tsne, dims = 3, perplexity = 30, verbose = FALSE)

conservation_rate <- function(original, reduced, k) {
  orig_nn <- get.knnx(data = original, query = original, k = k)
  red_nn <- get.knnx(data = reduced, query = reduced, k = k)
  overlap <- sapply(1:nrow(original), function(i) {
    length(intersect(orig_nn$nn.index[i,], red_nn$nn.index[i,]))
  })
  mean(overlap) / k
}

k <- 10
rate_2d <- conservation_rate(data_tsne, tsne_2d$Y, k)
rate_3d <- conservation_rate(data_tsne, tsne_3d$Y, k)
print(paste("t-SNE 2D Neighbor Conservation:", rate_2d))
print(paste("t-SNE 3D Neighbor Conservation:", rate_3d))

data_2d <- data.frame(tsne_2d$Y, Class = data$Class)
data_3d <- data.frame(tsne_3d$Y, Class = data$Class)

ggplot(data_2d, aes(x = V1, y = V2, color = Class)) +
  geom_point(alpha = 0.7) +
  ggtitle("t-SNE - 2 Dimensions")

plot_ly(data_3d, x = ~V1, y = ~V2, z = ~V3, color = ~Class) %>%
  add_markers() %>%
  layout(title = 't-SNE - 3 Dimensions')
