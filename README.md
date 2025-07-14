# Unsupervised Learning for RNA-seq Breast Cancer Data

This repository contains an analysis of RNA-seq gene expression data from breast cancer samples using unsupervised learning techniques. The goal was to identify patterns and latent structures in the data through dimensionality reduction methods.

## Project Structure

- `data/`: Contains the input data files (`gene_expression.csv`, `column_names.txt`, `classes.csv`).
- `scripts/unsupervised_analysis.R`: Fully commented R script implementing all four unsupervised learning methods.
- `README.md`: Project documentation.

## Objective

To compare the performance of four unsupervised dimensionality reduction techniques in capturing relevant patterns within high-dimensional biological data.

## Methods Implemented

1. **Principal Component Analysis (PCA)**
2. **Multidimensional Scaling (MDS)**
3. **Isomap**
4. **t-distributed Stochastic Neighbor Embedding (t-SNE)**

Each technique was applied in 2 and 3 dimensions. Performance was evaluated using:
- Variance explained (PCA)
- Stress value (MDS, Isomap)
- Neighbor conservation rate (t-SNE)


## Results Summary

- **PCA**:  
  The first two principal components explain **21.82%** of the variance, and the first three explain **30.52%**. While some variability is captured, important information may be lost when reducing to 2D or 3D. Most classes overlap, but samples from class **AGH** show a clearer separation.

- **MDS**:  
  Stress values are **0.566 (2D)** and **0.467 (3D)**, indicating moderate distortion. The AGH class is somewhat distinguishable, although other classes remain entangled. MDS performs slightly better than PCA in preserving global structure.

- **Isomap**:  
  High stress values (**2.32 in 2D**, **2.60 in 3D**) suggest poor distance preservation. However, visual separation of classes is more apparent, with some clusters clearly corresponding to class labels.

- **t-SNE**:  
  The 10-nearest neighbor preservation rate is **47.17% (2D)** and **50.96% (3D)**. While not optimal, t-SNE provides the most visually distinct class separations. However, its stochastic nature and lower preservation of local structures suggest limited reproducibility.


## Biological Interpretation

Although no method perfectly preserved the data structure, **t-SNE and MDS** revealed biologically relevant separations—particularly for class **AGH**, which consistently appeared isolated across methods. This suggests potential transcriptional uniqueness in that group.

The modest performance of PCA and Isomap indicates the complexity and non-linear nature of RNA-seq data, highlighting the importance of method selection in exploratory transcriptomic analyses.


## 👨‍💻 Author

**Oscar Saiz Gutierrez**  
MSc in Bioinformatics  

---

**Note:** This project was developed as part of the academic course *Algorithms and Artificial Intelligence* in the MSc in Bioinformatics.
