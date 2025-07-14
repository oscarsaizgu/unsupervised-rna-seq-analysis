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


## 📊 Results Summary

- **PCA** efficiently explained variance but failed to separate classes with nonlinear relationships.
- **MDS** revealed moderate separation based on distance preservation but suffered from stress.
- **Isomap** improved class separation when using an optimal number of neighbors (`k = 10`).
- **t-SNE** showed the best visual clustering and neighbor preservation, especially in 3D.

The results suggest that nonlinear methods like **t-SNE** and **Isomap** are more suitable for capturing biological structures in RNA-seq datasets.

## 🧠 Biological Interpretation

Some of the observed groupings in 2D/3D space corresponded with known biological subtypes in the `Class` labels, supporting the utility of unsupervised learning in exploratory transcriptomics.

## 👨‍💻 Author

**Oscar Saiz Gutierrez**  
MSc in Bioinformatics  

---

**Note:** This project was developed as part of the academic course *Algorithms and Artificial Intelligence* in the MSc in Bioinformatics.
