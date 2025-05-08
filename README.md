# README: Single-Cell RNA-seq Analysis Script

This repository contains a Python script for performing a standard single-cell RNA sequencing (scRNA-seq) data analysis workflow using the `scanpy` library. The script covers essential steps from data loading and quality control to dimensionality reduction, clustering, and initial visualization.

## Overview

The script `single_cell_analysis.py` provides a step-by-step analysis pipeline for scRNA-seq data stored in the AnnData (`.h5ad`) format. It is designed to be adaptable, with key parameters highlighted for user modification based on their specific dataset.

## Analysis Steps Covered

1.  **Setup:** Import necessary libraries (`scanpy`, `pandas`, `numpy`, `matplotlib`) and configure Scanpy settings.
2.  **Data Loading:** Load scRNA-seq data from an `.h5ad` file into an AnnData object.
3.  **Initial Data Inspection:** Examine cell metadata (`adata.obs`), gene/feature metadata (`adata.var`), and the expression data matrix (`adata.X`).
4.  **Quality Control (QC):**
    * Preserve raw data.
    * Calculate QC metrics (`total_counts`, `n_genes_by_counts`).
    * Visualize QC metrics using violin plots and scatter plots.
    * Apply filtering based on user-defined thresholds for detected genes and total counts per cell, and minimum cells per gene.
5.  **Preprocessing:**
    * Normalize data to a target sum of 10,000 counts per cell and apply log1p transformation.
    * Identify highly variable genes (HVGs) using the `seurat_v3` flavor. (Note: Requires `scikit-misc`. A warning might appear if data is already normalized/log-transformed).
    * Scale data to unit variance and zero mean.
6.  **Dimensionality Reduction:**
    * Perform Principal Component Analysis (PCA) on the scaled data.
    * Visualize PCA variance ratio ("elbow plot") to help select the number of PCs.
    * Compute a neighborhood graph based on selected PCs.
    * Compute UMAP embedding for 2D visualization.
7.  **Clustering:**
    * Perform Leiden clustering on the neighborhood graph using a specified resolution.
    * Visualize clusters on the UMAP embedding.
8.  **Visualization of Gene Expression (Attempted):**
    * The script attempts to visualize the expression of a few example genes on the UMAP plot alongside the clusters. (Note: An issue was noted in the analysis report where the UMAP coordinates were not detected for plotting despite being computed).

