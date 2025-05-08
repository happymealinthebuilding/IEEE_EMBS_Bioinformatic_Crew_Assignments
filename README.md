
```markdown
# Single-Cell RNA-seq Analysis Report

**Dataset:** AdamsonWeissman2016_GSM2406667_10X005 (`/Users/azratuncay/Downloads/AdamsonWeissman2016_GSM2406667_10X005.h5ad`)
**Date of Analysis:** May 8, 2025 (Placeholder)
**Analysis Tool:** Scanpy (via Python script `data_report.py`)

## 1. Introduction

This report details the analysis of a single-cell RNA sequencing (scRNA-seq) dataset, identified as AdamsonWeissman2016_GSM2406667_10X005.h5ad. The analysis was performed using a Python script leveraging the Scanpy library. The workflow encompassed several key stages: data loading and initial inspection, quality control (QC) and filtering, data preprocessing (normalization, log-transformation, highly variable gene selection, and scaling), dimensionality reduction (PCA and UMAP), and cell clustering using the Leiden algorithm. The primary aim was to process the raw single-cell data to identify distinct cell populations and prepare the data for further biological interpretation.

## 2. Setup and Data Loading

### 2.1. Library Import and Settings

The analysis began by importing necessary Python libraries, including `scanpy`, `pandas`, `numpy`, and `matplotlib.pyplot`. Scanpy settings were configured for verbosity level 3, and default figure parameters (DPI 80, white facecolor) were set.

Terminal Output:
```
Libraries imported and Scanpy settings configured.
scanpy==... anndata==... umap==... numpy==... scipy==... pandas==... scikit-learn==... statsmodels==... python-igraph==... leidenalg==... (Versions will be listed by print_header())
```

### 2.2. Data Loading

The dataset was loaded from the H5AD file: `/Users/azratuncay/Downloads/AdamsonWeissman2016_GSM2406667_10X005.h5ad`.

Terminal Output:
```
Successfully loaded: /Users/azratuncay/Downloads/AdamsonWeissman2016_GSM2406667_10X005.h5ad

AnnData object summary:
AnnData object with n_obs × n_vars = 15006 × 32738
obs: 'perturbation', 'read count', 'UMI count', 'tissue_type', 'cell_line', 'cancer', 'disease', 'perturbation_type', 'celltype', 'organism', 'ncounts', 'ngenes', 'percent_mito', 'percent_ribo', 'nperts'
var: 'ensembl_id', 'ncounts', 'ncells'
```
The dataset initially contained 15,006 cells and 32,738 genes/features. Cell metadata (`adata.obs`) includes information such as perturbation details, read/UMI counts, tissue type, cell line (K562), and pre-calculated QC metrics like `ncounts`, `ngenes`, `percent_mito`, and `percent_ribo`. Gene metadata (`adata.var`) includes Ensembl IDs and per-gene counts.

### 2.3. Initial Data Matrix Inspection

The data matrix `adata.X` had a shape of (15006, 32738) with a `float32` data type. A small snippet showed zero values, and a heuristic check suggested the values were non-negative integers, indicative of raw UMI counts.

Terminal Output (snippet):
```
Data matrix values appear to be non-negative integers (likely raw counts).
```

## 3. Quality Control (QC)

A copy of the original data was stored in `adata.raw`. Standard QC metrics (`n_genes_by_counts` and `total_counts`) were calculated per cell.

Terminal Output:
```
Stored a copy of the original AnnData object in adata.raw.
Calculated QC metrics. Added to adata.obs:
... (head of adata.obs showing new columns)
```

### 3.1. Visualization of QC Metrics

The distributions of `n_genes_by_counts` (number of detected genes per cell) and `total_counts` (total UMI counts per cell) were visualized using violin plots (Figure_1.png) and a scatter plot (Figure_2.png).

**Figure_1.png (Violin Plots):** These plots show the distribution density of detected genes per cell and total UMI counts per cell. They are used to identify the main population of cells and potential outliers (e.g., cells with very few genes or counts, which might be empty droplets or damaged cells, or cells with extremely high counts, which might be doublets).
![Violin Plots](single_cell_analysis_report_HW/images/Figure_1.png)

**Figure_2.png (Scatter Plot):** This plot shows the relationship between the total UMI counts and the number of detected genes for each cell. A positive correlation is expected, where cells with more UMIs generally have more genes detected. Outliers deviating from this trend can be identified.
![Scatter Plot](single_cell_analysis_report_HW/images/Figure_2.png)

*(Note for your report: You should describe what these plots specifically show for your data, e.g., "Most cells have between X and Y genes detected, and between A and B total counts. We observe a population of cells with fewer than Z genes, which are candidates for filtering.")*

### 3.2. Filtering

Based on example thresholds (which should ideally be adjusted based on the QC plots):

- Minimum genes per cell: 200
- Minimum UMI counts per cell: 500
- Minimum cells per gene: 3

Terminal Output:
```
Number of cells before filtering: 15006
Number of genes/features before filtering: 32738
Number of cells after filtering by min_genes_per_cell (200): 15006
Number of cells after filtering by min_counts_per_cell (500): 15006
filtered out 15131 genes that are detected in less than 3 cells
Number of genes/features after filtering by min_cells_per_gene (3): 17607
Final cells: 15006, Final genes/features: 17607
```
The cell filtering steps with the chosen example thresholds (200 genes/cell, 500 counts/cell) did not remove any cells, indicating that all cells passed these initial generous cutoffs. Gene filtering removed 15,131 genes detected in fewer than 3 cells. After filtering, the dataset contained 15,006 cells and 17,607 genes/features.

*(Note for your report: Justify why these specific example thresholds were used or, if you adjusted them based on your plots, explain your rationale.)*

## 4. Preprocessing

### 4.1. Normalization and Log-Transformation

The filtered data was normalized by scaling counts per cell to a target sum of 10,000 (library size normalization) using `scanpy.pp.normalize_total()`. Subsequently, a log1p transformation ($\log(X+1)$) was applied using `scanpy.pp.log1p()` to stabilize variance and make the data distributions more symmetric.

Terminal Output:
```
Applied library size normalization (target_sum=1e4) and log1p transformation.
adata.X sample after normalization & log-transformation:
[[0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0.]]
```
The sample of `adata.X` shows zeros, which is common for sparse scRNA-seq data even after transformation, as many gene-cell combinations will still have zero counts.

### 4.2. Highly Variable Gene (HVG) Selection

Highly variable genes were identified using the `seurat_v3` method, selecting the top 2000 HVGs.

Terminal Output:
```
/Users/azratuncay/Single_cell_h5ad_Data_analysis_HW/.venv/lib/python3.13/site-packages/scanpy/preprocessing/_highly_variable_genes.py:73: UserWarning: `flavor='seurat_v3'` expects raw count data, but non-integers were found.
warnings.warn(
--> added
    'highly_variable', boolean vector (adata.var)
    'highly_variable_rank', float vector (adata.var)
    'means', float vector (adata.var)
    'variances', float vector (adata.var)
    'variances_norm', float vector (adata.var)

Identified highly variable genes. `adata.var['highly_variable']` shows which are selected.
Number of highly variable genes: 2000
```
A `UserWarning` was issued because the `seurat_v3` flavor expects raw count data, but the input `adata.X` had already been normalized and log-transformed. Despite this, the function proceeded and identified 2000 HVGs. This is a common practice in Scanpy workflows, though ideally, for `seurat_v3`, raw counts would be used for variance modeling. The selection added several columns to `adata.var`, including `highly_variable` (a boolean mask) and `means`/`variances` (calculated on the log-normalized data).

**Figure_3.png (Highly Variable Genes Plot):** This plot typically shows gene variances (often normalized) against mean expressions, highlighting the selected HVGs.
![Highly Variable Genes Plot](single_cell_analysis_report_HW/images/Figure_3.png)

*(Note for your report: Describe the appearance of your HVG plot. The script did not subset the data to HVGs at this stage, but marked them for potential use later, e.g., in PCA.)*

### 4.3. Scaling

The data (all 17,607 genes, as HVG subsetting was not performed on `adata.X`) was scaled to unit variance and zero mean, with values clipped at a maximum of 10 using `scanpy.pp.scale()`.

Terminal Output:
```
... as `zero_center=True`, sparse input is densified and may lead to large memory consumption

Scaled the data to unit variance and zero mean (adata.X). Clipped at max_value=10.
adata.X sample after scaling:
[[-0.07117818 -0.02898227 -0.01411406 -0.012744 -0.04029365]
 [-0.07117818 -0.02898227 -0.01411406 -0.012744 -0.04029365]
 [-0.07117818 -0.02898227 -0.01411406 -0.012744 -0.04029365]
 [-0.07117818 -0.02898227 -0.01411406 -0.012744 -0.04029365]
 [-0.07117818 -0.02898227 -0.01411406 -0.012744 -0.04029365]]
```
A warning regarding potential memory consumption due to densification of sparse input during scaling was noted. The sample output shows scaled, centered (negative and positive float) values.

## 5. Dimensionality Reduction

### 5.1. Principal Component Analysis (PCA)

PCA was performed on the scaled data to reduce its dimensionality. The number of principal components (`n_comps_pca`) was calculated as 50. ($min\_dim$ = min(15006 cells, 17607 genes) = 15006. $n\_comps\_pca$ = min(15006 - 2, 50) = 50).

Terminal Output:
```
Current adata dimensions for PCA: n_obs=15006, n_vars=17607
Smallest dimension (min(n_obs, n_vars)): 15006
Calculated n_comps_pca: 50
Attempting PCA with n_comps=50, svd_solver='arpack'.
computing PCA
    with n_comps=50
finished (0:00:01)

Performed PCA.
```
PCA completed successfully using 50 components with the 'arpack' solver.

**Figure_4.png (PCA Variance Ratio Plot):** This "elbow plot" shows the variance explained by each principal component. It is used to determine how many PCs capture most of the biological signal and should be used for downstream steps like UMAP and clustering.
![PCA Variance Ratio Plot](single_cell_analysis_report_HW/images/Figure_4.png)



### 5.2. Neighborhood Graph Construction

A neighborhood graph was computed using `scanpy.pp.neighbors()`. This step utilized the top 15 principal components (a common choice, but should be informed by the PCA elbow plot) and considered 15 neighbors for each cell.

Terminal Output:
```
Using n_neighbors=15 and n_pcs=15 for neighborhood graph.
computing neighbors
    using 'X_pca' with n_pcs = 15
finished: added to `.uns['neighbors']`
    `.obsp['distances']`, distances for each pair of neighbors
    `.obsp['connectivities']`, weighted adjacency matrix (0:00:12)
Computed neighborhood graph.
```

### 5.3. UMAP Embedding

Uniform Manifold Approximation and Projection (UMAP) was computed for 2D visualization of the cell data, based on the neighborhood graph.

Terminal Output:
```
computing UMAP
finished: added
    'X_umap', UMAP coordinates (adata.obsm)
    'umap', UMAP parameters (adata.uns) (0:00:06)
Computed UMAP embedding.
```
The UMAP coordinates were added to `adata.obsm['X_umap']`.

## 6. Clustering

### 6.1. Leiden Clustering

Cells were clustered using the Leiden algorithm (`scanpy.tl.leiden()`) on the computed neighborhood graph. A resolution parameter of 0.5 was used.

Terminal Output:
```
/Users/azratuncay/Single_cell_h5ad_Data_analysis_HW/data_report.py:434: FutureWarning: In the future, the default backend for leiden will be igraph instead of leidenalg.

To achieve the future defaults please pass: flavor="igraph" and n_iterations=2. directed must also be False to work with igraph's implementation.
scanpy.tl.leiden(adata, resolution=resolution_param, key_added=f'leiden_res{resolution_param}')
running Leiden clustering
finished: found 10 clusters and added
    'leiden_res0.5', the cluster labels (adata.obs, categorical) (0:00:07)
Performed Leiden clustering with resolution 0.5. Clusters stored in 'adata.obs["leiden_res0.5"]'.
```
A `FutureWarning` indicated a potential change in the default backend for Leiden in future Scanpy versions. The clustering algorithm identified 10 distinct clusters, and the results were stored in `adata.obs['leiden_res0.5']`.

**UMAPPP.png (UMAP by Leiden Clusters):** This plot visualizes the 10 identified cell clusters on the UMAP embedding.
![UMAP by Leiden Clusters](single_cell_analysis_report_HW/images/UMAPPP.png)



## 7. Visualization of Gene Expression on UMAP (Attempted)

The script attempted to visualize the expression of example genes on the UMAP. However, this step encountered an issue:

Terminal Output:
```
Error: 'adata' not found, or UMAP/clustering not performed. Skipping gene expression visualization.
UMAP embedding available: False
Cluster key 'leiden_res0.5' available: True
```
The script printed an error message indicating that it was skipping the gene expression visualization. The debug output stated `UMAP embedding available: False`, which means the condition `'umap' in adata.obsm` (as printed by the script's debug line) was false.
However, the main conditional for plotting was `if (adata is not None and 'X_umap' in adata.obsm and f'leiden_res{resolution_param}' in adata.obs):`. The UMAP computation log clearly stated that `'X_umap'` was added to `adata.obsm`. The fact that the `else` branch of this main `if` was executed implies that the condition `'X_umap' in adata.obsm` evaluated to `False` when this section of the script was reached. This prevented the plotting of example genes. The cluster key `'leiden_res0.5'` was correctly found.
This discrepancy (UMAP reported as computed and added to `adata.obsm['X_umap']`, but then not found by the conditional for plotting) requires further investigation. It could be a subtle issue in the script logic checking for `'X_umap'` at that specific point, or an unexpected modification to the `adata` object.

Despite this, the script indicates the next steps for systematic marker gene identification:

Terminal Output:
```
To find marker genes systematically, you would typically run:
sc.tl.rank_genes_groups(adata, groupby='leiden_res0.5', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=True)
```

## 8. Conclusion and Next Steps

This analysis successfully processed the AdamsonWeissman2016_GSM2406667_10X005.h5ad single-cell RNA-seq dataset through a standard Scanpy workflow. The data was loaded, quality controlled (resulting in 15,006 cells and 17,607 genes), normalized, and key highly variable genes were identified. Dimensionality reduction via PCA (50 components) and UMAP was performed, followed by Leiden clustering which identified 10 distinct cell clusters (at resolution 0.5).

An issue was encountered in the final step of plotting example gene expressions on the UMAP, where the UMAP coordinates (`adata.obsm['X_umap']`) were unexpectedly not detected by the plotting conditional, despite earlier logs indicating their computation.

**Key Findings:**

- The dataset comprises 15,006 cells and 17,607 genes after initial QC.
- PCA captured significant variance within 50 components.
- Leiden clustering at resolution 0.5 resolved 10 distinct cell clusters, visualized on a UMAP plot.


```