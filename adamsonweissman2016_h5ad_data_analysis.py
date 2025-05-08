# # Single-Cell RNA-seq Analysis 
#
# This script outlines the steps for analyzing single-cell RNA sequencing data using Scanpy.
# It will cover data loading, quality control, preprocessing, dimensionality reduction,
# clustering, and visualization.

# I. Setup: Importing Libraries
import scanpy 
import pandas as pd
import numpy
import matplotlib.pyplot as plt

# Scanpy settings
scanpy.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
scanpy.logging.print_header()
scanpy.settings.set_figure_params(dpi=80, facecolor='white')

print("Libraries imported and Scanpy settings configured.")

# ## 1. Setup: Importing Libraries
#
# The first step in our analysis is to import the necessary Python libraries.
# * `scanpy` (as `sc`): This is the core library for single-cell analysis.
# * `pandas` (as `pd`): Used for data manipulation, especially for annotations.
# * `numpy` (as `np`): Fundamental package for numerical computation.
# * `matplotlib.pyplot` (as `plt`): A widely used plotting library.
#
# We also configure some basic Scanpy settings for verbosity and figure appearance.

# Initialize adata to None. This is important for subsequent checks.
adata = None
resolution_param = 0.5 # Define resolution_param globally or ensure it's defined before use in cell 15

# II. Data Loading
# Define the file path to your data - UPDATED PATH
file_path = "/Users/azratuncay/Downloads/AdamsonWeissman2016_GSM2406677_10X005.h5ad"

# Load the data into an AnnData object
try:
    adata = scanpy.read_h5ad(file_path)
    print(f"Successfully loaded: {file_path}")
except FileNotFoundError:
    print(f"Error: The file was not found at {file_path}")
    print("Please double-check the file path and ensure the file exists.")
except Exception as e:
    print(f"An error occurred while loading the file: {e}")

# Display the AnnData object to see a summary
# and exit if data loading failed
if adata is not None:
    print("\nAnnData object summary:")
    print(adata)
else:
    print("\nData loading failed. Exiting script.")
    exit() # Exit the script if adata could not be loaded

# ## 2. Data Loading
#
# We load our single-cell dataset from an `h5ad` file. This format efficiently stores
# the gene expression matrix, cell annotations (`obs`), gene annotations (`var`),
# and other data layers.
#
# The `AnnData` object summary shows its dimensions (number of cells × number of genes/features)
# and available annotations.

# III. Initial Data Inspection: Cell Metadata (adata.obs)
if adata is not None:
    print("\nFirst 5 rows of Cell Metadata (adata.obs):")
    print(adata.obs.head()) # Changed display() to print()

    print(f"\nTotal number of cells: {adata.n_obs}")
    print(f"Available cell annotations (columns in adata.obs): {list(adata.obs.columns)}")
else:
    # This part should not be reached if the script exits on loading failure
    print("Error: 'adata' not loaded.")

# ## 3. Initial Data Inspection: Cell Metadata (`adata.obs`)
#
# We inspect the cell metadata (`adata.obs`), which contains annotations for each cell
# (e.g., batch, cell type, QC metrics).

# IV. Initial Data Inspection: Gene/Feature Metadata (adata.var)
if adata is not None:
    print("\nFirst 5 rows of Gene/Feature Metadata (adata.var):")
    print(adata.var.head()) # Changed display() to print()

    print(f"\nTotal number of genes/features: {adata.n_vars}")
    print(f"Available gene/feature annotations (columns in adata.var): {list(adata.var.columns)}")
else:
    print("Error: 'adata' not loaded.")

# ## 4. Initial Data Inspection: Gene/Feature Metadata (`adata.var`)
#
# Next, we inspect the gene/feature metadata (`adata.var`), containing annotations
# for each gene or feature (e.g., gene IDs, feature types).

# V. Initial Data Inspection: Data Matrix (adata.X)
if adata is not None:
    print("\nData matrix (adata.X) information:")
    print(f"Shape of the data matrix (cells x features): {adata.X.shape}")
    print(f"Data type of the matrix: {adata.X.dtype}")

    print("\nSmall snippet of the data matrix (adata.X[:5, :5]):")
    if hasattr(adata.X, 'toarray'):
        try:
            print(adata.X[:5, :5].toarray())
        except AttributeError:
            print(adata.X[:5, :5]) # For sparse matrices that might not have .toarray() directly on slice
    else:
        print(adata.X[:5, :5]) # For dense numpy arrays

    # Heuristic check for data type (counts vs normalized)
    if adata.X.shape[0] > 0 and adata.X.shape[1] > 0: # Ensure matrix is not empty
        sample_dense = adata.X[:min(100, adata.X.shape[0]), :min(100, adata.X.shape[1])]
        if hasattr(sample_dense, 'toarray'):
            sample_dense = sample_dense.toarray()
        sample_values = sample_dense.flatten()

        if numpy.all(sample_values >= 0):
            is_float_present = numpy.any(numpy.abs(sample_values - numpy.round(sample_values)) > 1e-9) # Check for non-integer values
            if is_float_present or numpy.max(sample_values) > 1000: # Heuristic: floats or large integers
                print("\nData matrix values appear to be non-negative (possibly raw counts, TPMs/CPMs or already normalized).")
            elif numpy.all(numpy.equal(numpy.mod(sample_values, 1), 0)): # Check if all are integers
                 print("\nData matrix values appear to be non-negative integers (likely raw counts).")
        else:
            print("\nData matrix values include negatives (possibly log-transformed and centered/scaled data).")
    else:
        print("\nData matrix is empty or too small for heuristic check.")
else:
    print("Error: 'adata' not loaded.")

# ## 5. Initial Data Inspection: Data Matrix (`adata.X`)
#
# The primary data is stored in `adata.X` (cells × genes/features). We check its shape,
# data type, and a small snippet. This helps understand if it's raw counts, normalized,
# or log-transformed data, guiding preprocessing choices.

# VI. Quality Control (QC) - Step 1: Preserve Raw Data
if adata is not None:
    adata.raw = adata.copy() # Store a copy of the raw data before filtering/normalization
    print("\nStored a copy of the original AnnData object in adata.raw.")
else:
    print("Error: 'adata' not loaded.")

# # Quality Control (QC)
#
# Before analysis, we perform QC to remove low-quality cells or uninformative genes/features.
# For scRNA-seq, common metrics include the number of genes detected per cell,
# total counts per cell, and often mitochondrial gene percentage.
# *(Note: This workflow currently doesn't explicitly calculate mitochondrial/ribosomal
# gene percentages, which often requires specific gene annotations in `adata.var`.
# If available, these would be valuable additional QC metrics.)*

# ## 1. Preserve Raw Data (already done as VI)
# It's good practice to store a copy of the original data before filtering.

# QC - Step 2: Calculate QC Metrics (Scanpy calls this section VII in the notebook)
if adata is not None:
    scanpy.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    # For RNA data, 'n_genes_by_counts' is the number of detected genes.
    # 'total_counts' is the total UMI counts for that cell.
    print("\nCalculated QC metrics. Added to adata.obs:")
    print(adata.obs.head()) # Changed display() to print()
else:
    print("Error: 'adata' not loaded.")

# ## 2. Calculate QC Metrics
#
# We use `sc.pp.calculate_qc_metrics` to compute:
# * `total_counts`: Total UMI counts per cell.
# * `n_genes_by_counts`: Number of unique genes detected per cell.
# These are added to `adata.obs`.

# QC - Step 3: Visualize QC Metrics
if adata is not None:
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    scanpy.pl.violin(adata, 'n_genes_by_counts', jitter=0.4, ax=axs[0], show=False)
    axs[0].set_title('Number of Detected Genes')
    axs[0].set_ylabel('Genes per cell')

    scanpy.pl.violin(adata, 'total_counts', jitter=0.4, ax=axs[1], show=False)
    axs[1].set_title('Total UMI Counts')
    axs[1].set_ylabel('Counts per cell')
    plt.tight_layout()
    plt.show() # Shows the violin plots

    # New figure for the scatter plot
    plt.figure(figsize=(6, 4))
    plt.scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], alpha=0.5, s=5)
    plt.xlabel('Total UMI Counts per Cell')
    plt.ylabel('Number of Detected Genes per Cell')
    plt.title('Counts vs. Detected Genes')
    # Consider plt.xscale('log') and plt.yscale('log') if ranges are very large
    plt.show() # Shows the scatter plot
else:
    print("Error: 'adata' not loaded.")

# ## 3. Visualize QC Metrics
#
# Violin plots and scatter plots help visualize the distribution of these QC metrics,
# aiding in setting filtering thresholds to remove outliers (e.g., cells with very few
# genes/counts, potentially empty droplets or damaged cells).
# **(NOTE: Interpret your plots here. Describe the distributions and identify potential cutoffs.)**

# QC - Step 4: Apply Filters (ADJUST THRESHOLDS BASED ON YOUR DATA!)
if adata is not None:
    print(f"\nNumber of cells before filtering: {adata.n_obs}")
    print(f"Number of genes/features before filtering: {adata.n_vars}")

    # --- ADJUST THESE THRESHOLDS CAREFULLY BASED ON YOUR PLOTS ---
    min_genes_per_cell = 200     # Example: minimum number of genes detected per cell
    min_counts_per_cell = 500    # Example: minimum total UMI counts per cell
    min_cells_per_gene = 3       # Example: gene must be detected in at least this many cells
    # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

    # Apply filters
    scanpy.pp.filter_cells(adata, min_genes=min_genes_per_cell)
    print(f"Number of cells after filtering by min_genes_per_cell ({min_genes_per_cell}): {adata.n_obs}")

    scanpy.pp.filter_cells(adata, min_counts=min_counts_per_cell)
    print(f"Number of cells after filtering by min_counts_per_cell ({min_counts_per_cell}): {adata.n_obs}")

    scanpy.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    print(f"Number of genes/features after filtering by min_cells_per_gene ({min_cells_per_gene}): {adata.n_vars}")
    print(f"Final cells: {adata.n_obs}, Final genes/features: {adata.n_vars}")
else:
    print("Error: 'adata' not loaded.")

# ## 4. Apply Filters
#
# Based on the QC visualizations, we filter out low-quality cells and genes/features
# not detected in enough cells. **Thresholds are data-dependent and must be adjusted.**
# **(NOTE: Justify your chosen thresholds based on your data's distributions.)**

# VII. Preprocessing - Step 1: Normalization (This matches notebook cell numbering)
if adata is not None:
    # Normalize to a target sum of 10,000 counts per cell (library size normalization)
    scanpy.pp.normalize_total(adata, target_sum=1e4)

    # Log-transform the data
    scanpy.pp.log1p(adata)
    print("\nApplied library size normalization (target_sum=1e4) and log1p transformation.")
    print("adata.X sample after normalization & log-transformation:")
    if hasattr(adata.X, 'toarray'):
        print(adata.X[:5, :5].toarray())
    else:
        print(adata.X[:5, :5])
else:
    print("Error: 'adata' not loaded.")

# # Preprocessing
# After QC, we preprocess the data.

# ## 1. Normalization and Log-Transformation
# * **Library Size Normalization (`sc.pp.normalize_total`):** Scales the counts in each cell
#   to a common target sum (e.g., 10,000), correcting for differences in sequencing depth.
# * **Log Transformation (`sc.pp.log1p`):** Applies log(X+1) to the normalized counts.
#   This helps stabilize variance and make the data more symmetric, as scRNA-seq data
#   is often highly skewed.

# Preprocessing - Step 2: Identify Highly Variable Genes (HVGs)
if adata is not None:
    # Before running this, ensure 'scikit-misc' is installed if you use flavor='seurat_v3'
    # You can install it via: pip install scikit-misc (in your venv terminal)
    try:
        scanpy.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', subset=False)
        # subset=False adds 'highly_variable' column to adata.var
        # subset=True would subset adata to only HVGs
        print("\nIdentified highly variable genes. `adata.var['highly_variable']` shows which are selected.")
        print(f"Number of highly variable genes: {adata.var['highly_variable'].sum()}")
        scanpy.pl.highly_variable_genes(adata, show=True) # Shows HVG plot

        # Optional: If you want to proceed with only HVGs for PCA and downstream:
        # adata = adata[:, adata.var.highly_variable].copy()
        # print(f"Subsetting to {adata.n_vars} highly variable genes.")

    except ImportError as e:
        print(f"\nImportError for sc.pp.highly_variable_genes (flavor='seurat_v3'): {e}")
        print("Please ensure 'scikit-misc' is installed: run 'pip install scikit-misc' in your terminal with .venv activated, then re-run this script.")
    except Exception as e:
        print(f"\nAn error occurred during HVG selection: {e}")
else:
    print("Error: 'adata' not loaded.")

# ## 2. Identify Highly Variable Genes (HVGs)
#
# HVG selection focuses on genes that show the most biological variation across cells,
# reducing noise and computational load for downstream analyses. We use the `seurat_v3` method.
# *(Note: The `seurat_v3` method requires the `scikit-misc` package.
# If not installed, run `pip install scikit-misc` in your venv terminal and restart VS Code kernel/re-run script.)*
#
# For now, we mark HVGs in `adata.var['highly_variable']` but don't subset `adata` yet.
# You can choose to subset later if desired.
# **(NOTE: Comment on the HVG plot and your decision whether to subset or not.)**

# Preprocessing - Step 3: Scaling the Data
if adata is not None:
    # Scale data to unit variance and zero mean.
    # Apply only to HVGs if you subsetted adata to them in the previous step.
    # Otherwise, applies to all genes if no subsetting was done.
    scanpy.pp.scale(adata, max_value=10)
    print("\nScaled the data to unit variance and zero mean (adata.X). Clipped at max_value=10.")
    print("adata.X sample after scaling:")
    if hasattr(adata.X, 'toarray'):
        print(adata.X[:5, :5].toarray())
    else:
        print(adata.X[:5, :5])
else:
    print("Error: 'adata' not loaded.")

# ## 3. Scaling the Data
#
# Scaling (`sc.pp.scale`) shifts the expression of each gene to have a mean of zero
# and unit variance across cells. This is important before PCA to ensure genes with
# higher expression levels don't disproportionately influence the analysis.
# Values are often clipped (e.g., at `max_value=10`) to prevent extreme outliers
# from dominating. If you subsetted to HVGs, this step applies only to them.
# Otherwise, it applies to all genes.

# VIII. Dimensionality Reduction - Step 1: Principal Component Analysis (PCA)
if adata is not None:
    print(f"\nCurrent adata dimensions for PCA: n_obs={adata.n_obs}, n_vars={adata.n_vars}")

    # Determine the smallest dimension of the data
    min_dim = min(adata.n_obs, adata.n_vars)
    print(f"Smallest dimension (min(n_obs, n_vars)): {min_dim}")

    n_comps_pca = 0 # Initialize
    if min_dim <= 1: # If 0 or 1, PCA is not possible/meaningful in the usual sense
        n_comps_pca = 0
    elif min_dim == 2: # Can technically do 1 PC
        n_comps_pca = 1
    else: # min_dim > 2
        # For arpack, n_comps must be < rank. Max rank is min_dim.
        # If features became constant after scaling, effective rank might be min_dim - 1 (or less).
        # So, aim for n_comps <= (effective rank) - 1, which is roughly min_dim - 2.
        # Also cap at 50 components as a general heuristic.
        n_comps_pca = min(min_dim - 2, 50)
        if n_comps_pca < 1: # Ensure at least 1 if calculation results in 0 or less but min_dim > 2
             n_comps_pca = 1


    print(f"Calculated n_comps_pca: {n_comps_pca}")

    if n_comps_pca < 2: # UMAP/clustering usually needs at least 2 PCs, often more.
        print(f"Warning: n_comps_pca is {n_comps_pca}. Meaningful PCA for downstream UMAP/clustering typically requires more components.")
        print("Consider revisiting QC/filtering if too few cells/features remain, or if many features have zero variance.")
        # Depending on the workflow, you might want to exit() or skip subsequent steps.
        # For now, we'll allow it to proceed but PCA results might be trivial.
        if n_comps_pca == 0:
            print("PCA cannot be run with 0 components. Skipping PCA and dependent steps.")
        elif n_comps_pca > 0 : # If n_comps_pca is 1
             try:
                print(f"Attempting PCA with n_comps={n_comps_pca}, svd_solver='arpack'.")
                scanpy.tl.pca(adata, svd_solver='arpack', n_comps=n_comps_pca) # use_highly_variable=True if you want to use only HVGs without subsetting adata
                print("\nPerformed PCA.")
                n_pcs_for_plot = min(n_comps_pca, adata.obsm['X_pca'].shape[1])
                if n_pcs_for_plot > 0:
                    scanpy.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs_for_plot, show=True) # Shows PCA variance plot
                else:
                    print("No PCs available to plot variance ratio.")
             except ValueError as e:
                print(f"\nPCA with 'arpack' failed: {e}")
                print("Consider reducing n_comps_pca or using svd_solver='auto'.")

    else: # n_comps_pca >= 2
        print(f"Attempting PCA with n_comps={n_comps_pca}, svd_solver='arpack'.")
        try:
            # If you want to run PCA only on HVGs (and 'adata' was NOT subsetted to HVGs):
            # sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps_pca, use_highly_variable=True)
            # Otherwise, if 'adata' IS ALREADY SUBSETTED to HVGs, or you want to use all (scaled) genes:
            scanpy.tl.pca(adata, svd_solver='arpack', n_comps=n_comps_pca)
            print("\nPerformed PCA.")

            # Ensure n_pcs for plotting variance ratio doesn't exceed available components
            n_pcs_for_plot = min(n_comps_pca, adata.obsm['X_pca'].shape[1])
            if n_pcs_for_plot > 0 :
                 scanpy.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs_for_plot, show=True) # Shows PCA variance plot
            else:
                 print("No PCs available to plot variance ratio.")
        except ValueError as e:
            print(f"\nPCA with 'arpack' failed: {e}")
            print("This can happen if n_comps is not strictly less than the matrix rank (number of non-constant features).")
            print(f"You have n_vars={adata.n_vars}. If one or more features became constant after scaling, the effective rank is reduced.")
            print(f"Tried with n_comps_pca={n_comps_pca}. Ensure this is less than the effective rank.")
            print("Consider reducing n_comps_pca further or using svd_solver='auto'.")
else:
    print("Error: 'adata' not loaded.")

# # Dimensionality Reduction
#
# We reduce the dimensionality of the data to capture the main axes of variation.
# ## 1. Principal Component Analysis (PCA)
#
# PCA performs linear dimensionality reduction. We run it on the scaled data
# (potentially restricted to HVGs). The "variance ratio" plot (elbow plot) helps
# decide how many PCs to retain for downstream analysis, looking for an "elbow"
# where adding more PCs gives diminishing returns.
# *(Note on `n_comps_pca`: The number of components for PCA with the 'arpack' solver
# must be strictly less than the effective rank of the data matrix. If features
# become constant (e.g., all zeros) after scaling, this reduces the effective rank.
# The code attempts to set `n_comps_pca` conservatively.)*
# **(NOTE: Interpret your PCA variance plot: "The plot shows that the first X PCs
# explain a significant portion of the variance. We will use Y PCs for
# neighborhood graph construction...")**

# Dimensionality Reduction - Step 2: Neighborhood Graph, UMAP, and Clustering
if adata is not None and 'X_pca' in adata.obsm and adata.obsm['X_pca'].shape[1] > 0:
    # --- CHOOSE THE NUMBER OF PCs BASED ON YOUR PCA ELBOW PLOT ---
    n_neighbors_param = 15 # Default, adjust based on dataset size and desired resolution
    n_pcs_param = 15       # Example: choose based on your PCA variance plot
                           # (e.g., if n_comps_pca was 22, you might use 15-20 PCs)
    # --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

    # Ensure n_pcs_param is not more than available PCs
    available_pcs = adata.obsm['X_pca'].shape[1]
    if n_pcs_param > available_pcs:
        print(f"Warning: n_pcs_param ({n_pcs_param}) is greater than available PCs ({available_pcs}). Adjusting to max available.")
        n_pcs_param = available_pcs

    if n_pcs_param < 2 : # Need at least some PCs for neighborhood calculation
        print(f"Error: n_pcs_param is {n_pcs_param}, which is too low for neighborhood calculation. Ensure PCA ran successfully with enough components.")
    else:
        print(f"\nUsing n_neighbors={n_neighbors_param} and n_pcs={n_pcs_param} for neighborhood graph.")
        scanpy.pp.neighbors(adata, n_neighbors=n_neighbors_param, n_pcs=n_pcs_param)
        print("Computed neighborhood graph.")

        # UMAP for visualization
        scanpy.tl.umap(adata)
        print("Computed UMAP embedding.")

        # Clustering (e.g., Leiden algorithm)
        # resolution_param was defined globally earlier
        scanpy.tl.leiden(adata, resolution=resolution_param, key_added=f'leiden_res{resolution_param}')
        print(f"Performed Leiden clustering with resolution {resolution_param}. Clusters stored in 'adata.obs[\"leiden_res{resolution_param}\"]'.")

        # Visualize UMAP with clusters
        plt.figure() # Ensures a new figure for the UMAP plot
        scanpy.pl.umap(adata, color=[f'leiden_res{resolution_param}'], legend_loc='on data', title=f'Leiden Clusters (res {resolution_param})', show=True) # Shows UMAP plot
else:
    print("\nError: 'adata' object not found or PCA ('X_pca') not available or has 0 components. Skipping Neighborhood Graph, UMAP, and Clustering.")

# ## 2. Neighborhood Graph, UMAP Embedding, and Clustering

#
# Using the selected PCs:
# 1.  **Neighborhood Graph (`sc.pp.neighbors`):** Constructs a graph where cells are
#     connected if they are close in PCA space.
# 2.  **UMAP Embedding (`sc.tl.umap`):** Non-linear dimensionality reduction for 2D visualization.
# 3.  **Clustering (`sc.tl.leiden`):** Partitions cells into clusters based on the
#     neighborhood graph. The `resolution` parameter affects cluster granularity.
# **(NOTE: Interpret your UMAP: "The UMAP shows distinct cell clusters. We used X PCs
# and a Leiden resolution of Y, resulting in Z clusters...")**

# IX. Visualizing Gene Expression on UMAP (Example)
if (adata is not None and
    'X_umap' in adata.obsm and # Corrected check for UMAP coordinates
    f'leiden_res{resolution_param}' in adata.obs):
    
    example_genes_to_plot = []
    if adata.n_vars > 0:
        num_genes_to_sample = min(4, adata.n_vars) # Plot up to 4 example genes
        # Ensure var_names are strings, and handle potential non-unique indices if any
        if not pd.api.types.is_string_dtype(adata.var_names):
            adata.var_names = adata.var_names.astype(str)
        if not adata.var_names.is_unique:
            print("Warning: Gene names (adata.var_names) are not unique. Making them unique for plotting.")
            adata.var_names_make_unique()

        example_genes_to_plot = list(adata.var_names[:num_genes_to_sample])
    else:
        print("No genes available in adata.var to plot.")

    if example_genes_to_plot:
        print(f"\nPlotting example genes: {example_genes_to_plot}")

        # Add the cluster key to the list of things to color by
        color_by = example_genes_to_plot + [f'leiden_res{resolution_param}']

        # Create titles for each plot
        plot_titles = [f'{gene} Expression' for gene in example_genes_to_plot] + [f'Leiden Clusters (res {resolution_param})']

        # Ensure titles list matches color_by list length if Scanpy matches them by order
        # Scanpy's title argument can take a single string or a list matching 'color'
        final_titles = plot_titles if len(plot_titles) == len(color_by) else None


        scanpy.pl.umap(adata,
                    color=color_by,
                    ncols=min(3, len(color_by)), # Adjust ncols for better layout
                    title=final_titles, # Use the generated titles
                    show=True) # Shows UMAP plots with gene expression
    elif adata.n_vars > 0: # Genes exist, but example_genes_to_plot is empty
        print("\nCould not select example genes to plot, though adata.n_vars > 0.")
    # If adata.n_vars is 0, the earlier message "No genes available..." would have printed.

    # For a more systematic exploration of marker genes per cluster:
    print("\nTo find marker genes systematically, you would typically run:")
    print(f"sc.tl.rank_genes_groups(adata, groupby='leiden_res{resolution_param}', method='wilcoxon')")
    print(f"sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=True)")
else:
    print("\nError: 'adata' not found, or UMAP/clustering not performed. Skipping gene expression visualization.")
    if adata is not None: # More detailed error context
        print(f"UMAP embedding available: {'umap' in adata.obsm}")
        # Check if resolution_param was defined and if the key exists
        cluster_key = f'leiden_res{resolution_param}' if 'resolution_param' in globals() else 'leiden_res[undefined_resolution]'
        print(f"Cluster key '{cluster_key}' available: {cluster_key in adata.obs if 'resolution_param' in globals() else False}")

# ## 3. Visualizing Gene Expression on UMAP
#
# To understand cluster identities, we visualize the expression of individual genes
# on the UMAP. This example plots the first few genes from your dataset.
# More systematically, `sc.tl.rank_genes_groups` is used to find marker genes
# that differentiate clusters.
# **(NOTE: Select known marker genes or interesting genes from `adata.var_names` to plot.
# Interpret the plots: "Gene X shows high expression primarily in cluster A, suggesting...")**

print("\n--- End of Single-Cell Analysis Script ---")