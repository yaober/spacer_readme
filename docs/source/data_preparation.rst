Data Preparation
================

Before training SPACER, prepare your spatial transcriptomics dataset in a `.h5ad` format compatible with **Scanpy**.  
This section describes the preprocessing, annotation, and structure requirements for input data.

---

Preprocessing Workflow
----------------------

SPACER expects a **Scanpy AnnData** object (`adata`) as input.  
You can generate this file from raw count matrices using **Scanpy** following these preprocessing steps:

.. code-block:: python

   import scanpy as sc

   # Load raw data (example)
   adata = sc.read_10x_h5("sample_filtered_feature_bc_matrix.h5")

   # Step 1: Filter low-quality data
   sc.pp.filter_cells(adata, min_genes=n)
   sc.pp.filter_genes(adata, min_cells=m)

   # Step 2: Normalize and log-transform
   sc.pp.normalize_total(adata, target_sum=1e4)
   sc.pp.log1p(adata)

   # Step 3: Annotate cell types (user-defined)
   # e.g., adata.obs['cell_type'] = cell_type_annotation_vector

---

Required Data Structure
-----------------------

SPACER requires the `.h5ad` object to contain the following essential fields:

- **Expression matrix (`adata.X`)**: normalized and log-transformed gene expression values.
- **Metadata table (`adata.obs`)**: must contain the following columns:

  +----------------+--------------------------------------------------------------+
  | **Column**     | **Description**                                              |
  +================+==============================================================+
  | `X`, `Y`       | Spatial coordinates (in microns or pixel units) of each cell |
  +----------------+--------------------------------------------------------------+
  | `cell_type`    | Integer-encoded major cell class                             |
  |                |  - `0` → non-tumor, non-immune cell                          |
  |                |  - `1` → tumor cell                                          |
  |                |  - `2` → immune cell (customizable, e.g., T/B/macrophage)    |
  +----------------+--------------------------------------------------------------+
  | `<ImmuneTag>`  | Binary indicator for target immune cell type (1 = target)    |
  +----------------+--------------------------------------------------------------+

The `<ImmuneTag>` column is determined by your study focus.  
SPACER uses it to define the **center cell type** for neighborhood (“bag”) construction.  
Below is the default mapping used in SPACER:

.. code-block:: python

   mapping = {
       'tcell': 'T',
       'bcell': 'B',
       'macrophage': 'Macrophage',
       'neutrophil': 'Neutrophil',
       'fibroblast': 'Fibroblast',
       'endothelial': 'Endothelial',
   }

For example, if you are studying **T-cell recruitment**, the column name in `adata.obs` should be `"T"`,  
and its values should be `1` for T cells and `0` for all other cells.

