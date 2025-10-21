Data Preparation
================

Before training **spacer**, prepare your spatial transcriptomics dataset in a `.h5ad` format compatible with **Scanpy**.  
This section describes the preprocessing, annotation, and structure requirements for input data.

---

Preprocessing Workflow
----------------------

spacer expects a **Scanpy AnnData** object (`adata`) as input.  
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

spacer requires the `.h5ad` object to contain the following essential fields:

- **Expression matrix (`adata.X`)**: normalized and log-transformed gene expression values.  
- **Metadata table (`adata.obs`)**: must contain the following columns:

  +----------------+--------------------------------------------------------------+
  | **Column**     | **Description**                                              |
  +================+==============================================================+
  | `X`, `Y`       | Spatial coordinates (in microns or pixel units) of each cell |
  +----------------+--------------------------------------------------------------+
  | `cell_type`    | Integer-encoded major cell class                             |
  |                | - `0` → other cell                                           |
  |                | - `1` → recruiting cell                                      |
  |                | - `2` → engaging cell (customizable, e.g., T/B/macrophage)   |
  +----------------+--------------------------------------------------------------+
  | `<EngagingTag>`| Binary indicator for whether the cell at the center of the   |
  |                | bag belongs to the target engaging cell type (`1`) or not (`0`) |
  +----------------+--------------------------------------------------------------+


The `<EngagingTag>` column defines which cells will serve as the **center** for each neighborhood (“bag”) in spacer.  
Below is the default mapping used in our study:

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

---

Customizing the Mapping
-----------------------

In this work, we used the above mapping to ensure **consistent annotation across datasets** involving multiple stromal and immune cell types.  
Each key in the mapping corresponds to a general immune or stromal population, while the assigned value (e.g., `"T"`, `"B"`, `"Macrophage"`)  
serves as a compact label for downstream modeling and visualization.  

However, this mapping is **fully customizable**.  
Users can freely modify or extend it to match their experimental context or cell annotation schema.  
For instance, if you are analyzing **brain tissues**, you could define:

.. code-block:: python

   mapping = {
       'microglia': 'Microglia',
       'astrocyte': 'Astrocyte',
       'oligodendrocyte': 'Oligodendrocyte',
   }

spacer will automatically adapt its bag construction and learning process to your new mapping.  
The only requirement is that the corresponding binary column (e.g., `"Microglia"`) exists in `adata.obs`  
with values `1` for target cells and `0` otherwise.

