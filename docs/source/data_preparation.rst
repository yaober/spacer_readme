Data Preparation
================

Before training, prepare your spatial transcriptomics data in an `.h5ad` format compatible with Scanpy.

1. **Preprocess raw data:**
   - Filter low-quality cells and genes.
   - Normalize and log-transform expression.
   - Annotate cell types if available.

2. **Define neighborhoods:**
   - SPACER constructs multi-instance “bags” centered on target immune or stromal cells.
   - Use spatial radius settings appropriate for your dataset.

Example:
.. code-block:: python

   spacer.create_bags(adata, center_type="T_cell", radius=150)
