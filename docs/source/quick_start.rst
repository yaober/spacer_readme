Quick Start
===========

This guide provides a minimal example to run SPACER with your own spatial transcriptomics dataset.

.. code-block:: python

   import spacer
   adata = spacer.load_data("sample_visium_hd.h5ad")
   model = spacer.SPACERModel()
   model.fit(adata)
   model.predict(adata)

For detailed instructions, refer to the “Build Your Own SPACER Model” section.
