Build Your Own SPACER Model
===========================

This section introduces how to construct and apply your own **SPACER** model 
using spatial transcriptomics datasets.

SPACER supports flexible input formats, automated neighborhood construction, and 
multi-instance learning (MIL)-based prediction.

---

Neighborhood Construction
-------------------------

SPACER organizes spatial transcriptomics data into **multi-instance “bags”**,  
where each bag is centered on a target (e.g., immune) cell and includes its spatial neighbors within a defined radius.

You can create these bags using the ``spacer.create_bags()`` function:

.. code-block:: python

   from spacer import create_bags

   bags = create_bags(
       adata,
       center_type="T_cell",   # or any immune/stromal cell type
       radius=150              # spatial radius in microns
   )

   from spacer.data import BagsDataset

   dataset = BagsDataset(
       input_data="sample_visium_hd.h5ad",
       immune_cell="tcell",
       max_instances=500,
       n_genes=3000,
       k=2   # Ensure 'k' matches the number of bags per batch
   )

Recommended radius values depend on dataset resolution:

- **Visium (standard):** 150 μm  
- **Visium HD:** 75 μm  
- **Slide-tag or low-density tissue:** 150 μm  
- **Cardiac tissue:** up to 150 μm (to account for large cardiomyocytes)

SPACER will automatically exclude cells without sufficient neighbors  
(default minimum = 10 instances per bag).

---

Output
------

The output of ``create_bags()`` is a processed dataset ready for model training or prediction,  
containing:

- **Center cell indices** (target immune or stromal cells)  
- **Neighbor indices** within the defined radius  
- **Expression matrices** for each instance within a bag  

Each “bag” represents a **microenvironmental unit** for downstream multi-instance learning.

---

Prediction
----------

Once a model has been trained, you can load it and perform prediction on your processed dataset.

.. code-block:: python

   from spacer import SPACERModel

   # Load a pretrained SPACER model
   model = SPACERModel.load("trained_model.pt")

   # Predict recruitment or engagement scores
   predictions = model.predict(dataset)
   print(predictions.head())

For details on how to train a model, please refer to the  
:doc:`Quick Start <quick_start>` section.

---

Jupyter Notebook Example
------------------------

You can explore the entire workflow interactively in a Jupyter Notebook:

.. code-block:: python

   import spacer
   from spacer import SPACERModel

   # 1. Load data
   adata = spacer.load_data("sample_visium_hd.h5ad")

   # 2. Create neighborhood bags
   bags = spacer.create_bags(adata, center_type="T_cell", radius=150)

   # 3. Load pretrained model and predict
   model = SPACERModel.load("trained_model.pt")
   predictions = model.predict(bags)

   print(predictions.head())

.. note::

   - SPACER supports GPU acceleration via PyTorch for efficient inference.  
   - You can specify the device using ``device='cuda'`` during initialization.  
   - A runnable version of this workflow is available at  
     ``docs/notebooks/SPACER_QuickStart.ipynb``.

---

Next Steps
----------

After prediction, you can:

- Visualize predicted recruitment scores spatially on tissue sections  
- Correlate predictions with immunopeptidomics or TCR data  
- Perform gene-level importance analysis via SPACER’s interpretability module
