Build Your Own SPACER Model
===========================

This section introduces how to build and train your own **SPACER** model 
using your spatial transcriptomics datasets.

SPACER supports flexible data input, model training, and prediction steps.

---

Data Loading
------------

To start, load your processed spatial transcriptomics data (e.g., Xenium, Visium, Slide-tag)
into an AnnData object compatible with SPACER.

.. code-block:: python

   import spacer
   adata = spacer.load_data("sample_visium_hd.h5ad")

   # Create MIL bags centered on immune cells
   bags = spacer.create_bags(adata, center_type="T_cell", radius=150)

---

Training
--------

SPACER models are trained using a multi-instance learning (MIL) strategy.
Each "bag" represents a focal cell and its local neighborhood.

.. code-block:: python

   from spacer import SPACERModel

   model = SPACERModel()
   model.train(bags, epochs=50, batch_size=32)

   # Save the trained model
   model.save("trained_model.pt")

---

Predicting
----------

Once the model is trained, you can apply it to predict new spatial recruitment scores or
evaluate unseen datasets.

.. code-block:: python

   model = SPACERModel.load("trained_model.pt")

   predictions = model.predict(bags)
   print(predictions.head())

---

.. note::

   SPACER supports GPU acceleration through PyTorch for efficient large-scale training.  
   You can control GPU usage by setting ``device='cuda'`` when initializing the model.
