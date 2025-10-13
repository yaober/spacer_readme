Installation
============

.. _installation:

Setting up the Environment
--------------------------

To use **SPACER**, it is recommended to create a clean Conda environment.

.. code-block:: console

   # Create environment
   $ conda create -n spacer python=3.8
   $ conda activate spacer

Install the required Python dependencies:

.. code-block:: console

   (spacer) $ pip install torch==2.3.1 numpy==1.26.4 pandas==2.2.2 \
                    scanpy==1.10.2 scikit-learn==1.5.1 scipy==1.13.1 tqdm==4.66.4

---

Installing SPACER
-----------------

Clone the official SPACER repository and install the package in editable mode:

.. code-block:: console

   (spacer) $ git clone https://github.com/yaober/SPACER.git
   (spacer) $ cd SPACER

After downloading , SPACER can be imported:

.. code-block:: python

   from model.dataset import BagsDataset, custom_collate_fn
   from model.model import MIL, EarlyStopping

   dataset = BagsDataset( ... )
   model = MIL( ... )
---


.. note::

   SPACER requires Python ≥3.8 and PyTorch ≥2.3.  
   GPU acceleration is recommended for model training.
