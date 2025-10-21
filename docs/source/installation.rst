Installation
============

.. _installation:

Setting up the Environment
--------------------------

To use **spacer**, it is recommended to create a clean Conda environment.

.. code-block:: console

   # Create environment
   $ conda create -n spacer python=3.8
   $ conda activate spacer

Install the required Python dependencies:

.. code-block:: console

   (spacer) $ pip install torch==2.3.1 numpy==1.26.4 pandas==2.2.2 \
                    scanpy==1.10.2 scikit-learn==1.5.1 scipy==1.13.1 tqdm==4.66.4

---

Installing Spacer
-----------------

Clone the official spacer repository and install the package in editable mode:

.. code-block:: console

   (spacer) $ git clone https://github.com/yaober/SPACER.git
   (spacer) $ cd SPACER

After downloading, spacer can be imported:

.. code-block:: python

   from model.dataset import BagsDataset, custom_collate_fn
   from model.model import MIL, EarlyStopping

   dataset = BagsDataset( ... )
   model = MIL( ... )

---

.. note::

   spacer requires Python ≥3.8 and PyTorch ≥2.3.  
   GPU acceleration is recommended for model training. However, GPU usage is **not mandatory** — the model can also be trained on CPUs.  
   We recommend having at least **256 GB of system memory (RAM)** if running on CPU, to ensure smooth data loading and model optimization during training.
