Usage
=====

.. _installation:

Installation
------------

To use **SPACER**, we recommend setting up a clean Conda environment.

.. code-block:: console

   # Create environment
   $ conda create -n spacer python=3.8
   $ conda activate spacer

   # Install core dependencies
   (spacer) $ pip install torch==2.3.1 numpy==1.26.4 pandas==2.2.2 \
                    scanpy==1.10.2 scikit-learn==1.5.1 scipy==1.13.1 tqdm==4.66.4

Next, clone the SPACER repository and install the package:

.. code-block:: console

   (spacer) $ git clone https://github.com/yaober/SPACER.git
   (spacer) $ cd SPACER
   (spacer) $ pip install -e .

After installation, SPACER can be imported in Python as a standard module.

.. code-block:: python

   import spacer

   model = spacer.load_model("SPACER-v1")
   model.fit(data)

---

Reference
---------

For detailed methodological explanations and biological case studies, 
please refer to the publication section of the documentation.

.. note::

   SPACER is under active development.
   The current release (v1.0) supports mouse and human high-definition spatial transcriptomics datasets.
