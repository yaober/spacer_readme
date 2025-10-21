Quick Start
===========

This guide shows how to train the **spacer** model directly from the command line 
using the provided `train.py` script.  
This approach is ideal for batch experiments or server environments without notebooks.

---

Running spacer from the Command Line
------------------------------------

After installation, navigate to your spacer directory and run:

.. code-block:: console

   (spacer) $ python train.py \
       --data path/to/training_data.h5ad \
       --reference_gene path/to/reference_genes.csv \
       --output_dir results/ \
       --immune_cell tcell \
       --learning_rate 0.0001 \
       --num_epochs 1000 \
       --patience 5 \
       --delta 0.001 \
       --max_instances 500 \
       --n_genes 10000 \
       --selection positive

The script will automatically detect and use a GPU if available (CUDA).

---

Argument Descriptions
---------------------

The following command-line arguments are supported by `train.py`:

+----------------------+---------------------------------------------------------------+
| **Argument**         | **Description**                                               |
+======================+===============================================================+
| `--data`             | Path to the input dataset (e.g., `.h5ad` or `.csv`).          |
+----------------------+---------------------------------------------------------------+
| `--reference_gene`   | Path to a CSV file listing all reference genes.               |
+----------------------+---------------------------------------------------------------+
| `--output_dir`       | Directory where models, metrics, and spacer scores are saved. |
+----------------------+---------------------------------------------------------------+
| `--immune_cell`      | Immune cell type used as bag centers (default: `tcell`).      |
+----------------------+---------------------------------------------------------------+
| `--learning_rate`    | Learning rate for the optimizer (default: 0.0001).            |
+----------------------+---------------------------------------------------------------+
| `--num_epochs`       | Total number of training epochs (default: 1000).              |
+----------------------+---------------------------------------------------------------+
| `--patience`         | Early stopping patience for validation loss (default: 5).     |
+----------------------+---------------------------------------------------------------+
| `--delta`            | Minimum improvement to reset early stopping (default: 0.001). |
+----------------------+---------------------------------------------------------------+
| `--max_instances`    | Maximum number of instances per bag (optional).               |
+----------------------+---------------------------------------------------------------+
| `--n_genes`          | Top N of recruiting cell highly expressed genes to include.   |
+----------------------+---------------------------------------------------------------+
| `--selection`        | Select `"positive (induce)"` or `"negative (repel)"` training.|
+----------------------+---------------------------------------------------------------+

---

Script Overview
---------------

The `train.py` script performs the following steps:

1. **Load Reference Genes**  
   Reads the list of all genes from the specified reference file (`reference_gene`).In this study we use all human/mouse genes as our reference geneset. All human geens are provided in the `data/` folder of the repository:
   - `data/human_reference_genes.csv`

2. **Initialize the Model**  
   Builds a `MIL` (Multi-Instance Learning) model with modules for:
   - distance attention  
   - gene expression weighting  
   - spacer moudle scoring

3. **Create Dataset and DataLoaders**  
   Loads the bag-level dataset via `BagsDataset`, then splits it into 70% training and 30% validation.

4. **Train the Model**  
   Optimizes binary cross-entropy (BCE) loss using the **AdamW** optimizer.  
   Early stopping monitors validation loss (`patience`, `delta`).

5. **Validate and Save Best Model**  
   Evaluates validation AUROC each epoch and saves the best performing weights as `best_model.pth`.

6. **Log Training Metrics**  
   Saves epoch-level metrics (`train_loss`, `val_loss`, `val_AUROC`) to `training_metrics.csv`.

7. **Track spacer Scores**  
   For each epoch, saves `spacer_score_changes_epoch_X.csv`,  
   showing gene-level spacer scores before and after training.

8. **Final Model Output**  
   The fully trained model is stored as `final_model.pth` in your output directory.

---

Example Outputs
---------------

After training completes, your `output_dir` will contain:

.. code-block:: text

   results/
   ├── best_model.pth
   ├── final_model.pth
   ├── training_metrics.csv
   ├── spacer_score_changes_epoch_1.csv
   ├── spacer_score_changes_epoch_2.csv
   └── ...

Each `spacer_score_changes_epoch_X.csv` file summarizes gene-specific immunogenicity
score shifts during training, sorted by magnitude.

---

Tips
----

- **GPU Acceleration**: spacer automatically uses CUDA if available.  
  You can verify this in the log output.
