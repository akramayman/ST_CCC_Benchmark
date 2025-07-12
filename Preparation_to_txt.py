#!/usr/bin/env python
# coding: utf-8

# # CSV Gene Expression and Metadata Processor with CPM Normalization

# ### Import packages 

# In[10]:


import os
import pandas as pd
import numpy as np


# # This function reads gene expression and metadata CSV files

# In[ ]:


def process_sample(gene_expr_path, metadata_path, output_dir, sample_name="sample"):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load the csv files 
    expr = pd.read_csv(gene_expr_path, index_col=0)
    meta_raw = pd.read_csv(metadata_path, index_col=0)
    expr.index = expr.index.astype(str)
    meta_raw.index = meta_raw.index.astype(str)

    # Normalize gene expression to CPM
    cpm = expr.div(expr.sum(axis=1), axis=0) * 1e6
    log_cpm = np.log1p(cpm)

    # Rename coordinates columns as 'X' and 'Y' in metadata file
    meta_raw = meta_raw.rename(columns={"center_x": "X", "center_y": "Y"})
    # Validate required columns exist
    if "X" not in meta_raw.columns or "Y" not in meta_raw.columns:
        raise ValueError("Metadata file must contain 'center_x' and 'center_y' columns.")
    meta = meta_raw[["X", "Y"]].copy()
    
    # Check if 'cell_type' exists in original metadata
    if "cell_type" in meta_raw.columns:
        meta["cell_type"] = meta_raw["cell_type"]
    else:
        # If not present, add a default column with value "All"
        meta["cell_type"] = "All"


    # Define output file paths
    counts_file = os.path.join(output_dir, f"{sample_name}_counts.txt")
    meta_file   = os.path.join(output_dir, f"{sample_name}_meta.txt")


    # Save gene expression and metadata as tab-delimited txt files
    expr.to_csv(counts_file, sep="\t", index=True, index_label=None)
    meta.to_csv(meta_file, sep="\t", index=True, index_label=None)

    print(f"Saved: {counts_file} and {meta_file}")


# In[ ]:


if __name__ == "__main__":
    gene_expr_path = input("Enter path to gene expression CSV: ").strip()
    metadata_path = input("Enter path to metadata CSV: ").strip()
    output_dir = input("Enter directory path / the new name to save results: ").strip()
    sample_name = input("Output file label:").strip()

    process_sample(gene_expr_path, metadata_path, output_dir, sample_name)

