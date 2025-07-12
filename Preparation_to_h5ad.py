#!/usr/bin/env python
# coding: utf-8

# # Preprocessing to create H5ad files from csv

# In[2]:


import os
import pandas as pd
import scanpy as sc


# ### This fuction to read the gene expression and metadata CSV files

# In[ ]:


def processing(gene_expr_path, metadata_path, output_dir, sample_name="sample"):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # CSVs Paths
    counts_file    = pd.read_csv(gene_expr_path, index_col=0)
    meta_file      = pd.read_csv(metadata_path, index_col=0)

    # AnnData
    adata = sc.AnnData(
        X=counts_file.values,
        obs=meta_file,
        var=pd.DataFrame(index=counts_file.columns)
    )

    # Add spatial coordinates
    if {'center_x','center_y'}.issubset(meta_file.columns):
        adata.obsm['spatial'] = meta_file[['center_x','center_y']].values

    # Save
    out_path = os.path.join(output_dir, f"{sample_name}.h5ad")
    adata.write_h5ad(out_path)
    print(f"Saved {out_path}: {adata.n_obs} cells × {adata.n_vars} genes")


# In[ ]:


if __name__ == "__main__":
    gene_expr_path = input("Enter path to gene expression CSV: ").strip()
    metadata_path = input("Enter path to metadata CSV: ").strip()
    output_dir = input("Enter directory path / the new name to save results: ").strip()
    sample_name = input("Output file label:").strip()

    processing(gene_expr_path, metadata_path, output_dir, sample_name)

