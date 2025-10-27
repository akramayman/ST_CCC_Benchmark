#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import argparse


# ### Input mode is "csv" expression and meta data, and the output is txt files with log1 normalization

# In[ ]:


def process_csv_to_txt(gene_expr_path, metadata_path, output_dir, sample_name="sample", normalized=False):
    # Load the csv files 
    expr = pd.read_csv(gene_expr_path, index_col=0)
    meta_raw = pd.read_csv(metadata_path, index_col=0)
    expr.index = expr.index.astype(str)
    meta_raw.index = meta_raw.index.astype(str)

    # --- Normalize if not already normalized ---
    if normalized:
        print("🔄 Normalizing expression data (CPM + log1p)...")
        cpm = expr.div(expr.sum(axis=1), axis=0) * 1e6
        expr = np.log1p(cpm)
    else:
        print("✅ Skipping normalization (already normalized)")

    # Rename coordinate columns to 'X' and 'Y' depending on which exist
    if 'center_x' in meta_raw.columns and 'center_y' in meta_raw.columns:
        meta_raw = meta_raw.rename(columns={'center_x': 'X', 'center_y': 'Y'})
    elif 'x' in meta_raw.columns and 'y' in meta_raw.columns:
        meta_raw = meta_raw.rename(columns={'x': 'X', 'y': 'Y'})

    # Validate required columns exist
    if "X" not in meta_raw.columns or "Y" not in meta_raw.columns:
        raise ValueError("Metadata file must contain 'center_x' and 'center_y' columns.")
    
    # Check if 'cell_type' exists in original metadata
    if "celltype" in meta_raw.columns:
        meta_raw["cell_type"] = meta_raw["celltype"]
    elif "cell_type" in meta_raw.columns:
        meta_raw["cell_type"] = meta_raw["cell_type"]
    else:
        raise ValueError("❌ Metadata file must contain a 'cell_type' or 'celltype' column.")

    meta = meta_raw[["X", "Y", "cell_type"]].copy()
    # Define output file paths
    counts_file = os.path.join(output_dir, f"{sample_name}_counts.txt")
    meta_file   = os.path.join(output_dir, f"{sample_name}_meta.txt")

    
    expr.to_csv(counts_file, sep="\t", index=True, index_label=None)
    meta.to_csv(meta_file, sep="\t", index=True, index_label=None)

    print(f"Saved: {counts_file} and {meta_file}")


# ### Input mode is "csv" expression and meta data, and the output is h5ad files

# In[ ]:


def CSV_to_h5ad(gene_expr_path, metadata_path, output_dir, sample_name="sample"):
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


# ### Input mode is "csv" expression and meta data, and the output is csv of spaCI

# In[ ]:


def process_csv_to_csvspaCI(gene_expr_path, metadata_path, output_dir, sample_name="sample"):
    # Load the csv files 
    expr = pd.read_csv(gene_expr_path, index_col=0)
    meta_raw = pd.read_csv(metadata_path, index_col=0)
    expr.index = expr.index.astype(str)
    meta_raw.index = meta_raw.index.astype(str)

    # If expression rows (index) don't match metadata index (cells), and columns do → transpose
    if set(expr.index).intersection(meta_raw.index):
        print("🔄 Detected cells as rows and genes as columns — transposing expression matrix.")
        expr = expr.T
        expr.index = expr.index.str.upper()
    else:
        print("✅ Expression matrix orientation looks correct (genes as rows, cells as columns).")

    # --- Handle coordinate renaming ---
    if 'center_x' in meta_raw.columns and 'center_y' in meta_raw.columns:
        meta_raw = meta_raw.rename(columns={'center_x': 'X', 'center_y': 'Y'})
    elif 'x' in meta_raw.columns and 'y' in meta_raw.columns:
        meta_raw = meta_raw.rename(columns={'x': 'X', 'y': 'Y'})

    # Validate required columns exist
    if "X" not in meta_raw.columns or "Y" not in meta_raw.columns:
        raise ValueError("❌ Metadata file must contain 'center_x' and 'center_y' columns.")

    meta = meta_raw[["X", "Y"]].copy()

    # --- Check for cell type ---
    if "celltype" in meta_raw.columns:
        meta["cell_type"] = meta_raw["celltype"]
    elif "cell_type" in meta_raw.columns:
        meta["cell_type"] = meta_raw["cell_type"]
    else:
        raise ValueError("❌ Metadata file must contain a 'cell_type' or 'celltype' column.")

    meta = meta.rename(columns={
        "X": "x",
        "Y": "y",
        "cell_type": "type"
    })

    # --- Save results ---
    counts_file = os.path.join(output_dir, f"{sample_name}_counts.csv")
    meta_file   = os.path.join(output_dir, f"{sample_name}_meta.csv")

    expr.to_csv(counts_file, index=True)
    meta.to_csv(meta_file, index=True)

    print(f"✅ Saved: {counts_file} and {meta_file}")


# ### Input mode is "h5ad" data, and the output is txt files with log1 normalization 

# In[ ]:


def process_h5ad_to_txt(h5ad_path, output_dir, sample_name="sample", normalized=False):
    """
    Process a .h5ad file to extract gene expression and spatial metadata as txt for Spacia input.
    If normalized=False, perform CPM + log1p normalization; otherwise, keep counts as-is.
    """
    # Load the h5ad file
    adata = anndata.read_h5ad(h5ad_path)

    # Extract counts matrix and convert to DataFrame
    expr = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )

    # --- Normalize expression to CPM + log1p ---
     # --- Normalize if not already normalized ---
    if normalized:
        print("🔄 Normalizing expression data (CPM + log1p)...")
        cpm = expr.div(expr.sum(axis=1), axis=0) * 1e6
        expr = np.log1p(cpm)
    else:
        print("✅ Skipping normalization (already normalized)")

    # Convert index (cell IDs) to string
    expr.index = expr.index.astype(str)

    # Copy metadata
    meta = adata.obs.copy()

    # --- Detect spatial coordinates ---
    if {"X","Y"}.issubset(meta.columns):
        print("✅ Using existing X/Y in obs")
    elif {"center_x","center_y"}.issubset(meta.columns):
        print("🔄 Renaming center_x/center_y → X/Y")
        meta = meta.rename(columns={'center_x':'X','center_y':'Y'})
    elif "spatial" in adata.obsm:
        print("🔄 Extracting spatial coordinates from obsm['spatial']")
        spatial_coords = adata.obsm['spatial']
        meta['X'] = spatial_coords[:,0]
        meta['Y'] = spatial_coords[:,1]
    else:
        raise ValueError("Spatial coordinates not found in obs or obsm['spatial'].")

    # Check for cell type annotation, otherwise default to "All"
    if "celltype" in meta.columns:
        meta["cell_type"] = meta["celltype"]
    elif "cell_type" in meta.columns:
        meta["cell_type"] = meta["cell_type"]
    else:
        raise ValueError("❌ Metadata file must contain a 'cell_type' or 'celltype' column.")

    # Keep only necessary columns
    meta = meta[["X", "Y", "cell_type"]]

    # Convert meta index (cell IDs) to string explicitly
    meta.index = meta.index.astype(str)

    # Define output file paths
    counts_file = os.path.join(output_dir, f"{sample_name}_counts.txt")
    meta_file = os.path.join(output_dir, f"{sample_name}_meta.txt")

    # Remove index name to avoid printing index header
    expr.index.name = None
    meta.index.name = None

    # Save expression and metadata files
    expr.to_csv(counts_file, sep="\t", index=True, index_label=None)
    meta.to_csv(meta_file, sep="\t", index=True, index_label=None)

    print(f"✅ Saved: {counts_file} and {meta_file}")


# ### Input mode is "h5ad" data, and the output is csv

# In[ ]:


def process_h5ad_to_csv(h5ad_path, output_dir, sample_name="sample"):
    """
    Process a .h5ad file to extract gene expression and spatial metadata as csv.
    """
    # Load the h5ad file
    adata = anndata.read_h5ad(h5ad_path)

    # Extract counts matrix and convert to DataFrame
    expr = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )

    # Copy metadata
    meta = adata.obs.copy()

    # --- Detect spatial coordinates ---
    if {"X","Y"}.issubset(meta.columns):
        print("✅ Using existing X/Y in obs")
    elif {"center_x","center_y"}.issubset(meta.columns):
        print("🔄 Renaming center_x/center_y → X/Y")
        meta = meta.rename(columns={'center_x':'X','center_y':'Y'})
    elif "spatial" in adata.obsm:
        print("🔄 Extracting spatial coordinates from obsm['spatial']")
        spatial_coords = adata.obsm['spatial']
        meta['X'] = spatial_coords[:,0]
        meta['Y'] = spatial_coords[:,1]
    else:
        raise ValueError("Spatial coordinates not found in obs or obsm['spatial'].")

    # Check for cell type annotation, otherwise default to "All"
    if "celltype" in meta.columns:
        meta["cell_type"] = meta["celltype"]
    elif "cell_type" in meta.columns:
        meta["cell_type"] = meta["cell_type"]
    else:
        raise ValueError("❌ Metadata file must contain a 'cell_type' or 'celltype' column.")

    # Keep only necessary columns
    meta = meta[["X", "Y", "cell_type"]]

    # Convert meta index (cell IDs) to string explicitly
    meta.index = meta.index.astype(str)

    # Define output file paths
    counts_file = os.path.join(output_dir, f"{sample_name}_counts.csv")
    meta_file   = os.path.join(output_dir, f"{sample_name}_meta.csv")

    # Save expression and metadata files
    expr.to_csv(counts_file, sep="\t", index=True, index_label=None)
    meta.to_csv(meta_file, sep="\t", index=True, index_label=None)

    print(f"✅ Saved: {counts_file} and {meta_file}")


# ### Input mode is "h5ad" data, and the output is csv spaCI

# In[ ]:


def process_h5ad_to_csvspaCI(h5ad_path, output_dir, sample_name="sample"):
    """
    Process a .h5ad file to extract gene expression and spatial metadata as csv for spaCI.
    """
    # Load the h5ad file
    adata = anndata.read_h5ad(h5ad_path)

    # Extract counts matrix and convert to DataFrame
    expr = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )
    expr.index = expr.index.str.upper()
    # --- Check orientation: if rows = cells, columns = genes → transpose ---
    if set(expr.index).intersection(adata.obs_names):
        print("🔄 Detected cells as rows and genes as columns — transposing expression matrix.")
        expr = expr.T
    else:
        print("✅ Expression matrix orientation looks correct (genes as rows, cells as columns).")


    # Copy metadata
    meta = adata.obs.copy()

    # --- Detect spatial coordinates ---
    if {"X","Y"}.issubset(meta.columns):
        print("✅ Using existing X/Y in obs")
    elif {"center_x","center_y"}.issubset(meta.columns):
        print("🔄 Renaming center_x/center_y → X/Y")
        meta = meta.rename(columns={'center_x':'X','center_y':'Y'})
    elif "spatial" in adata.obsm:
        print("🔄 Extracting spatial coordinates from obsm['spatial']")
        spatial_coords = adata.obsm['spatial']
        meta['X'] = spatial_coords[:,0]
        meta['Y'] = spatial_coords[:,1]
    else:
        raise ValueError("Spatial coordinates not found in obs or obsm['spatial'].")

    # Check for cell type annotation, otherwise default to "All"
    if "celltype" in meta.columns:
        meta["cell_type"] = meta["celltype"]
    elif "cell_type" in meta.columns:
        meta["cell_type"] = meta["cell_type"]
    else:
        raise ValueError("❌ Metadata file must contain a 'cell_type' or 'celltype' column.")

    # Keep only necessary columns
    meta = meta[["X", "Y", "cell_type"]]
    meta = meta.rename(columns={
        "X": "x",
        "Y": "y",
        "cell_type": "type"
    })
    # Convert meta index (cell IDs) to string explicitly
    meta.index = meta.index.astype(str)

    # Define output file paths
    counts_file = os.path.join(output_dir, f"{sample_name}_counts.csv")
    meta_file   = os.path.join(output_dir, f"{sample_name}_meta.csv")

    # Save expression and metadata files
    expr.to_csv(counts_file, sep="\t", index=True, index_label=None)
    meta.to_csv(meta_file, sep="\t", index=True, index_label=None)

    print(f"✅ Saved: {counts_file} and {meta_file}")


# ## Main function for convert the files depend on the selections

# In[ ]:


def convert_file(input_format, output_format, gene_expr_path=None, metadata_path=None, h5ad_path=None,
                 output_dir="./output", sample_name="sample",normalized=False):
    """
    General dispatcher for converting spatial transcriptomics data between formats:
      - from csv → h5ad / txt / csv_spaci
      - from h5ad → csv / txt / csv_spaci
    """
    # Normalize formats input and output
    from_format = input_format.lower()
    to_format = output_format.lower()
    os.makedirs(output_dir, exist_ok=True)

    # { CSV → TXT }
    if from_format == "csv" and to_format == "txt":
        process_csv_to_txt(gene_expr_path, metadata_path, output_dir, sample_name,normalized)
        print(f"✅ Saved txt file for {sample_name}")

    # { CSV → H5AD }
    elif from_format == "csv" and to_format == "h5ad":
        CSV_to_h5ad(gene_expr_path, metadata_path, output_dir, sample_name)
        print(f"✅ Saved h5ad file for {sample_name}")

    # { CSV → CSV SpaCI }
    elif from_format == "csv" and to_format == "csv_spaci":
        process_csv_to_csvspaCI(gene_expr_path, metadata_path, output_dir, sample_name)
        print(f"✅ Saved SpaCI-compatible CSV files for {sample_name}")

    # { H5AD → TXT }
    elif from_format == "h5ad" and to_format == "txt":
        process_h5ad_to_txt(h5ad_path, output_dir, sample_name,normalized)
        print(f"✅ Saved txt file for {sample_name}")

    # { H5AD → CSV }
    elif from_format == "h5ad" and to_format == "csv":
        process_h5ad_to_csv(h5ad_path, output_dir, sample_name)
        print(f"✅ Saved csv file for {sample_name}")

    # { H5AD → CSV SpaCI }
    elif from_format == "h5ad" and to_format == "csv_spaci":
        process_h5ad_to_csvspaCI(h5ad_path, output_dir, sample_name)
        print(f"✅ Saved SpaCI-compatible CSV files for {sample_name}")

    else:
        raise ValueError(f"Unsupported conversion: {input_format} → {output_format}")


# # Main function

# In[ ]:


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert spatial transcriptomics data between formats.")
    parser.add_argument("--from", dest="input_format", required=True, choices=["csv", "h5ad"], help="Input format")
    parser.add_argument("--to", dest="output_format", required=True, choices=["h5ad", "txt", "csv_spaCI", "csv"], help="Output format")
    parser.add_argument("--expr", dest="gene_expr_path", help="Path to gene expression CSV")
    parser.add_argument("--meta", dest="metadata_path", help="Path to metadata CSV")
    parser.add_argument("--h5ad", dest="h5ad_path", help="Path to H5AD file")
    parser.add_argument("--out", dest="output_dir", default="./output", help="Output directory (default: ./output)")
    parser.add_argument("--name", dest="sample_name", default="sample", help="Sample name prefix")

    parser.add_argument(
    "--normalized",
    dest="normalized",
    action="store_true",
    help="Used when the input is not normalized (only relevant for TXT output), -log1p(cpm) normalization-"
    )
    args = parser.parse_args()

    convert_file(
        input_format=args.input_format,
        output_format=args.output_format,
        gene_expr_path=args.gene_expr_path,
        metadata_path=args.metadata_path,
        h5ad_path=args.h5ad_path,
        output_dir=args.output_dir,
        sample_name=args.sample_name,
        normalized=args.normalized
    )

