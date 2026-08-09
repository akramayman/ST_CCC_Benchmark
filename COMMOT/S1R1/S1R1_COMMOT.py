#!/usr/bin/env python3
import commot as ct
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np


# ---------------- paths ----------------
h5ad_path = "/home/akram/share/data/Image_Anndata/S1R1_Full_Atlas.h5ad"
output_dir = Path("/home/akram/share/Result/COMMOT_S1R1")
output_dir.mkdir(exist_ok=True, parents=True)

# ---------------- data -----------------
adata = sc.read_h5ad(h5ad_path)
print(adata)

# ---------------- plots ----------------
sc.pl.scatter(adata, x="X", y="Y", color="cell_type", size=1)
plt.savefig(output_dir / "spatial_cell_type_Full.png", dpi=200, bbox_inches="tight")
plt.close()


# ---------------- normalize ------------
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

print("Normalization anf lop is done")

# ---------------- commot ---------------
df_ligrec = ct.pp.ligand_receptor_database(species="mouse",
                                           signaling_type="Secreted Signaling",
                                           database="CellChat")
df_cellchat_filtered = ct.pp.filter_lr_database(df_ligrec, adata, min_cell_pct=0.05)
df_cellchat_filtered.to_csv(output_dir / "lr_filtered_Full.csv", index=False)
print("LR fitered is done and spatial communication is starting")

print("COMMOT required obsm .spatial contain the coordinate")
# Make sure coordinates are numeric
adata.obs['X'] = adata.obs['X'].astype(float)
adata.obs['Y'] = adata.obs['Y'].astype(float)

# Put into .obsm['spatial']
adata.obsm['spatial'] = np.vstack([adata.obs['X'], adata.obs['Y']]).T

print("create the spatial section is done")
ct.tl.spatial_communication(adata,
    database_name='cellchat',
    df_ligrec=df_cellchat_filtered,
    dis_thr=500,
    heteromeric=True,
    pathway_sum=True)

adata.write(output_dir / "adata_spatial_communication_Full.h5ad")
print("All done! Results saved in", output_dir)
