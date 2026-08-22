# cd /Users/rzhuang/Documents/PycharmProjects/deconvATAC
# source .venv/bin/activate
# python run_cell2location.py

import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')

import scanpy as sc
import numpy as np
import pandas as pd
from deconvatac.tl import cell2location, jsd, rmse
from deconvatac.pp import highly_variable_peaks

print("Loading data...")
russell_st = sc.read_h5ad("data/example_notebooks/cell2location/russell_250_atac.h5ad")
russell_sc = sc.read_h5ad("data/example_notebooks/cell2location/russel_ref_atac.h5ad")
# russell_st: spatial ATAC (360 spots x 53451 peaks)
#   .obsm['proportions'] = ground-truth cell-type fractions
#   .obsm['spatial']     = spatial coordinates
#   .layers['log_norm'], .layers['tfidf_normalized']
#
# russell_sc: reference scATAC (2535 cells x 53451 peaks)
#   .obs['cell_type']    = cell-type labels
#   .layers['log_norm'], .layers['tfidf_normalized']
# --- Feature selection (subset to 20k highly variable peaks) ---
print(f"  Spatial: {russell_st.shape[0]} spots x {russell_st.shape[1]} peaks")
print(f"  Reference: {russell_sc.shape[0]} cells x {russell_sc.shape[1]} peaks")

print("Selecting highly variable peaks...")
highly_variable_peaks(adata=russell_sc, cluster_key="cell_type")
# russell_sc.var["highly_variable"] is a boolean vector of length n_features, 
# True for highly variable features, False otherwise
russell_st = russell_st[:, russell_sc.var["highly_variable"]].copy()
russell_sc = russell_sc[:, russell_sc.var["highly_variable"]].copy()

print(f"russell_st shape: {russell_st.shape[0]} x {russell_st.shape[1]}")
print(f"russell_sc shape: {russell_sc.shape[0]} x {russell_sc.shape[1]}")
print(f"  After HVP selection: {russell_st.shape[1]} peaks")

n_cell_types = russell_sc.obs["cell_type"].nunique()
print(f"Number of cell types in russell_sc: {n_cell_types}")


print("Running Cell2Location (this may take a while)...")
russell_st, russell_sc = cell2location(
    adata_spatial=russell_st,
    adata_ref=russell_sc,
    N_cells_per_location=8,
    detection_alpha=20,
    labels_key="cell_type",  # cell type key in russell_sc.obs, there are 10 types of cells
    use_gpu=False,
    max_epochs_spatial=400,
    max_epochs_ref=400,
    return_adatas=True,
    plots=False,
)
print("Cell2Location complete!")

print("Saving visualization...")
russell_st.obs[russell_st.uns["mod"]["factor_names"]] = russell_st.obsm["q05_cell_abundance_w_sf"]
sc.pl.embedding(russell_st, basis="spatial", color=russell_st.uns["mod"]["factor_names"], show=False)
import matplotlib.pyplot as plt
plt.savefig("cell2location_results/deconvolution_plot.png", dpi=150, bbox_inches="tight")
plt.close()
print("Done! Results saved to cell2location_results/")

# --- Calculate metrics ---
# print("Calculating metrics...")
# targets = pd.DataFrame(russell_st.obsm["proportions"], columns=russell_st.uns["proportion_names"], index=russell_st.obs_names)
# abundances = russell_st.obsm["q05_cell_abundance_w_sf"]
# predictions = abundances.div(abundances.sum(axis=1), axis=0)
# predictions = predictions.loc[targets.index, targets.columns]
# print(f"  JSD:  {jsd(predictions, targets)}")
# print(f"  RMSE: {rmse(predictions, targets)}")
