import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')

import os
import mudata as mu
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from deconvatac.tl import rctd, jsd, rmse
from deconvatac.pp import highly_variable_peaks

RESULTS_PATH = "./rctd_results"
R_LIB_PATH = os.path.expanduser("~/R/library")

# --- Load data ---
print("Loading data...")
heart_st = mu.read_h5mu("data/example_notebooks/simulation/Heart_heterogeneous_4zones.h5mu").mod["atac"]
heart_sc = sc.read_h5ad("data/example_notebooks/rctd/human_cardiac_niches_atac.h5ad")
print(f"  Spatial: {heart_st.shape[0]} spots x {heart_st.shape[1]} peaks")
print(f"  Reference: {heart_sc.shape[0]} cells x {heart_sc.shape[1]} peaks")

# --- Visualize ground truth ---
print("Plotting ground truth...")
heart_st.obs = heart_st.obs.reset_index().join(
    pd.DataFrame(heart_st.obsm["proportions"], columns=heart_st.uns["proportion_names"]).reset_index(drop=True))
heart_st.obs['ground_truth'] = heart_st.obs.iloc[:, 2:].idxmax(axis=1)
sc.pl.embedding(heart_st, basis="spatial", color=['cell_count', 'ground_truth'], show=False)
os.makedirs(RESULTS_PATH, exist_ok=True)
plt.savefig(f"{RESULTS_PATH}/ground_truth.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {RESULTS_PATH}/ground_truth.png")

# --- Feature selection ---
print("Selecting highly variable peaks...")
highly_variable_peaks(adata=heart_sc, cluster_key="cell_type")
heart_st = heart_st[:, heart_sc.var["highly_variable"]].copy()
heart_sc = heart_sc[:, heart_sc.var["highly_variable"]].copy()
print(f"  After HVP selection: {heart_st.shape[1]} peaks")

# --- Run RCTD ---
print("Running RCTD (this may take a while)...")
rctd(
    adata_spatial=heart_st,
    adata_ref=heart_sc,
    labels_key="cell_type",
    r_lib_path=R_LIB_PATH,
    doublet_mode='full',
    create_rctd_kwargs={
        "CELL_MIN_INSTANCE": 0,
        "gene_cutoff": 0,
        "fc_cutoff": 0,
        "gene_cutoff_reg": 0,
        "fc_cutoff_reg": 0,
        "UMI_min": 0,
    },
)
print("RCTD complete!")

# --- Load and attach results ---
print("Visualizing results...")
deconv_results = pd.read_csv(f"{RESULTS_PATH}/estimated_proportions.csv", index_col=0)
deconv_results.index = heart_st.obs.index
heart_st.obsm["rctd_proportions"] = deconv_results

# --- Visualize RCTD vs ground truth ---
max_prob_cluster = np.argmax(heart_st.obsm["rctd_proportions"], axis=1)
cluster_id = deconv_results.columns.to_numpy()
heart_st.obs["rctd_max_abundant"] = cluster_id[max_prob_cluster]
heart_st.obs["rctd_max_abundant"] = pd.Categorical(
    heart_st.obs["rctd_max_abundant"],
    categories=heart_sc.obs.cell_type.cat.categories,
)
heart_st.uns["rctd_max_abundant_colors"] = heart_st.uns["ground_truth_colors"].copy()
sc.pl.embedding(heart_st, basis="spatial", color=["rctd_max_abundant", "ground_truth"], wspace=0.4, show=False)
plt.savefig(f"{RESULTS_PATH}/rctd_vs_ground_truth.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved {RESULTS_PATH}/rctd_vs_ground_truth.png")

# --- Calculate metrics ---
print("Calculating metrics...")
targets = pd.DataFrame(heart_st.obsm["proportions"], columns=heart_st.uns["proportion_names"], index=heart_st.obs_names)
predictions = heart_st.obsm["rctd_proportions"].loc[targets.index, targets.columns]
print(f"  JSD:  {jsd(predictions, targets)}")
print(f"  RMSE: {rmse(predictions, targets)}")
print("Done!")