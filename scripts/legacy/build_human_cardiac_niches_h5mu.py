from pathlib import Path

import scanpy as sc
import muon as mu
import deconvatac as dv

ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIR = ROOT / "data" / "raw" / "references" / "human_cardiac_niches"

adata_rna = sc.read_h5ad(REFERENCE_DIR / "rna" / "reference.h5ad", backed="r")
adata_atac = sc.read_h5ad(REFERENCE_DIR / "atac" / "reference.h5ad")

print(adata_rna.shape) # (704296, 32732)
print(adata_atac.shape) # (139835, 429828)

exit()
# Process ATAC
dv.pp.reads_to_fragments(adata_atac)
adata_atac.X = adata_atac.layers["fragments"].copy()

dv.pp.highly_variable_peaks(adata_atac, cluster_key="cell_type")
dv.pp.highly_accessible_peaks(adata_atac)
mu.atac.pp.tfidf(adata_atac, to_layer="tfidf_normalized")

adata_atac.layers["counts"] = adata_atac.X.copy()
sc.pp.normalize_total(adata_atac)
sc.pp.log1p(adata_atac)
adata_atac.layers["log_norm"] = adata_atac.X.copy()
adata_atac.X = adata_atac.layers["counts"].copy()
del adata_atac.layers["counts"]

# Keep only RNA cells that also have ATAC
adata_rna = adata_rna[adata_rna.obs_names.isin(adata_atac.obs_names)].to_memory()

# Process RNA
adata_rna.layers["counts"] = adata_rna.X.copy()
sc.pp.normalize_total(adata_rna)
sc.pp.log1p(adata_rna)
sc.pp.highly_variable_genes(adata_rna, n_top_genes=4000)
adata_rna.layers["log_norm"] = adata_rna.X.copy()
adata_rna.X = adata_rna.layers["counts"].copy()
del adata_rna.layers["counts"]

# Combine
mdata = mu.MuData({"atac": adata_atac, "rna": adata_rna})
mdata.write(REFERENCE_DIR / "human_cardiac_niches.h5mu")
