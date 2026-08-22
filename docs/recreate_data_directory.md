# Recreate the `data/` directory on another machine

This guide rebuilds the current operational `data/` tree from public downloads and repository scripts. It was verified on **2026-08-22**. Run every command from the repository root.

The completed directory is about **29.4 GiB**. Budget at least **40 GiB of free disk space** for the final files plus staged ZIP archives. The downloads themselves total about 18.5 GB; the Zenodo notebook archive expands from 1.85 GB to about 7.7 GiB.

## Recovery status

| Directory | Recovery route | Fidelity |
|---|---|---|
| `data/raw/references/human_cardiac_niches/` | Heart Cell Atlas direct downloads | Byte-exact |
| `data/example_notebooks/` | Zenodo record `15089738` | Authoritative archive; logically exact |
| `data/raw/sources/10x_genomics/` | 10x Genomics CDN | Byte-exact |
| `data/raw/sources/snapatac2/` | Official SnapATAC2/scverse mirror | Labels are exact; see the PBMC caveat below |
| `data/raw/references/russell_250/` | Zenodo source plus repository feature selection | Reproducible with the pinned environment |
| `data/raw/references/pbmc_*` | Repository preparation scripts | Reproducible; reference H5ADs were verified byte-exact |
| `data/processed/datasets/` | Repository simulation/preparation scripts | Scientifically reproducible; HDF5 bytes can depend on package versions |
| `data/archive/legacy_results/deconvolution_results/` | Zenodo record `15089738` | Byte-exact |
| Legacy Cell2location/RCTD top-level results | Transfer from the original machine or regenerate | Not publicly downloadable |

The top-level `data/` directory is intentionally git-ignored, including its README and manifest files. A repository clone therefore does not contain any of these files. This guide recreates the operational data and the metadata needed by the unified runner; copy documentation-only files from the original machine if an exact directory snapshot is required.

> **Safety:** Run this workflow against a fresh clone or an empty/new `data/` directory. It initializes `data/registry/datasets.yaml`, extracts archives with overwrite enabled, and invokes generation scripts with `--overwrite`. Back up an existing `data/` directory before using these commands on it.

## 1. Prepare the repository and environment

Clone or copy the repository at the same commit used on the original machine, then create the environment:

```bash
cd /path/to/deconvATAC

python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[simulation]"

export PYTHONPATH="$PWD/src"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/deconvatac-mpl"
export DECONVATAC_DOWNLOAD_DIR="${TMPDIR:-/tmp}/deconvatac-data-downloads"

mkdir -p "$MPLCONFIGDIR" "$DECONVATAC_DOWNLOAD_DIR" data/registry
printf '{}\n' > data/registry/datasets.yaml
```

The normal project install includes `rpy2` and therefore expects a working R installation. If that install fails because R is unavailable and the goal is only to rebuild and validate the maintained data, use this data-only environment instead:

```bash
python -m pip install \
  numpy==1.25.2 pandas==2.3.3 scipy==1.11.3 \
  anndata==0.9.2 scanpy==1.9.5 muon==0.1.5 h5py==3.16.0 \
  matplotlib PyYAML==6.0.3 session-info \
  igraph==1.0.0 leidenalg==0.12.0 pyarrow
python -m pip install -e . --no-deps
```

This fallback is sufficient for the downloads, reference preparation, simulations, and validation below. It does not support the optional legacy RCTD regeneration in section 12.

The current files were produced with these important versions:

```text
Python      3.11.15
NumPy       1.25.2
Pandas      2.3.3
SciPy       1.11.3
AnnData     0.9.2
Scanpy      1.9.5
Muon        0.1.5
h5py        3.16.0
igraph      1.0.0
leidenalg   0.12.0
PyYAML      6.0.3
```

Use those versions if byte-level agreement of regenerated HDF5 files is important. Later compatible versions should still produce equivalent scientific inputs, but PCA, Leiden clustering, tied feature rankings, YAML dates, or HDF5 serialization can differ.

Define a resumable download helper:

```bash
fetch_data() {
  destination="$1"
  url="$2"
  mkdir -p "$(dirname "$destination")"
  curl --fail --location \
    --retry 5 --retry-all-errors \
    --continue-at - \
    --output "$destination" \
    "$url"
}
```

If a server rejects resuming a damaged partial file, remove only that partial destination and rerun the corresponding command without `--continue-at -`.

## 2. Download the canonical Heart references

These are the original Heart Cell Atlas v2 files, renamed into the repository's canonical layout. No transformation is applied.

Source collection: [Spatially resolved multiomics of human cardiac niches](https://cellxgene.cziscience.com/collections/3116d060-0a8e-4767-99bb-e866badea1ed).

```bash
fetch_data \
  data/raw/references/human_cardiac_niches/atac/reference.h5ad \
  https://cellgeni.cog.sanger.ac.uk/heartcellatlas/v2/Adult_Peaks.h5ad

fetch_data \
  data/raw/references/human_cardiac_niches/rna/reference.h5ad \
  https://cellgeni.cog.sanger.ac.uk/heartcellatlas/v2/Global_raw.h5ad
```

Expected files:

| Path | Bytes | SHA-256 |
|---|---:|---|
| `atac/reference.h5ad` | 7,252,088,744 | `f5fd7c42cdc89ac19879f852b0635a8c6a59b4bb75bf1e9c2d5bc354041cd442` |
| `rna/reference.h5ad` | 9,125,578,296 | `2280a34361a67ac6b24a787ddc037d0f7f78d7db91f5dcfb8ee822d7e80c7ed2` |

Create the small reference manifest:

```bash
.venv/bin/python - <<'PY'
from pathlib import Path

path = Path("data/raw/references/human_cardiac_niches/reference.yaml")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(
    """reference_id: human_cardiac_niches
description: Human cardiac niches single-cell reference atlas used for regenerated Heart simulations.
labels_key: cell_type
modalities:
  atac:
    path: data/raw/references/human_cardiac_niches/atac/reference.h5ad
    source_filename: Adult_Peaks.h5ad
  rna:
    path: data/raw/references/human_cardiac_niches/rna/reference.h5ad
    source_filename: Global_raw.h5ad
"""
)
PY
```

## 3. Download the published Zenodo archives

Use the pinned publication record [Zenodo 15089738](https://zenodo.org/records/15089738), DOI `10.5281/zenodo.15089738`. Do not replace the record number with a moving `latest` URL.

```bash
fetch_data \
  "$DECONVATAC_DOWNLOAD_DIR/example_notebooks.zip" \
  https://zenodo.org/api/records/15089738/files/example_notebooks.zip/content

fetch_data \
  "$DECONVATAC_DOWNLOAD_DIR/deconvolution_results.zip" \
  https://zenodo.org/api/records/15089738/files/deconvolution_results.zip/content
```

Verify the archives before extracting them:

```bash
.venv/bin/python - <<'PY'
import hashlib
import os
from pathlib import Path

root = Path(os.environ["DECONVATAC_DOWNLOAD_DIR"])
expected = {
    root / "example_notebooks.zip": "c0e961e518e1f5c5dadc879792565d3b",
    root / "deconvolution_results.zip": "f28734f2d98c98906c88bb90e001b8d5",
}

for path, wanted in expected.items():
    with path.open("rb") as handle:
        observed = hashlib.file_digest(handle, "md5").hexdigest()
    if observed != wanted:
        raise SystemExit(f"MD5 mismatch for {path}: {observed} != {wanted}")
    print(f"OK  {observed}  {path}")
PY
```

Extract them into their expected locations:

```bash
unzip -q -o "$DECONVATAC_DOWNLOAD_DIR/example_notebooks.zip" -d data

mkdir -p data/archive/legacy_results
unzip -q -o "$DECONVATAC_DOWNLOAD_DIR/deconvolution_results.zip" \
  -d data/archive/legacy_results
```

The first archive creates:

```text
data/example_notebooks/
  cell2location/russel_ref_atac.h5ad
  cell2location/russell_250_atac.h5ad
  rctd/human_cardiac_niches_atac.h5ad
  simulation/Heart_heterogeneous_4zones.h5mu
  simulation/Heart_homogeneous_1zone.h5mu
```

The three H5AD files are byte-identical to the local copies. The two H5MU files in this workspace were later reserialized, so their bytes differ from Zenodo, but an exhaustive HDF5 comparison found identical groups, attributes, dtypes, shapes, and values.

The second archive creates `data/archive/legacy_results/deconvolution_results/`; all 249 result files were verified byte-identical to the current local archive.

## 4. Download the 10x PBMC Multiome sources

Official dataset page: [PBMC from a healthy donor, granulocytes removed through cell sorting (10k), Cell Ranger ARC 2.0.0](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0).

```bash
PBMC_10X_DIR=data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0
PBMC_10X_BASE=https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k

fetch_data \
  "$PBMC_10X_DIR/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5" \
  "$PBMC_10X_BASE/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"

fetch_data \
  "$PBMC_10X_DIR/pbmc_granulocyte_sorted_10k_atac_peaks.bed" \
  "$PBMC_10X_BASE/pbmc_granulocyte_sorted_10k_atac_peaks.bed"

fetch_data \
  "$PBMC_10X_DIR/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv" \
  "$PBMC_10X_BASE/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv"

fetch_data \
  "$PBMC_10X_DIR/pbmc_granulocyte_sorted_10k_web_summary.html" \
  "$PBMC_10X_BASE/pbmc_granulocyte_sorted_10k_web_summary.html"
```

The large fragments file is intentionally not part of the current directory and is not required for these peak-matrix benchmarks.

Expected checksums:

| File | Bytes | SHA-256 |
|---|---:|---|
| `filtered_feature_bc_matrix.h5` | 192,125,528 | `f6824171378787baab244f559b8b438f79db2eb39f78d17b2196f7ecd2c03549` |
| `atac_peaks.bed` | 3,439,145 | `3975a4057f9caa3fb69ddaecc6ae9e530e77551717a1464c2d93ac9d73cb60ab` |
| `per_barcode_metrics.csv` | 88,822,764 | `fd3e069b83e152145af234667b419c982968aca0df322a92adb71284d0d902cd` |
| `web_summary.html` | 6,227,259 | `4f443bc6908c3345326cac73e11d7b16e0adc279f8023de2867e3f5d85f86ec5` |

## 5. Download the SnapATAC2 PBMC label source

Download the current official SnapATAC2/scverse object:

```bash
fetch_data \
  data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad \
  https://exampledata.scverse.org/snapatac2/10x-Multiome-Pbmc10k-RNA.h5ad
```

Expected source:

```text
bytes:   45,486,209
sha256:  a25327acff48b20b295c12221a84fd00f8f3f486ff3e7bd090fdef241b996a22
shape:   9,631 cells x 29,095 genes
labels:  19 cell-type categories
```

### PBMC provenance caveat

The original machine contains a historical 23,147,214-byte, 3,000-HVG derivative of this RNA object with SHA-256 `72505189d3c03cdebc299e31e80aafed29af538dcef34a6708f76b7998ed6dcb`. Its old Mendeley draft URL now returns HTTP 403, so that exact serialization is not publicly recoverable.

The official object above has the same 9,631 barcodes in the same order and exactly the same `cell_type` labels. It produces a byte-identical `cell_type_mapping.csv` and byte-identical downstream PBMC ATAC/RNA references. It is therefore the recommended functional replacement. Transfer the historical 23 MB file from the original machine only if its exact bytes or 3,000-gene contents are independently important.

## 6. Verify all directly downloaded inputs

Hashing the two large Heart references can take a few minutes:

```bash
.venv/bin/python - <<'PY'
import hashlib
import os
from pathlib import Path

download_root = Path(os.environ["DECONVATAC_DOWNLOAD_DIR"])
checks = [
    ("sha256", Path("data/raw/references/human_cardiac_niches/atac/reference.h5ad"), "f5fd7c42cdc89ac19879f852b0635a8c6a59b4bb75bf1e9c2d5bc354041cd442"),
    ("sha256", Path("data/raw/references/human_cardiac_niches/rna/reference.h5ad"), "2280a34361a67ac6b24a787ddc037d0f7f78d7db91f5dcfb8ee822d7e80c7ed2"),
    ("sha256", Path("data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"), "f6824171378787baab244f559b8b438f79db2eb39f78d17b2196f7ecd2c03549"),
    ("sha256", Path("data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/pbmc_granulocyte_sorted_10k_atac_peaks.bed"), "3975a4057f9caa3fb69ddaecc6ae9e530e77551717a1464c2d93ac9d73cb60ab"),
    ("sha256", Path("data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv"), "fd3e069b83e152145af234667b419c982968aca0df322a92adb71284d0d902cd"),
    ("sha256", Path("data/raw/sources/10x_genomics/pbmc_granulocyte_sorted_10k/cellranger_arc_2.0.0/pbmc_granulocyte_sorted_10k_web_summary.html"), "4f443bc6908c3345326cac73e11d7b16e0adc279f8023de2867e3f5d85f86ec5"),
    ("sha256", Path("data/raw/sources/snapatac2/pbmc10k_multiome/rna.h5ad"), "a25327acff48b20b295c12221a84fd00f8f3f486ff3e7bd090fdef241b996a22"),
    ("md5", download_root / "example_notebooks.zip", "c0e961e518e1f5c5dadc879792565d3b"),
    ("md5", download_root / "deconvolution_results.zip", "f28734f2d98c98906c88bb90e001b8d5"),
]

for algorithm, path, wanted in checks:
    with path.open("rb") as handle:
        observed = hashlib.file_digest(handle, algorithm).hexdigest()
    if observed != wanted:
        raise SystemExit(f"{algorithm} mismatch for {path}: {observed} != {wanted}")
    print(f"OK  {algorithm:6s}  {path}")
PY
```

Do not begin regeneration if any checksum fails.

## 7. Rebuild the Russell reference and processed dataset

The Zenodo archive supplies the original Russell reference and spatial example. The current canonical files add two 20,000-feature boolean annotations and export truth to CSV.

Run the complete rebuild and metadata bootstrap:

```bash
PYTHONPATH=src MPLCONFIGDIR="$MPLCONFIGDIR" .venv/bin/python - <<'PY'
from pathlib import Path

import anndata as ad
import pandas as pd
import yaml

from deconvatac.pp import highly_accessible_peaks, highly_variable_peaks

reference = ad.read_h5ad(
    "data/example_notebooks/cell2location/russel_ref_atac.h5ad"
)
spatial = ad.read_h5ad(
    "data/example_notebooks/cell2location/russell_250_atac.h5ad"
)

highly_variable_peaks(
    reference,
    cluster_key="cell_type",
    n_top_features=20000,
)
highly_accessible_peaks(
    reference,
    n_top_features=20000,
)

for column in ("highly_variable", "highly_accessible"):
    selected = reference.var_names[reference.var[column].astype(bool)]
    spatial.var[column] = spatial.var_names.isin(selected)

reference_path = Path(
    "data/raw/references/russell_250/atac/reference.h5ad"
)
spatial_path = Path(
    "data/processed/datasets/russell_250/atac/spatial.h5ad"
)
truth_path = Path(
    "data/processed/datasets/russell_250/truth/proportions.csv"
)
reference_manifest_path = Path(
    "data/raw/references/russell_250/reference.yaml"
)
dataset_config_path = Path(
    "data/processed/datasets/russell_250/dataset.yaml"
)
registry_path = Path("data/registry/datasets.yaml")

for path in (
    reference_path,
    spatial_path,
    truth_path,
    reference_manifest_path,
    dataset_config_path,
    registry_path,
):
    path.parent.mkdir(parents=True, exist_ok=True)

reference.write_h5ad(reference_path)
spatial.write_h5ad(spatial_path)

pd.DataFrame(
    spatial.obsm["proportions"],
    index=spatial.obs_names,
    columns=spatial.uns["proportion_names"],
).to_csv(truth_path)

reference_manifest = {
    "reference_id": "russell_250",
    "description": "Russell ATAC single-cell reference used by the Russell spatial ATAC simulation example.",
    "labels_key": "cell_type",
    "modalities": {
        "atac": {
            "path": str(reference_path),
            "source_filename": "russel_ref_atac.h5ad",
            "notes": "Contains processed feature annotations used by the Russell dataset feature sets.",
        }
    },
}
reference_manifest_path.write_text(
    yaml.safe_dump(reference_manifest, sort_keys=False)
)

dataset_config = {
    "dataset_id": "russell_250",
    "source": "processed_feature_annotations",
    "description": "Russell ATAC simulation example already present in this workspace.",
    "labels_key": "cell_type",
    "spatial_key": "spatial",
    "modalities": {
        "atac": {
            "reference": {"path": str(reference_path)},
            "spatial": {"path": str(spatial_path)},
            "labels_key": "cell_type",
            "spatial_key": "spatial",
            "truth": {"path": str(truth_path)},
            "feature_sets": {
                "highly_variable": {"var_column": "highly_variable"},
                "highly_accessible": {"var_column": "highly_accessible"},
                "all": {"mode": "all"},
            },
        }
    },
}
dataset_config_path.write_text(
    yaml.safe_dump(dataset_config, sort_keys=False)
)

registry = yaml.safe_load(registry_path.read_text()) or {}
registry["russell_250"] = {"config": str(dataset_config_path)}
registry_path.write_text(yaml.safe_dump(registry, sort_keys=False))

print(reference_path)
print(spatial_path)
print(truth_path)
PY
```

With the pinned environment, the rebuilt reference was verified as:

```text
bytes:   98,712,149
sha256:  34094c568e0267b3f16a98401059439679934458f7ee35c43bcda7914d578bf8
shape:   2,535 cells x 53,451 peaks
```

Do not use `scripts/prepare_feature_annotations.py` as a substitute on a clean tree: it expects the target spatial file to exist already and creates a different processed layout.

## 8. Regenerate all four Heart simulations

```bash
PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH=src \
MPLCONFIGDIR="$MPLCONFIGDIR" \
.venv/bin/python scripts/regenerate_heart_simulations.py \
  --datasets \
    human_cardiac_niches_sim_1zone_3ct_low_density \
    human_cardiac_niches_sim_1zone_10ct \
    human_cardiac_niches_sim_4zone_stripes \
    human_cardiac_niches_sim_4zone_circles \
  --overwrite
```

This writes each dataset's ATAC/RNA `spatial.h5ad`, truth CSV, source-cell JSONL, `dataset.yaml`, and registry entry. It creates feature directories but does not populate them.

Compute the three shared Heart feature lists once, then copy them into all four datasets:

```bash
export DECONVATAC_HEART_FEATURE_DIR="${TMPDIR:-/tmp}/deconvatac-heart-features"

PYTHONPATH=src .venv/bin/python scripts/prepare_shared_feature_sets.py \
  --dataset human_cardiac_niches_sim_1zone_10ct \
  --modalities atac rna \
  --feature-set-id human_cardiac_niches \
  --output-root "$DECONVATAC_HEART_FEATURE_DIR" \
  --n-top-features 20000 \
  --chunk-size 2048 \
  --overwrite

for dataset_id in \
  human_cardiac_niches_sim_1zone_3ct_low_density \
  human_cardiac_niches_sim_1zone_10ct \
  human_cardiac_niches_sim_4zone_stripes \
  human_cardiac_niches_sim_4zone_circles
do
  cp "$DECONVATAC_HEART_FEATURE_DIR/human_cardiac_niches/atac/highly_variable.txt" \
    "data/processed/datasets/$dataset_id/atac/features/"
  cp "$DECONVATAC_HEART_FEATURE_DIR/human_cardiac_niches/atac/highly_accessible.txt" \
    "data/processed/datasets/$dataset_id/atac/features/"
  cp "$DECONVATAC_HEART_FEATURE_DIR/human_cardiac_niches/rna/highly_variable.txt" \
    "data/processed/datasets/$dataset_id/rna/features/"
done
```

Important details:

- Do not pass `--skip-leiden` if matching the current schemas.
- The script resets the NumPy seed to zero for each dataset.
- The requested `num_spots: 1000` becomes a 31 by 31 grid, so each current Heart spatial file has 961 spots.
- Each feature list contains 20,000 names. The corresponding lists are byte-identical across all four Heart datasets.

## 9. Build PBMC labels and canonical references

Run these commands in order:

```bash
PYTHONPATH=src .venv/bin/python scripts/prepare_pbmc_multiome_labels.py \
  --overwrite

PYTHONPATH=src .venv/bin/python scripts/prepare_pbmc_multiome_references.py \
  --overwrite
```

Do **not** pass `--allow-celltypist-fallback` for the canonical build. The expected results are:

```text
10x input cells:        11,898
labeled retained cells:  9,627
dropped cells:            2,271
ATAC features:           143,887
RNA features:             36,601
```

The official replacement SnapATAC2 source produces the same mapping and these byte-exact references:

| Path | Bytes | SHA-256 |
|---|---:|---|
| `pbmc_.../atac/reference.h5ad` | 720,924,465 | `cdaefffbfd5dd3cb36318f68158d4bec0df1d2b2bf63562b14065f1055cf6ee6` |
| `pbmc_.../rna/reference.h5ad` | 156,284,068 | `c78cae2f49fcf31a84af8fc2c75a397c1467c7eb097f08073c8375bba3fd9458` |

The generated source manifest will record the larger official SnapATAC2 file and the current preparation date, so that small YAML file will not be byte-identical to the historical manifest.

## 10. Regenerate both PBMC simulations

Run both dataset IDs together and in this exact order:

```bash
PYTHONPATH=src .venv/bin/python scripts/regenerate_pbmc_simulations.py \
  --datasets \
    pbmc_granulocyte_sorted_10k_sim_equal_celltype \
    pbmc_granulocyte_sorted_10k_sim_observed_abundance \
  --num-spots 1024 \
  --mean-cells-per-spot 10 \
  --min-cell-type-cells 100 \
  --n-top-features 20000 \
  --chunk-size 2048 \
  --seed 0 \
  --overwrite
```

The script adds the dataset index to the base seed. Running the two together as shown assigns seed 0 to `equal_celltype` and seed 1 to `observed_abundance`, matching the current dataset YAMLs. Running the observed dataset by itself with `--seed 0` would produce a different simulation.

This command writes both modality H5ADs, feature lists, truth CSVs, source-cell JSONL files, dataset YAMLs, and registry entries.

## 11. Optional: recreate the deprecated CellTypist fallback

The canonical PBMC build does not use this branch. Recreate it only to mirror the current directory or preserve fallback provenance:

```bash
python -m pip install celltypist==1.7.1

PYTHONPATH=src .venv/bin/python scripts/prepare_pbmc_multiome_labels.py \
  --allow-celltypist-fallback \
  --overwrite
```

This downloads the `Immune_All_High.pkl` and `Immune_All_Low.pkl` CellTypist models and recreates the fallback mapping, summary, and manifest under:

```text
data/raw/sources/celltypist/pbmc_granulocyte_sorted_10k/
```

Do not feed that fallback mapping into `prepare_pbmc_multiome_references.py` unless intentionally changing the canonical PBMC labels.

## 12. Local-only archived results

Zenodo does not contain these six small files:

```text
data/archive/legacy_results/cell2location_results/
  deconvolution_plot.png
  means_cell_abundance_w_sf.csv
  q05_cell_abundance_w_sf.csv

data/archive/legacy_results/rctd_results/
  estimated_proportions.csv
  ground_truth.png
  rctd_vs_ground_truth.png
```

Together they occupy about 2.5 MiB and are not inputs to the maintained unified runner. For exact preservation, transfer those two directories from the original machine, for example:

```bash
rsync -av \
  original-machine:/absolute/path/to/deconvATAC/data/archive/legacy_results/cell2location_results/ \
  data/archive/legacy_results/cell2location_results/

rsync -av \
  original-machine:/absolute/path/to/deconvATAC/data/archive/legacy_results/rctd_results/ \
  data/archive/legacy_results/rctd_results/
```

They can instead be regenerated with the legacy scripts, but exact reproduction is not guaranteed: Cell2location is stochastic and the RCTD script requires `spacexr`, `rpy2`, and adjustment of its R-library path.

```bash
python -m pip install -e ".[cell2location]"
PYTHONPATH=src .venv/bin/python scripts/legacy/run_cell2location.py
PYTHONPATH=src .venv/bin/python scripts/legacy/run_rctd.py
PYTHONPATH=src .venv/bin/python scripts/migrate_legacy_results.py --execute
```

If desired, recreate the current archive manifest after transfer:

```bash
.venv/bin/python - <<'PY'
from pathlib import Path

path = Path("data/archive/legacy_results/manifest.tsv")
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(
    "source\tdestination\tstatus\n"
    "cell2location_results\tdata/archive/legacy_results/cell2location_results\tarchived\n"
    "rctd_results\tdata/archive/legacy_results/rctd_results\tarchived\n"
    "data/deconvolution_results\tdata/archive/legacy_results/deconvolution_results\tarchived\n"
)
PY
```

## 13. Validate the rebuilt tree

First inspect disk usage and the registry:

```bash
du -sh data
sed -n '1,200p' data/registry/datasets.yaml
```

The registry should contain these seven IDs:

```text
russell_250
human_cardiac_niches_sim_1zone_3ct_low_density
human_cardiac_niches_sim_1zone_10ct
human_cardiac_niches_sim_4zone_stripes
human_cardiac_niches_sim_4zone_circles
pbmc_granulocyte_sorted_10k_sim_equal_celltype
pbmc_granulocyte_sorted_10k_sim_observed_abundance
```

Validate the expected H5AD shapes without loading matrices into memory:

```bash
PYTHONPATH=src .venv/bin/python - <<'PY'
from pathlib import Path

import anndata as ad

expected = {
    "data/raw/references/human_cardiac_niches/atac/reference.h5ad": (139835, 429828),
    "data/raw/references/human_cardiac_niches/rna/reference.h5ad": (704296, 32732),
    "data/raw/references/russell_250/atac/reference.h5ad": (2535, 53451),
    "data/raw/references/pbmc_granulocyte_sorted_10k_multiome/atac/reference.h5ad": (9627, 143887),
    "data/raw/references/pbmc_granulocyte_sorted_10k_multiome/rna/reference.h5ad": (9627, 36601),
    "data/processed/datasets/russell_250/atac/spatial.h5ad": (360, 53451),
}

heart_ids = [
    "human_cardiac_niches_sim_1zone_3ct_low_density",
    "human_cardiac_niches_sim_1zone_10ct",
    "human_cardiac_niches_sim_4zone_stripes",
    "human_cardiac_niches_sim_4zone_circles",
]
for dataset_id in heart_ids:
    expected[f"data/processed/datasets/{dataset_id}/atac/spatial.h5ad"] = (961, 429828)
    expected[f"data/processed/datasets/{dataset_id}/rna/spatial.h5ad"] = (961, 32732)

pbmc_ids = [
    "pbmc_granulocyte_sorted_10k_sim_equal_celltype",
    "pbmc_granulocyte_sorted_10k_sim_observed_abundance",
]
for dataset_id in pbmc_ids:
    expected[f"data/processed/datasets/{dataset_id}/atac/spatial.h5ad"] = (1024, 143887)
    expected[f"data/processed/datasets/{dataset_id}/rna/spatial.h5ad"] = (1024, 36601)

for filename, wanted in expected.items():
    path = Path(filename)
    if not path.is_file():
        raise SystemExit(f"missing: {path}")
    obj = ad.read_h5ad(path, backed="r")
    observed = obj.shape
    obj.file.close()
    if observed != wanted:
        raise SystemExit(f"shape mismatch for {path}: {observed} != {wanted}")
    print(f"OK  {observed!s:20s}  {path}")
PY
```

Finally, perform the maintained loader's feature-alignment and truth validation for every registered modality. This is more expensive because it reads each selected 20,000-feature view:

```bash
PYTHONPATH=src MPLCONFIGDIR="$MPLCONFIGDIR" .venv/bin/python - <<'PY'
import gc
from pathlib import Path

import yaml

from deconvatac.data import load_deconvolution_input

registry_path = Path("data/registry/datasets.yaml")
registry = yaml.safe_load(registry_path.read_text()) or {}

for dataset_id, entry in registry.items():
    config_path = Path(entry["config"] if isinstance(entry, dict) else entry)
    config = yaml.safe_load(config_path.read_text())
    for modality in config["modalities"]:
        loaded = load_deconvolution_input(
            dataset_id=dataset_id,
            modality=modality,
            feature_set="highly_variable",
            registry_path=registry_path,
            project_root=Path.cwd(),
        )
        print(
            "OK",
            dataset_id,
            modality,
            loaded.reference.shape,
            loaded.spatial.shape,
        )
        del loaded
        gc.collect()
PY
```

After successful extraction and validation, the two staged Zenodo ZIP files are no longer needed. Keep them as an offline backup or delete those two specific files to recover about 1.74 GiB.

## Authoritative source links

- Heart Cell Atlas / CELLxGENE collection: <https://cellxgene.cziscience.com/collections/3116d060-0a8e-4767-99bb-e866badea1ed>
- Heart Cell Atlas v2 files: <https://cellgeni.cog.sanger.ac.uk/heartcellatlas/v2/>
- deconvATAC publication archive: <https://zenodo.org/records/15089738>
- 10x PBMC Multiome dataset: <https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0>
- SnapATAC2 PBMC helper documentation: <https://scverse.org/SnapATAC2/version/2.9/api/_autosummary/snapatac2.datasets.pbmc10k_multiome.html>
- GET Foundation PBMC preparation notebook: <https://github.com/GET-Foundation/get_model/blob/master/tutorials/prepare_pbmc.ipynb>
