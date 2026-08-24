# Russell Simulation Provenance Plan

## Purpose

Heart simulation datasets contain:

```text
simulation/source_cells_by_spot.jsonl
```

This file records which reference cells were used to create each synthetic spatial spot. Russell currently does not have this file. This document explains why Heart and Russell differ, and how to reconstruct an equivalent Russell simulation-provenance file for review.

## Current Difference

### Heart Simulations

Heart spatial spots were created by a random simulation algorithm.

For each spot, the simulation:

1. Chooses a spatial region.
2. Chooses which cell types are allowed in that region.
3. Draws a number of cells for the spot.
4. Randomly samples reference cells from allowed cell types.
5. Sums the sampled source cells' ATAC and RNA profiles to create the spot.
6. Computes ground-truth cell-type proportions from those sampled source cells.
7. Writes the source cell IDs and cell types to `simulation/source_cells_by_spot.jsonl`.

So for Heart, `simulation/source_cells_by_spot.jsonl` is direct provenance from the random simulation.

Example row:

```json
{
  "spot_id": "0",
  "region": [0, 0, 0, 0, 0],
  "cell_id": ["cell_a", "cell_b", "cell_c"],
  "cell_type": ["Mural cell", "Mural cell", "Mesothelial cell"]
}
```

### Russell

Russell was not created by the same random simulation algorithm.

The reproduction notebook shows that Russell spatial spots were created by coordinate binning:

1. Start from a single-cell ATAC reference with original spatial coordinates in `reference.obsm["spatial"]`.
2. Divide the coordinate plane into square windows.
3. For every non-empty window, collect all reference cells inside that window.
4. Sum those cells' ATAC profiles to create one synthetic spatial spot.
5. Store the window coordinate as the spot's spatial coordinate.
6. Compute ground-truth cell-type proportions from the grouped cells.

So for Russell, the equivalent simulation provenance is not "randomly sampled cells"; it is "cells grouped into each spatial bin."

## Existing Evidence

The relevant notebook is:

```text
notebooks/reproduction/simulation_artificial_visium.ipynb
```

The important function is `Simulated(adata, window, layer=None)`. Its logic is:

```python
spatial_loc = pd.DataFrame(adata.obsm["spatial"], index=adata.obs_names, columns=["x", "y"])

for x in np.arange(spatial_loc["x"].min() // window, spatial_loc["x"].max() // window + 1):
    for y in np.arange(spatial_loc["y"].min() // window, spatial_loc["y"].max() // window + 1):
        tmp_loc = spatial_loc[
            (x * window < spatial_loc["x"])
            & (spatial_loc["x"] < (x + 1) * window)
            & (y * window < spatial_loc["y"])
            & (spatial_loc["y"] < (y + 1) * window)
        ]
        if len(tmp_loc) > 0:
            combined_spot_loc.append([x, y])
            combined_spot.append(tmp_loc.index.to_list())
```

The notebook uses:

```python
window = 250
```

I verified that reconstructing Russell with this window size gives:

```text
reconstructed spots: 360
current spatial spots: 360
cell_count match: True
spatial coordinate match: True
truth proportion max absolute difference: 0.0
```

This means the current Russell reference and spatial files contain enough information to reconstruct the per-spot grouped-cell simulation provenance exactly.

## Proposed Output

Create:

```text
data/processed/datasets/russell_250/simulation/source_cells_by_spot.jsonl
```

This uses the same file name as Heart for schema consistency across generated datasets. The method is recorded separately so it is clear that Russell uses coordinate binning rather than random sampling.

Each row would contain:

```json
{
  "spot_id": "0",
  "grid_window": [0, 0],
  "cell_id": ["AAGGATTAGTTGTCTT-1"],
  "cell_type": ["tumour_1"]
}
```

Field meanings:

- `spot_id`: row ID from `data/processed/datasets/russell_250/atac/spatial.h5ad`.
- `grid_window`: integer `[x, y]` bin coordinate used to create the synthetic spot.
- `cell_id`: reference cell IDs grouped into that bin.
- `cell_type`: cell-type labels for those reference cells.

I recommend `grid_window` instead of `region` because Russell bins are geometric coordinate windows, not simulated biological regions.

## Proposed Algorithm

Inputs:

```text
data/processed/references/russell_250/atac/reference.h5ad
data/processed/datasets/russell_250/atac/spatial.h5ad
```

Parameters:

```text
window = 250
labels_key = "cell_type"
reference_spatial_key = "spatial"
spatial_key = "spatial"
```

Algorithm:

1. Load Russell reference ATAC:

```python
reference = ad.read_h5ad("data/processed/references/russell_250/atac/reference.h5ad")
```

2. Load Russell spatial ATAC:

```python
spatial = ad.read_h5ad("data/processed/datasets/russell_250/atac/spatial.h5ad")
```

3. Build a coordinate table from reference cells:

```python
reference_coords = pd.DataFrame(
    reference.obsm["spatial"],
    index=reference.obs_names,
    columns=["x", "y"],
)
```

4. Recreate bins using the notebook logic:

```python
rows = []
spot_idx = 0

for x in np.arange(reference_coords["x"].min() // window, reference_coords["x"].max() // window + 1):
    for y in np.arange(reference_coords["y"].min() // window, reference_coords["y"].max() // window + 1):
        in_window = reference_coords[
            (x * window < reference_coords["x"])
            & (reference_coords["x"] < (x + 1) * window)
            & (y * window < reference_coords["y"])
            & (reference_coords["y"] < (y + 1) * window)
        ]

        if len(in_window) == 0:
            continue

        cell_ids = in_window.index.astype(str).tolist()
        rows.append(
            {
                "spot_id": str(spatial.obs_names[spot_idx]),
                "grid_window": [int(x), int(y)],
                "cell_id": cell_ids,
                "cell_type": reference.obs.loc[cell_ids, labels_key].astype(str).tolist(),
            }
        )
        spot_idx += 1
```

5. Write JSONL:

```python
simulation_dir = Path("data/processed/datasets/russell_250/simulation")
simulation_dir.mkdir(parents=True, exist_ok=True)

with (simulation_dir / "source_cells_by_spot.jsonl").open("w") as handle:
    for row in rows:
        handle.write(json.dumps(row) + "\n")
```

6. Update `dataset.yaml` with:

```yaml
simulation:
  source_cells_by_spot: data/processed/datasets/russell_250/simulation/source_cells_by_spot.jsonl
  source_cells_method: coordinate_binning
  coordinate_window: 250
```

## Validation Plan

After writing the file, validate:

1. The number of JSONL rows equals `spatial.n_obs`.

2. Per-row number of `cell_id` entries equals `spatial.obs["cell_count"]`.

3. Recomputed cell-type proportions from `cell_type` equal:

```text
data/processed/datasets/russell_250/truth/proportions.csv
```

4. Reconstructed `grid_window` coordinates equal:

```python
spatial.obsm["spatial"]
```

5. All `cell_id` values exist in:

```python
reference.obs_names
```

## Naming Decision

Use this canonical path for generated spot provenance:

```text
simulation/source_cells_by_spot.jsonl
```

Use this dataset config key:

```yaml
simulation:
  source_cells_by_spot: data/processed/datasets/<dataset_id>/simulation/source_cells_by_spot.jsonl
```

This avoids the generic `metadata/` name and avoids saying "sampled" for Russell, where the source cells were grouped by coordinate windows rather than randomly sampled.

## Caveats

- This reconstruction depends on the notebook's `window = 250` convention.
- It assumes the current Russell `spatial.h5ad` was generated from the current Russell reference using the notebook binning logic.
- It does not create new biological information; it reconstructs provenance that is already implied by the reference coordinates, spatial coordinates, cell counts, and truth proportions.
