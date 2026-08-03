# How the Simulation Algorithm Works

## Key Terms

### Spot

A **spot** is one spatial observation.

In the data files, one spot is one row in:

```text
<modality>/spatial.h5ad
```

This is the unit that deconvolution methods operate on. For each spot, a method estimates the mixture of cell types.

### Spatial Region

A **spatial region** is not the same as a spot.

A spatial region is a larger simulated zone that contains many spots. The simulation uses regions to control which cell types appear in different parts of the synthetic tissue.

For example, a Heart simulation can use:

```yaml
n_regions: 4
region_type: circles
```

or:

```yaml
n_regions: 4
region_type: stripes
```

This means the simulation divides the spatial grid into 4 larger regions. Each individual spot belongs to one of those regions. The region then helps determine which cell types can be sampled for that spot.

In short:

```text
region = larger simulated spatial zone
spot = individual spatial location inside a region
```

## Heart Simulation Workflow

For each simulated spot, the algorithm:

1. Chooses the spot's spatial region.
2. Chooses which cell types are allowed in that region.
3. Draws the number of cells to include in the spot.
4. Randomly samples reference cells from the allowed cell types.
5. Sums the sampled source cells' ATAC and RNA profiles to create the spot.
6. Computes ground-truth cell-type proportions from those sampled source cells.
7. Writes the sampled cell IDs and cell types to:

```text
simulation/source_cells_by_spot.jsonl
```

## Ground Truth

For simulated data, ground truth is based on the simulation.

It is not experimentally measured truth. It is known because the simulation records exactly which reference cells were used to create each synthetic spot.

Example:

If a spot is created from 5 sampled source cells:

```text
3 Mural cells
2 Mesothelial cells
```

then the ground-truth proportions are:

```text
Mural cell: 3 / 5 = 0.6
Mesothelial cell: 2 / 5 = 0.4
all other cell types: 0
```

These values are stored in:

```text
truth/proportions.csv
```

## Simulation Provenance

The simulation provenance file records which source cells were used to create each synthetic spot.

Example row:

```json
{
  "spot_id": "0",
  "region": [0, 0, 0, 0, 0],
  "cell_id": ["cell_a", "cell_b", "cell_c"],
  "cell_type": ["Mural cell", "Mural cell", "Mesothelial cell"]
}
```

Field meanings:

- `spot_id`: the spatial spot ID.
- `region`: the simulated region assignment for the source cells in that spot.
- `cell_id`: reference cell IDs sampled to create the spot.
- `cell_type`: cell-type labels for the source cells.

## Why This Matters

The simulation creates three connected outputs:

```text
spatial.h5ad
truth/proportions.csv
simulation/source_cells_by_spot.jsonl
```

They are connected:

- `spatial.h5ad` contains the synthetic spatial profiles to deconvolve.
- `truth/proportions.csv` contains the known answer for benchmarking.
- `simulation/source_cells_by_spot.jsonl` explains how each spot was generated.

This lets us benchmark methods like Cell2location, RCTD, Tangram, DestVI, SpatialDWLS, and NNLS because we know the true cell-type proportions for every simulated spot.
