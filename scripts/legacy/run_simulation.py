# This file is not used

from deconvatac.tl import generate_spatial_data
import muon as mu
mdata = mu.read("data/raw/references/human_cardiac_niches/human_cardiac_niches.h5mu")
params = {
    "n_regions": 4,
    "cell_type_number": [10, 5, 7, 3],
    "cell_number_nu": [20, 20, 20, 20],
    "cell_number_mean": [15, 10, 12, 5],
    "region_type": "circles",
    "num_spots": 1000,
    "balance": "balanced",
}
simulation, samples = generate_spatial_data(mdata, cell_type_key="cell_type", **params)
simulation.write("my_simulated_spatial.h5mu")
