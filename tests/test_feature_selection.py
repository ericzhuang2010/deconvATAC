import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from deconvatac.pp import highly_accessible_peaks, highly_variable_peaks


def test_highly_variable_peaks_adds_boolean_var_column():
    adata = ad.AnnData(
        X=sparse.csr_matrix(
            np.array(
                [
                    [10, 0, 1],
                    [8, 0, 1],
                    [0, 9, 1],
                    [0, 7, 1],
                ],
                dtype=float,
            )
        ),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B", "B"]}),
        var=pd.DataFrame(index=["peak1", "peak2", "peak3"]),
    )

    highly_variable_peaks(adata, cluster_key="cell_type", n_top_features=2)

    assert "highly_variable" in adata.var
    assert adata.var["highly_variable"].dtype == bool
    assert int(adata.var["highly_variable"].sum()) == 2


def test_highly_accessible_peaks_adds_boolean_var_column():
    adata = ad.AnnData(
        X=sparse.csr_matrix(
            np.array(
                [
                    [1, 0, 1],
                    [1, 0, 0],
                    [1, 1, 0],
                ],
                dtype=float,
            )
        ),
        var=pd.DataFrame(index=["peak1", "peak2", "peak3"]),
    )

    highly_accessible_peaks(adata, n_top_features=1)

    assert "highly_accessible" in adata.var
    assert adata.var["highly_accessible"].dtype == bool
    assert list(adata.var_names[adata.var["highly_accessible"]]) == ["peak1"]
