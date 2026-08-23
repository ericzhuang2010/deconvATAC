import os
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from deconvatac.data import DeconvolutionInput, DeconvolutionResult
from deconvatac.methods import get_method, list_methods


def _toy_input():
    reference = ad.AnnData(
        X=np.array([[5, 0], [4, 1], [0, 5]], dtype=float),
        obs=pd.DataFrame({"cell_type": ["A", "A", "B"]}, index=["c1", "c2", "c3"]),
        var=pd.DataFrame(index=["f1", "f2"]),
    )
    spatial = ad.AnnData(
        X=np.array([[9, 1], [1, 8]], dtype=float),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=["f1", "f2"]),
    )
    spatial.obsm["spatial"] = np.array([[0, 0], [1, 0]], dtype=float)
    return DeconvolutionInput(
        dataset_id="toy",
        modality="atac",
        feature_set="all",
        spatial=spatial,
        reference=reference,
        labels_key="cell_type",
    )


def test_builtin_methods_are_registered():
    methods = list_methods()
    assert "cell2location" in methods
    assert "destvi" in methods
    assert "nnls" in methods
    assert "rctd" in methods
    assert "spatialdwls" in methods
    assert "tangram" in methods
    assert "shapemix" in methods
    assert get_method("nnls").method_name == "nnls"
    assert get_method("shapemix").method_name == "shapemix"


def test_listing_and_running_nnls_do_not_import_shapemix_dependencies():
    root = Path(__file__).resolve().parents[1]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(root / "src")
    command = """
import sys
from deconvatac.methods import get_method, list_methods
assert not any(name.startswith('deconvatac.shapemix') for name in sys.modules)
assert 'pysam' not in sys.modules
assert 'shapemix' in list_methods()
assert 'deconvatac.methods.shapemix' not in sys.modules
get_method('nnls')
assert not any(name.startswith('deconvatac.shapemix') for name in sys.modules)
assert 'pysam' not in sys.modules
assert 'deconvatac.methods.shapemix' not in sys.modules
"""
    completed = subprocess.run(
        [sys.executable, "-c", command],
        cwd=root,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr


def test_nnls_baseline_remains_available():
    result = get_method("nnls")().run(_toy_input())

    assert isinstance(result, DeconvolutionResult)
    assert result.method == "nnls"
    assert list(result.proportions.columns) == ["A", "B"]
    assert np.allclose(result.proportions.sum(axis=1).values, 1.0)
