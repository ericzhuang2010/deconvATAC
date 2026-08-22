import pandas as pd

from deconvatac.data import normalize_proportions
from deconvatac.methods.cell2location import clean_cell2location_columns
from deconvatac.methods.destvi import read_destvi_proportions
from deconvatac.methods.rctd import read_rctd_proportions
from deconvatac.methods.spatialdwls import read_spatialdwls_proportions
from deconvatac.methods.tangram import read_tangram_proportions


def test_normalize_proportions_handles_zero_rows():
    values = pd.DataFrame({"A": [2.0, 0.0], "B": [2.0, 0.0]}, index=["s1", "s2"])

    proportions = normalize_proportions(values)

    assert proportions.loc["s1", "A"] == 0.5
    assert proportions.loc["s1", "B"] == 0.5
    assert proportions.loc["s2"].sum() == 0.0


def test_cell2location_columns_are_cleaned():
    columns = clean_cell2location_columns(
        pd.Index(["meanscell_abundance_w_sf_T_cell", "q05_cell_abundance_w_sf_B_cell", "other"])
    )

    assert list(columns) == ["T_cell", "B_cell", "other"]


def test_rctd_reader_normalizes_csv(tmp_path):
    path = tmp_path / "estimated_proportions.csv"
    pd.DataFrame({"A": [2.0], "B": [2.0]}, index=["s1"]).to_csv(path)

    proportions = read_rctd_proportions(path)

    assert proportions.loc["s1", "A"] == 0.5
    assert proportions.loc["s1", "B"] == 0.5


def test_tangram_reader_normalizes_prediction_csv(tmp_path):
    path = tmp_path / "tangram_ct_pred.csv"
    pd.DataFrame({"A": [2.0], "B": [6.0]}, index=["s1"]).to_csv(path)

    proportions = read_tangram_proportions(path)

    assert proportions.loc["s1", "A"] == 0.25
    assert proportions.loc["s1", "B"] == 0.75


def test_destvi_reader_normalizes_prediction_csv(tmp_path):
    path = tmp_path / "predicted_proportions.csv"
    pd.DataFrame({"A": [0.2], "B": [0.2]}, index=["s1"]).to_csv(path)

    proportions = read_destvi_proportions(path)

    assert proportions.loc["s1", "A"] == 0.5
    assert proportions.loc["s1", "B"] == 0.5


def test_spatialdwls_reader_handles_wide_giotto_csv(tmp_path):
    path = tmp_path / "proportions.csv"
    pd.DataFrame(
        {
            "cell_ID": ["s1", "s2"],
            "A": [2.0, 0.0],
            "B": [6.0, 4.0],
        }
    ).to_csv(path)

    proportions = read_spatialdwls_proportions(path)

    assert proportions.loc["s1", "A"] == 0.25
    assert proportions.loc["s1", "B"] == 0.75
    assert proportions.loc["s2", "B"] == 1.0


def test_spatialdwls_reader_handles_long_giotto_csv(tmp_path):
    path = tmp_path / "proportions.csv"
    pd.DataFrame(
        {
            "cell_ID": ["s1", "s1", "s2", "s2"],
            "cell_type": ["A", "B", "A", "B"],
            "proportion": [2.0, 6.0, 0.0, 4.0],
        }
    ).to_csv(path)

    proportions = read_spatialdwls_proportions(path)

    assert proportions.loc["s1", "A"] == 0.25
    assert proportions.loc["s1", "B"] == 0.75
    assert proportions.loc["s2", "B"] == 1.0
