from __future__ import annotations

import pytest

from scripts.build_shapemix_gse246791_reference import (
    CELL_TYPES,
    MIN_CELLS_PER_TYPE,
    broad_label,
)


def test_gse246791_broad_ontology_is_frozen_in_output_order() -> None:
    assert CELL_TYPES == (
        "Excitatory neurons",
        "Inhibitory neurons",
        "Other neurons",
        "Astroglia/ependymal",
        "Oligodendrocytes",
        "OPCs",
        "Microglia/macrophages",
        "Other immune",
        "Vascular/stromal",
    )
    assert broad_label("Glut", "1 Example") == "Excitatory neurons"
    assert broad_label("Gaba", "2 Example") == "Inhibitory neurons"
    assert broad_label("Gly-Gaba", "3 Example") == "Inhibitory neurons"
    assert broad_label("Dopa", "4 Example") == "Other neurons"


def test_gse246791_support_gate_preserves_the_audited_rare_immune_class() -> None:
    assert MIN_CELLS_PER_TYPE["Other immune"] == 20
    assert {MIN_CELLS_PER_TYPE[cell_type] for cell_type in CELL_TYPES[:-2]} == {100}
    assert MIN_CELLS_PER_TYPE["Vascular/stromal"] == 100


@pytest.mark.parametrize(
    "identifier, expected",
    [
        (316, "Astroglia/ependymal"),
        (325, "Astroglia/ependymal"),
        (326, "OPCs"),
        (327, "Oligodendrocytes"),
        (328, "Oligodendrocytes"),
        (329, "Vascular/stromal"),
        (333, "Vascular/stromal"),
        (334, "Microglia/macrophages"),
        (335, "Microglia/macrophages"),
        (336, "Other immune"),
        (338, "Other immune"),
        (339, None),
    ],
)
def test_gse246791_non_neuronal_author_subclasses_are_frozen(
    identifier: int, expected: str | None
) -> None:
    assert broad_label("NN", f"{identifier} Author subclass") == expected


def test_gse246791_non_neuronal_mapping_fails_closed_outside_frozen_ranges() -> None:
    with pytest.raises(ValueError, match="unexpected subclass"):
        broad_label("NN", "315 Unexpected")
