from pathlib import Path

from scripts.download_gse194122 import (
    load_config,
    pilot_sample,
    resources_from_config,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/data_sources/shapemix_gse194122.yaml"


def test_frozen_atac_inventory_is_complete_and_injective():
    config = load_config(CONFIG)
    samples = config["atac_samples"]

    assert len(samples) == 13
    for key in ("gsm", "biosample", "experiment", "run", "bam_url", "bam_md5"):
        assert len({sample[key] for sample in samples}) == 13
    assert {(sample["site"], sample["donor"]) for sample in samples} == {
        (1, 1), (1, 2), (1, 3),
        (2, 1), (2, 4), (2, 5),
        (3, 10), (3, 3), (3, 6), (3, 7),
        (4, 1), (4, 8), (4, 9),
    }
    assert sum(sample["bam_bytes"] for sample in samples) == 531_234_571_942


def test_default_scope_does_not_cross_the_pilot_gate():
    config = load_config(CONFIG)
    core = resources_from_config(config, include_pilot_bam=False)
    gated = resources_from_config(config, include_pilot_bam=True)

    assert [resource.payload for resource in core] == ["gzip", "gzip"]
    assert len(gated) == 3
    assert gated[-1].gsm == config["pilot_gsm"] == "GSM5828489"
    assert gated[-1].expected_md5 == pilot_sample(config)["bam_md5"]
    assert gated[-1].expected_bytes == 16_961_779_971


def test_author_fragment_scope_contains_all_13_complete_trios():
    config = load_config(CONFIG)
    resources = resources_from_config(
        config, include_pilot_bam=False, include_fragments=True
    )
    fragments = resources[2:]

    assert len(fragments) == 39
    assert {resource.sample_key for resource in fragments} == {
        f"s{site}d{donor}"
        for site, donor in {
            (1, 1), (1, 2), (1, 3),
            (2, 1), (2, 4), (2, 5),
            (3, 10), (3, 3), (3, 6), (3, 7),
            (4, 1), (4, 8), (4, 9),
        }
    }
    assert sum(resource.expected_bytes for resource in fragments) == 31_687_758_209
    assert all(resource.expected_etag for resource in fragments)
