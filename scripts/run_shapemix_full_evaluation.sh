#!/usr/bin/env bash
set -euo pipefail

readonly ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

log() {
    printf '%s utc=%s\n' "$*" "$(date -u +%FT%TZ)"
}

require_file() {
    local path="$1"
    local label="$2"
    if [[ ! -s "$path" ]]; then
        printf 'Missing %s: %s\n' "$label" "$path" >&2
        return 1
    fi
}

wait_for_file() {
    local path="$1"
    local label="$2"
    while [[ ! -s "$path" ]]; do
        log "full_evaluation_waiting label=$label"
        sleep 60
    done
}

run_guarded() {
    local label="$1"
    shift
    local attempt=0
    local status=0
    while :; do
        attempt=$((attempt + 1))
        log "full_evaluation_stage_start label=$label attempt=$attempt"
        set +e
        scripts/run_shapemix_low_impact.sh "$@"
        status=$?
        set -e
        if [[ "$status" -eq 0 ]]; then
            log "full_evaluation_stage_complete label=$label"
            return 0
        fi
        if [[ "$status" -ne 75 ]]; then
            log "full_evaluation_stage_failed label=$label status=$status"
            return "$status"
        fi
        log "full_evaluation_resource_wait label=$label"
        sleep 60
    done
}

run_guarded_until_file() {
    local path="$1"
    local label="$2"
    shift 2
    local attempt=0
    local failures=0
    local status=0
    while [[ ! -s "$path" ]]; do
        attempt=$((attempt + 1))
        log "full_evaluation_stage_start label=$label attempt=$attempt"
        set +e
        scripts/run_shapemix_low_impact.sh "$@"
        status=$?
        set -e
        if [[ "$status" -eq 0 ]]; then
            require_file "$path" "$label"
            log "full_evaluation_stage_complete label=$label"
            return 0
        fi
        if [[ "$status" -eq 75 ]]; then
            log "full_evaluation_resource_wait label=$label"
            sleep 60
            continue
        fi
        failures=$((failures + 1))
        log "full_evaluation_stage_retry label=$label status=$status failures=$failures"
        if (( failures >= 10 )); then
            log "full_evaluation_stage_failed label=$label status=$status"
            return "$status"
        fi
        sleep 300
    done
    log "full_evaluation_stage_reused label=$label"
}

readonly ADULT_GATE="data/processed/shapemix/gse246791_mouse_brain_reference/manifests/adult_reference_preacquisition_gate.yaml"
readonly ADULT_READ_AUDIT="data/processed/shapemix/gse246791_mouse_brain_reference/source_audit/fragment_read_downloads.yaml"
readonly ADULT_REFERENCE="data/processed/references/gse246791_adult_mouse_brain_broad9_v1/reference.yaml"
readonly CUDA_REPORT="results/development/shapemix_gpu_qualification_v2/qualification_report.yaml"
readonly HUMAN_REFERENCE="data/processed/references/gse244618_human_hippocampus_donor3_region3_v1/reference.yaml"
readonly GSE194_SUMMARY="results/external_validation/shapemix_gse194122_lodo_v1/shapemix_gse194122_broad7_lodo_protocol_v1_cuda/evidence_summary.yaml"
readonly PBMC_SUMMARY="results/sensitivity/shapemix_pbmc_stress_v1/shapemix_pbmc_stress_protocol_v1_cuda/evidence_summary.yaml"
readonly EMBRYO_REFERENCE="data/processed/references/gse216371_mouse_embryo_e13_major_types_v1/reference.yaml"
readonly GSE205055_SUMMARY="results/real_spatial/shapemix_gse205055_v1/shapemix_gse205055_real_spatial_protocol_v1_cuda/evidence_summary.yaml"
readonly GSE263333_SUMMARY="results/real_spatial/shapemix_gse263333_v1/shapemix_gse263333_real_spatial_protocol_v1_cuda/evidence_summary.yaml"

wait_for_file "$ADULT_GATE" adult_preacquisition_gate
run_guarded_until_file "$ADULT_READ_AUDIT" adult_fragment_read_acquisition \
    .venv/bin/python scripts/download_shapemix_spatial.py \
    --config configs/data_sources/shapemix_gse246791_fragment_reads.yaml \
    --workers 4 \
    --timeout 1800
run_guarded adult_alignment_tests .venv/bin/python -m pytest -q tests/test_align_shapemix_gse246791.py tests/test_construct_shapemix_gse246791_fragments.py
run_guarded adult_alignment_compile_check .venv/bin/python -m py_compile scripts/align_shapemix_gse246791.py scripts/sort_shapemix_bam_stream.py

adult_samples=(
    GSM7877011
    GSM7877041
    GSM7876942
    GSM7877006
    GSM7877102
    GSM7877014
    GSM7877013
    GSM7877084
    GSM7877017
    GSM7876953
    GSM7876902
    GSM7877047
)
for gsm in "${adult_samples[@]}"; do
    fragment_manifest="data/processed/shapemix/gse246791_mouse_brain_reference/normalized_fragments/$gsm/manifest.yaml"
    if [[ ! -s "$fragment_manifest" ]]; then
        run_guarded "adult_align_$gsm" taskset --cpu-list 6,7 .venv-shapemix-fragments/bin/python scripts/align_shapemix_gse246791.py --gsm "$gsm"
        run_guarded "adult_fragments_$gsm" .venv-shapemix-fragments/bin/python scripts/construct_shapemix_gse246791_fragments.py --gsm "$gsm" --cleanup-bam
    else
        log "full_evaluation_stage_reused label=adult_fragments_$gsm"
    fi
    run_guarded "adult_shape_cache_$gsm" .venv/bin/python scripts/materialize_shapemix_gse246791_reference.py --stage fragment-cache --gsm "$gsm"
done
run_guarded adult_reference_assembly .venv/bin/python scripts/materialize_shapemix_gse246791_reference.py --stage reference
require_file "$ADULT_REFERENCE" adult_reference

if [[ ! -s "$CUDA_REPORT" ]]; then
    run_guarded cuda_v2_layout .venv/bin/python scripts/validate_shapemix_file_layout.py --experiment-config configs/experiments/shapemix_cuda_full_qualification_v2.yaml --allow-existing-results
    run_guarded cuda_v2_full .venv/bin/python scripts/run_deconvolution.py --experiment-config configs/experiments/shapemix_cuda_full_qualification_v2.yaml --overwrite
    run_guarded cuda_v2_summary .venv/bin/python scripts/summarize_shapemix_cuda_qualification.py --root results/development/shapemix_gpu_qualification_v2
fi
require_file "$CUDA_REPORT" cuda_qualification_report

if [[ ! -s "$HUMAN_REFERENCE" ]]; then
    run_guarded human_reference_tests .venv/bin/python -m pytest -q tests/test_prepare_shapemix_gse244618.py
    run_guarded human_reference_build .venv/bin/python scripts/prepare_shapemix_gse244618.py --stage all
fi
require_file "$HUMAN_REFERENCE" human_reference

if [[ ! -s "$GSE194_SUMMARY" ]]; then
    run_guarded gse194122_tests .venv/bin/python -m pytest -q tests/test_prepare_shapemix_gse194122.py tests/test_summarize_shapemix_gse194122.py
    run_guarded gse194122_materialize .venv/bin/python scripts/prepare_shapemix_gse194122.py all
    run_guarded gse194122_layout .venv/bin/python scripts/validate_shapemix_file_layout.py --experiment-config configs/experiments/shapemix_gse194122_lodo.yaml --allow-existing-results
    run_guarded gse194122_run .venv/bin/python scripts/run_deconvolution.py --experiment-config configs/experiments/shapemix_gse194122_lodo.yaml --resume
    run_guarded gse194122_summary .venv/bin/python scripts/summarize_shapemix_gse194122.py
fi
require_file "$GSE194_SUMMARY" gse194122_summary

if [[ ! -s "$PBMC_SUMMARY" ]]; then
    run_guarded pbmc_stress_tests .venv/bin/python -m pytest -q tests/test_prepare_shapemix_pbmc_sensitivity.py tests/test_summarize_shapemix_pbmc_sensitivity.py
    run_guarded pbmc_stress_materialize .venv/bin/python scripts/prepare_shapemix_pbmc_sensitivity.py all
    run_guarded pbmc_stress_layout .venv/bin/python scripts/validate_shapemix_file_layout.py --experiment-config configs/experiments/shapemix_pbmc_stress_v1.yaml --allow-existing-results
    run_guarded pbmc_stress_run .venv/bin/python scripts/run_deconvolution.py --experiment-config configs/experiments/shapemix_pbmc_stress_v1.yaml --resume
    run_guarded pbmc_stress_summary .venv/bin/python scripts/summarize_shapemix_pbmc_sensitivity.py
fi
require_file "$PBMC_SUMMARY" pbmc_stress_summary

if [[ ! -s "$EMBRYO_REFERENCE" ]]; then
    run_guarded embryo_reference_tests .venv/bin/python -m pytest -q tests/test_shapemix_gse216371_stream.py tests/test_materialize_shapemix_gse216371_reference.py
    run_guarded embryo_reference_build .venv/bin/python scripts/materialize_shapemix_gse216371_reference.py --stage all
fi
require_file "$EMBRYO_REFERENCE" embryo_reference

if [[ ! -s "$GSE205055_SUMMARY" || ! -s "$GSE263333_SUMMARY" ]]; then
    run_guarded real_spatial_tests .venv/bin/python -m pytest -q tests/test_shapemix_real_spatial_validation.py tests/test_shapemix_spatial_preprocessing.py
    run_guarded real_spatial_marker_features .venv/bin/python scripts/prepare_shapemix_reference_marker_features.py
    run_guarded real_spatial_materialize .venv/bin/python scripts/materialize_shapemix_real_spatial.py
    run_guarded real_spatial_layout .venv/bin/python scripts/validate_shapemix_file_layout.py --experiment-config configs/experiments/shapemix_gse205055_real_spatial_v1.yaml --experiment-config configs/experiments/shapemix_gse263333_real_spatial_v1.yaml --allow-existing-results
    run_guarded gse205055_real_spatial_run .venv/bin/python scripts/run_deconvolution.py --experiment-config configs/experiments/shapemix_gse205055_real_spatial_v1.yaml --resume
    run_guarded gse263333_real_spatial_run .venv/bin/python scripts/run_deconvolution.py --experiment-config configs/experiments/shapemix_gse263333_real_spatial_v1.yaml --resume
    run_guarded real_spatial_summary .venv/bin/python scripts/summarize_shapemix_real_spatial.py --experiment-config configs/experiments/shapemix_gse205055_real_spatial_v1.yaml --batch-dir results/real_spatial/shapemix_gse205055_v1/shapemix_gse205055_real_spatial_protocol_v1_cuda --experiment-config configs/experiments/shapemix_gse263333_real_spatial_v1.yaml --batch-dir results/real_spatial/shapemix_gse263333_v1/shapemix_gse263333_real_spatial_protocol_v1_cuda
fi
require_file "$GSE205055_SUMMARY" gse205055_real_spatial_summary
require_file "$GSE263333_SUMMARY" gse263333_real_spatial_summary

run_guarded full_evaluation_synthesis .venv/bin/python scripts/summarize_shapemix_full_evaluation.py
log full_evaluation_complete
