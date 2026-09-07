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

report_is_complete() {
    local path="$1"
    [[ -s "$path" ]] && awk '$1 == "status:" && ($2 == "complete" || $2 == "passed") {ok=1} END {exit !ok}' "$path"
}

run_guarded() {
    local label="$1"
    shift
    local status=0
    while :; do
        log "full_evaluation_v2_stage_start label=$label"
        set +e
        scripts/run_shapemix_low_impact.sh "$@"
        status=$?
        set -e
        if [[ "$status" -eq 0 ]]; then
            log "full_evaluation_v2_stage_complete label=$label"
            return 0
        fi
        if [[ "$status" -ne 75 ]]; then
            log "full_evaluation_v2_stage_failed label=$label status=$status"
            return "$status"
        fi
        log "full_evaluation_v2_resource_wait label=$label"
        sleep 60
    done
}

readonly PRIMARY_SUMMARY="results/primary/shapemix_primary_ablation_protocol_v1_summary/summary.yaml"
readonly CUDA_REPORT="results/development/shapemix_gpu_qualification_v2/qualification_report.yaml"
readonly GSE129_SUMMARY="results/external_validation/shapemix_gse129785_v2/shapemix_gse129785_external_protocol_v2_cpu_log_abundance/evidence_summary.yaml"
readonly GSE194_SUMMARY="results/external_validation/shapemix_gse194122_lodo_v1/shapemix_gse194122_broad7_lodo_protocol_v1_cuda/evidence_summary.yaml"
readonly PBMC_SUMMARY="results/sensitivity/shapemix_pbmc_stress_v2/shapemix_pbmc_stress_protocol_v2_cuda/evidence_summary.yaml"
readonly EMBRYO_REFERENCE="data/processed/references/gse216371_mouse_embryo_e13_major_types_v1/reference.yaml"
readonly GSE205055_SUMMARY="results/real_spatial/shapemix_gse205055_v1/shapemix_gse205055_real_spatial_protocol_v1_cuda/evidence_summary.yaml"
readonly GSE263333_SUMMARY="results/real_spatial/shapemix_gse263333_v1/shapemix_gse263333_real_spatial_protocol_v1_cuda/evidence_summary.yaml"
readonly FULL_SUMMARY="results/external_validation/shapemix_full_evaluation_v2/evidence_summary.yaml"

require_file "$PRIMARY_SUMMARY" primary_summary
require_file "$CUDA_REPORT" cuda_qualification
require_file "$GSE129_SUMMARY" gse129785_summary
require_file "$GSE194_SUMMARY" gse194122_summary

if ! report_is_complete "$PBMC_SUMMARY"; then
    run_guarded pbmc_stress_v2_tests .venv/bin/python -m pytest -q \
        tests/test_prepare_shapemix_pbmc_sensitivity.py \
        tests/test_prepare_shapemix_pbmc_sensitivity_v2.py \
        tests/test_summarize_shapemix_pbmc_sensitivity.py
    run_guarded pbmc_stress_v2_materialize .venv/bin/python \
        scripts/prepare_shapemix_pbmc_sensitivity_v2.py all
    run_guarded pbmc_stress_v2_layout .venv/bin/python \
        scripts/validate_shapemix_file_layout.py \
        --experiment-config configs/experiments/shapemix_pbmc_stress_v2.yaml \
        --allow-existing-results
    run_guarded pbmc_stress_v2_run .venv/bin/python scripts/run_deconvolution.py \
        --experiment-config configs/experiments/shapemix_pbmc_stress_v2.yaml \
        --resume
    run_guarded pbmc_stress_v2_summary .venv/bin/python \
        scripts/summarize_shapemix_pbmc_sensitivity_v2.py
fi
require_file "$PBMC_SUMMARY" pbmc_stress_v2_summary

if [[ ! -s "$EMBRYO_REFERENCE" ]]; then
    run_guarded embryo_reference_tests .venv/bin/python -m pytest -q \
        tests/test_shapemix_gse216371_stream.py \
        tests/test_materialize_shapemix_gse216371_reference.py
    run_guarded embryo_reference_build .venv/bin/python \
        scripts/materialize_shapemix_gse216371_reference.py --stage all
fi
require_file "$EMBRYO_REFERENCE" embryo_reference

if ! report_is_complete "$GSE205055_SUMMARY" || ! report_is_complete "$GSE263333_SUMMARY"; then
    run_guarded real_spatial_tests .venv/bin/python -m pytest -q \
        tests/test_materialize_shapemix_real_spatial.py \
        tests/test_prepare_shapemix_reference_marker_features.py \
        tests/test_shapemix_spatial_preprocessing.py \
        tests/test_summarize_shapemix_real_spatial.py
    run_guarded real_spatial_marker_features .venv/bin/python \
        scripts/prepare_shapemix_reference_marker_features.py
    run_guarded real_spatial_materialize .venv/bin/python \
        scripts/materialize_shapemix_real_spatial.py
    run_guarded real_spatial_layout .venv/bin/python \
        scripts/validate_shapemix_file_layout.py \
        --experiment-config configs/experiments/shapemix_gse205055_real_spatial_v1.yaml \
        --experiment-config configs/experiments/shapemix_gse263333_real_spatial_v1.yaml \
        --allow-existing-results
    run_guarded gse205055_real_spatial_run .venv/bin/python scripts/run_deconvolution.py \
        --experiment-config configs/experiments/shapemix_gse205055_real_spatial_v1.yaml \
        --resume
    run_guarded gse263333_real_spatial_run .venv/bin/python scripts/run_deconvolution.py \
        --experiment-config configs/experiments/shapemix_gse263333_real_spatial_v1.yaml \
        --resume
    run_guarded real_spatial_summary .venv/bin/python scripts/summarize_shapemix_real_spatial.py \
        --experiment-config configs/experiments/shapemix_gse205055_real_spatial_v1.yaml \
        --batch-dir results/real_spatial/shapemix_gse205055_v1/shapemix_gse205055_real_spatial_protocol_v1_cuda \
        --experiment-config configs/experiments/shapemix_gse263333_real_spatial_v1.yaml \
        --batch-dir results/real_spatial/shapemix_gse263333_v1/shapemix_gse263333_real_spatial_protocol_v1_cuda
fi
require_file "$GSE205055_SUMMARY" gse205055_real_spatial_summary
require_file "$GSE263333_SUMMARY" gse263333_real_spatial_summary

if [[ ! -s "$FULL_SUMMARY" ]]; then
    run_guarded full_evaluation_v2_synthesis .venv/bin/python \
        scripts/summarize_shapemix_full_evaluation.py \
        --config configs/experiments/shapemix_full_evaluation_v2.yaml
fi
require_file "$FULL_SUMMARY" full_evaluation_v2_summary
log full_evaluation_v2_complete
