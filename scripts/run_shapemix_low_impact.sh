#!/usr/bin/env bash
set -euo pipefail

readonly MAX_ONE_MINUTE_LOAD="6.0"
readonly MAX_WORKERS=2
readonly MIN_AVAILABLE_MEMORY_KIB=4194304
readonly MAX_GPU_MEMORY_USED_MIB=2048
readonly MAX_GPU_TEMPERATURE_C=79
readonly MAX_ALLOWED_DISPLAY_PROCESS_MEMORY_MIB=512
readonly STABILITY_SECONDS=60
readonly LOCK_PATH="/tmp/deconvatac-shapemix-resource.lock"

usage() {
    echo "Usage: $0 <command> [args ...]" >&2
}

unrelated_gpu_processes() {
    local raw_processes="$1"
    local pid process_name used_memory_mib executable
    while IFS=',' read -r pid process_name used_memory_mib; do
        pid="${pid//[[:space:]]/}"
        process_name="${process_name#${process_name%%[![:space:]]*}}"
        process_name="${process_name%${process_name##*[![:space:]]}}"
        used_memory_mib="${used_memory_mib//[[:space:]]/}"
        [[ -z "$pid" ]] && continue
        executable="${process_name##*/}"
        if [[ "$executable" == "gnome-remote-desktop-daemon" \
            && "$used_memory_mib" =~ ^[0-9]+$ \
            && "$used_memory_mib" -le "$MAX_ALLOWED_DISPLAY_PROCESS_MEMORY_MIB" ]]; then
            continue
        fi
        printf '%s, %s, %s\n' "$pid" "$process_name" "$used_memory_mib"
    done <<< "$raw_processes"
}

if (( $# == 0 )); then
    usage
    exit 64
fi

arguments=("$@")
for ((index = 0; index < ${#arguments[@]}; index++)); do
    argument="${arguments[index]}"
    worker_value=""
    if [[ "$argument" == "--workers" ]]; then
        if (( index + 1 >= ${#arguments[@]} )); then
            echo "--workers requires an integer value." >&2
            exit 64
        fi
        worker_value="${arguments[index + 1]}"
    elif [[ "$argument" == --workers=* ]]; then
        worker_value="${argument#--workers=}"
    fi
    if [[ -n "$worker_value" ]]; then
        if [[ ! "$worker_value" =~ ^[0-9]+$ ]] || (( worker_value < 1 || worker_value > MAX_WORKERS )); then
            echo "Refusing --workers=$worker_value; the active co-tenant limit is 1-$MAX_WORKERS." >&2
            exit 64
        fi
    fi
done

exec 9>"$LOCK_PATH"
if ! flock -n 9; then
    echo "Another deconvATAC task holds $LOCK_PATH; leaving this job queued." >&2
    exit 75
fi

one_minute_load="$(awk '{print $1}' /proc/loadavg)"
if ! awk -v observed="$one_minute_load" -v maximum="$MAX_ONE_MINUTE_LOAD" \
    'BEGIN { exit !(observed < maximum) }'; then
    echo "Load gate closed: one-minute load $one_minute_load is >= $MAX_ONE_MINUTE_LOAD." >&2
    exit 75
fi

if ! gpu_processes="$(nvidia-smi --query-compute-apps=pid,process_name,used_memory \
    --format=csv,noheader,nounits 2>/dev/null)"; then
    echo "Unable to query GPU process state; refusing a fail-open launch." >&2
    exit 69
fi
unrelated_processes="$(unrelated_gpu_processes "$gpu_processes")"
if [[ -n "${unrelated_processes//[$'\t\r\n ']/}" ]]; then
    echo "GPU gate closed; unrelated compute process(es) are active:" >&2
    echo "$unrelated_processes" >&2
    exit 75
fi

if ! gpu_state="$(nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu,temperature.gpu \
    --format=csv,noheader,nounits 2>/dev/null)"; then
    echo "Unable to query GPU capacity; refusing a fail-open launch." >&2
    exit 69
fi
IFS=',' read -r gpu_used_mib gpu_total_mib gpu_utilization gpu_temperature <<< "$gpu_state"
gpu_used_mib="${gpu_used_mib//[[:space:]]/}"
gpu_total_mib="${gpu_total_mib//[[:space:]]/}"
gpu_utilization="${gpu_utilization//[[:space:]]/}"
gpu_temperature="${gpu_temperature//[[:space:]]/}"
if (( gpu_used_mib > MAX_GPU_MEMORY_USED_MIB || gpu_temperature > MAX_GPU_TEMPERATURE_C )); then
    echo "GPU capacity gate closed: used=${gpu_used_mib}MiB temperature=${gpu_temperature}C." >&2
    exit 75
fi

available_kib="$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)"
if (( available_kib < MIN_AVAILABLE_MEMORY_KIB )); then
    echo "Memory gate closed: available=${available_kib}KiB is below ${MIN_AVAILABLE_MEMORY_KIB}KiB." >&2
    exit 75
fi

# Snakemake can leave a short, misleadingly idle gap between large jobs. Hold
# the project lock and require the host to remain safe across that gap before
# launching any ShapeMix command.
sleep "$STABILITY_SECONDS"
stable_one_minute_load="$(awk '{print $1}' /proc/loadavg)"
stable_available_kib="$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)"
if ! stable_gpu_processes="$(nvidia-smi --query-compute-apps=pid,process_name,used_memory \
    --format=csv,noheader,nounits 2>/dev/null)"; then
    echo "Unable to repeat the GPU process query; refusing a fail-open launch." >&2
    exit 69
fi
stable_unrelated_processes="$(unrelated_gpu_processes "$stable_gpu_processes")"
if ! stable_gpu_state="$(nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu,temperature.gpu \
    --format=csv,noheader,nounits 2>/dev/null)"; then
    echo "Unable to repeat the GPU capacity query; refusing a fail-open launch." >&2
    exit 69
fi
IFS=',' read -r stable_gpu_used_mib stable_gpu_total_mib stable_gpu_utilization stable_gpu_temperature <<< "$stable_gpu_state"
stable_gpu_used_mib="${stable_gpu_used_mib//[[:space:]]/}"
stable_gpu_total_mib="${stable_gpu_total_mib//[[:space:]]/}"
stable_gpu_utilization="${stable_gpu_utilization//[[:space:]]/}"
stable_gpu_temperature="${stable_gpu_temperature//[[:space:]]/}"
if ! awk -v observed="$stable_one_minute_load" -v maximum="$MAX_ONE_MINUTE_LOAD" \
    'BEGIN { exit !(observed < maximum) }' \
    || (( stable_available_kib < MIN_AVAILABLE_MEMORY_KIB )) \
    || [[ -n "${stable_unrelated_processes//[$'\t\r\n ']/}" ]] \
    || (( stable_gpu_used_mib > MAX_GPU_MEMORY_USED_MIB \
        || stable_gpu_temperature > MAX_GPU_TEMPERATURE_C )); then
    echo "Stability gate closed after ${STABILITY_SECONDS}s: load=${stable_one_minute_load} available_memory_kib=${stable_available_kib} gpu_used_mib=${stable_gpu_used_mib} gpu_temperature_c=${stable_gpu_temperature}." >&2
    exit 75
fi
one_minute_load="$stable_one_minute_load"
available_kib="$stable_available_kib"
gpu_used_mib="$stable_gpu_used_mib"
gpu_total_mib="$stable_gpu_total_mib"
gpu_utilization="$stable_gpu_utilization"
gpu_temperature="$stable_gpu_temperature"
echo "resource_preflight one_minute_load=$one_minute_load available_memory_kib=$available_kib gpu_used_mib=$gpu_used_mib gpu_total_mib=$gpu_total_mib gpu_utilization=$gpu_utilization gpu_temperature_c=$gpu_temperature workers_max=$MAX_WORKERS cpu_threads=1"

project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CUDA_VISIBLE_DEVICES=0
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export PYTHONPATH="$project_root/src${PYTHONPATH:+:$PYTHONPATH}"

exec nice -n 10 ionice -c 2 -n 7 "${arguments[@]}"
