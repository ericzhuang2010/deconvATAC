#!/usr/bin/env bash
set -euo pipefail

readonly ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

readonly STATE_ROOT="data/work/preprocessing/shapemix_full_evaluation_v2"
readonly DRIVER_PID_FILE="$STATE_ROOT/driver.pid"
readonly DRIVER_LOG="$STATE_ROOT/driver.log"
readonly PBMC_RUNS="results/sensitivity/shapemix_pbmc_stress_v2/shapemix_pbmc_stress_protocol_v2_cuda/runs.csv"
readonly INTERVAL_SECONDS="${1:-300}"

if [[ ! -s "$DRIVER_PID_FILE" ]]; then
    printf 'Missing driver PID file: %s\n' "$DRIVER_PID_FILE" >&2
    exit 1
fi
readonly DRIVER_PID="$(tr -d '[:space:]' < "$DRIVER_PID_FILE")"
if [[ ! "$DRIVER_PID" =~ ^[0-9]+$ ]]; then
    printf 'Invalid driver PID: %s\n' "$DRIVER_PID" >&2
    exit 1
fi

while :; do
    timestamp="$(date --iso-8601=seconds)"
    alive=0
    if kill -0 "$DRIVER_PID" 2>/dev/null; then
        alive=1
    fi
    rows=0
    if [[ -s "$PBMC_RUNS" ]]; then
        rows="$(( $(wc -l < "$PBMC_RUNS") - 1 ))"
    fi
    load1="$(awk '{print $1}' /proc/loadavg)"
    available_kib="$(awk '$1 == "MemAvailable:" {print $2}' /proc/meminfo)"
    gpu="$(nvidia-smi --query-gpu=memory.used,utilization.gpu,temperature.gpu --format=csv,noheader,nounits | tr -d ' ' | head -n 1)"
    stage="$(tail -n 1 "$DRIVER_LOG" 2>/dev/null || true)"
    process_totals="$(ps -s "$DRIVER_PID" -o %cpu=,rss= 2>/dev/null | awk '{cpu += $1; rss += $2} END {printf "%.1f/%d", cpu, rss}')"
    printf 'monitor_v2 timestamp=%s driver_pid=%s alive=%s cpu_pct_rss_kib=%s pbmc_rows=%s load1=%s mem_available_kib=%s gpu_used_util_temp=%s stage=%q\n' \
        "$timestamp" "$DRIVER_PID" "$alive" "$process_totals" "$rows" "$load1" "$available_kib" "$gpu" "$stage"
    if [[ "$alive" -eq 0 ]]; then
        break
    fi
    sleep "$INTERVAL_SECONDS"
done
