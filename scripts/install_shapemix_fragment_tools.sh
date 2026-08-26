#!/usr/bin/env bash
set -euo pipefail

readonly SNAPATAC2_VERSION="2.4.0"
readonly SNAPATAC2_CP311_WHEEL_SHA256="06ced8baef6f6eaef5685bb5b538caa09f040de8fea5e4a2e2100c3c594932a4"
readonly BWA_TAG="v0.7.17"
readonly PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
readonly PYTHON_BASE="${PROJECT_ROOT}/.venv/bin/python"
readonly ENVIRONMENT="${PROJECT_ROOT}/.venv-shapemix-fragments"
readonly TOOL_ROOT="${PROJECT_ROOT}/data/work/preprocessing/gse246791_mouse_brain_reference/tool_source"
readonly WHEEL_ROOT="${PROJECT_ROOT}/data/work/downloads/shapemix_fragment_tools"
readonly WHEEL="${WHEEL_ROOT}/snapatac2-2.4.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
readonly BWA_SOURCE="${TOOL_ROOT}/bwa"

if [[ "${DECONVATAC_RESOURCE_GUARD:-}" != "1" ]]; then
    echo "Run this installer through scripts/run_shapemix_low_impact.sh." >&2
    exit 64
fi

if [[ ! -x "${PYTHON_BASE}" ]]; then
    echo "Tested project Python is missing: ${PYTHON_BASE}" >&2
    exit 1
fi
python_base_version="$("${PYTHON_BASE}" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"
if [[ "${python_base_version}" != "3.11" ]]; then
    echo "SnapATAC2 2.4.0 requires the pinned CPython 3.11 ABI; found ${python_base_version}." >&2
    exit 1
fi

if [[ ! -x "${ENVIRONMENT}/bin/python" ]]; then
    "${PYTHON_BASE}" -m venv "${ENVIRONMENT}"
fi

mkdir -p "${WHEEL_ROOT}" "${TOOL_ROOT}"
if [[ ! -f "${WHEEL}" ]]; then
    "${ENVIRONMENT}/bin/python" -m pip download --no-deps --only-binary=:all: --dest "${WHEEL_ROOT}" "snapatac2==${SNAPATAC2_VERSION}"
fi
observed_wheel_sha256="$(sha256sum "${WHEEL}" | awk '{print $1}')"
if [[ "${observed_wheel_sha256}" != "${SNAPATAC2_CP311_WHEEL_SHA256}" ]]; then
    echo "SnapATAC2 wheel SHA-256 mismatch: ${observed_wheel_sha256}" >&2
    exit 1
fi

"${ENVIRONMENT}/bin/python" -m pip install \
    "${WHEEL}" \
    "pysam==0.24.0" \
    "PyYAML==6.0.3"

if [[ ! -d "${BWA_SOURCE}/.git" ]]; then
    if [[ -e "${BWA_SOURCE}" ]]; then
        echo "Non-Git BWA source path already exists: ${BWA_SOURCE}" >&2
        exit 1
    fi
    git clone --branch "${BWA_TAG}" --depth 1 https://github.com/lh3/bwa.git "${BWA_SOURCE}"
fi
if [[ "$(git -C "${BWA_SOURCE}" describe --tags --exact-match)" != "${BWA_TAG}" ]]; then
    echo "BWA source is not pinned to ${BWA_TAG}." >&2
    exit 1
fi

make -C "${BWA_SOURCE}" -j1
install -m 0755 "${BWA_SOURCE}/bwa" "${ENVIRONMENT}/bin/bwa"

snapatac2_observed="$("${ENVIRONMENT}/bin/python" -c 'import snapatac2; print(snapatac2.__version__)')"
if [[ "${snapatac2_observed}" != "${SNAPATAC2_VERSION}" ]]; then
    echo "Unexpected SnapATAC2 version: ${snapatac2_observed}" >&2
    exit 1
fi
"${ENVIRONMENT}/bin/python" -c 'import h5py, pysam, scipy, yaml'
bwa_help="$("${ENVIRONMENT}/bin/bwa" 2>&1 || true)"
bwa_observed="$(awk '/^Version:/ {print $2}' <<< "${bwa_help}")"
if [[ "${bwa_observed}" != "0.7.17-r1188" ]]; then
    echo "Unexpected BWA version: ${bwa_observed}" >&2
    exit 1
fi

echo "fragment_tools_ready snapatac2=${snapatac2_observed} bwa=${bwa_observed}"
echo "bwa_commit=$(git -C "${BWA_SOURCE}" rev-parse HEAD)"
