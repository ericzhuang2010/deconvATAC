# Recreate ignored hidden directories on Ubuntu

This guide recreates the ignored hidden directories at the `deconvATAC` repository root after the Git-tracked files have been synchronized to an Ubuntu machine. It was audited on **2026-08-23**. Run all commands from the repository root unless a section says otherwise.

This is a native rebuild, not a byte-for-byte copy. The source directories contain macOS ARM binaries and absolute macOS paths that cannot run on Ubuntu.

## Scope and migration decisions

The repository currently has four ignored hidden working directories. They are all ignored by the broad `/*/` rule on line 2 of `.gitignore`.

| Directory | Current size | Ubuntu action | Purpose |
|---|---:|---|---|
| `.venv/` | 1.55 GiB | Rebuild | Python environment |
| `.r-lib/` | 700 MiB | Rebuild if using local RCTD | Project-local R package library |
| `.pip-cache/` | 116 MiB | Skip; optionally regenerate | Disposable pip download cache |
| `.pytest_cache/` | 40 KiB | Skip; pytest regenerates it | Disposable test state |

Confirm the ignore rules with:

```bash
git check-ignore -v .venv .r-lib .pip-cache .pytest_cache
```

Other root entries are outside this guide:

- `.git/` is repository metadata created by `git clone`. Do not copy it from the Mac over an existing Ubuntu clone.
- `.DS_Store` is a tracked macOS metadata file, not an environment directory. Git may synchronize it, but there is no Ubuntu recreation step.
- Git-tracked directories arrive through the normal repository sync.
- The ignored, non-hidden `data/` tree has its own [data-recreation guide](<recreate_data_directory (important).md>).

> **Safety:** Do not archive `.venv/` or `.r-lib/` on macOS and extract them on Ubuntu. The current Python interpreter link points into `/Users/.../cpython-3.11-macos-aarch64`, and both environments contain compiled Darwin/Mach-O libraries. Build new directories on the destination machine. If either name already exists there, move it aside and inspect the backup before rebuilding.

## Known source environment

The most important versions on the source machine are:

| Component | Source version |
|---|---|
| Python | CPython 3.11.15 |
| Virtual-environment creator | uv 0.11.7 |
| R | 4.6.1 |
| Bioconductor | 3.23 |
| `SingleCellExperiment` | 1.34.0 |
| `S4Vectors` | 0.50.1 |
| `spacexr` | 2.2.1, commit `9f5dc33c8060f946c6072a138b70e189636e1435` |

Use Python 3.11 on Ubuntu. The package metadata allows older Python versions, but the current environment and ShapeMix workflows were tested with Python 3.11.

For the closest R match, install R 4.6.x and Bioconductor 3.23. A different R release can provide a functional environment with its matching Bioconductor release, but its transitive package versions will differ. The repository does not currently contain an `renv.lock`, so the complete 104-package R library cannot be reconstructed byte-for-byte from a lock file.

## 1. Check the Ubuntu target

Start with a clone or synchronized checkout at the same commit as the source:

```bash
cd /path/to/deconvATAC

git rev-parse HEAD
uname -srm
```

If the machine will use a GPU, also record its driver state before installing JAX or PyTorch:

```bash
nvidia-smi
```

The tracked Cell2location and DestVI configurations currently request a GPU. On a CPU-only machine, use CPU-specific copies of those configurations or change `use_gpu` to `false` in the configuration used for the run.

## 2. Install Ubuntu build prerequisites

R must be installed before the Python project because the base Python dependencies include `rpy2`.

The following packages cover the compilers, R headers, numerical libraries, and common system headers needed by this dependency tree:

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential gfortran pkg-config git curl \
  r-base r-base-dev \
  libblas-dev liblapack-dev \
  libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev \
  libfontconfig1-dev libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev \
  libharfbuzz-dev libfribidi-dev
```

Ubuntu's default `r-base` version depends on the Ubuntu release. If close R/Bioconductor parity matters, follow the [official CRAN Ubuntu instructions](https://stat.ethz.ch/CRAN/bin/linux/ubuntu/) for that Ubuntu release and install R 4.6.x instead of accepting a different distribution version.

If the configured Ubuntu repositories provide Python 3.11, install its virtual-environment and development packages as well:

```bash
sudo apt-get install -y python3.11 python3.11-venv python3.11-dev
```

Package availability varies by Ubuntu release. The `uv` route in the next section can install Python 3.11.15 when those packages are unavailable.

Confirm that R and its development configuration are visible before continuing:

```bash
R --version
R RHOME
```

## 3. Rebuild `.venv/`

Create the environment with Python 3.11.15. If `uv` is not already available, install it using an organization-approved method from the [official uv installation guide](https://docs.astral.sh/uv/getting-started/installation/). The source environment used uv 0.11.7; its versioned standalone installer is:

```bash
curl -LsSf https://astral.sh/uv/0.11.7/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
```

Then create the environment:

```bash
uv python install 3.11.15
uv venv --python 3.11.15 --seed .venv
```

Alternatively, use an Ubuntu-provided or otherwise managed Python 3.11 interpreter:

```bash
python3.11 -m venv .venv
```

Activate the environment:

```bash
source .venv/bin/activate
python -m pip install --upgrade pip
export PYTHONPATH="$PWD/src"
```

To use a project-local pip cache, set the following **before** the package installation commands. Otherwise, omit it and let pip use its normal Ubuntu cache location.

```bash
export PIP_CACHE_DIR="$PWD/.pip-cache"
```

At the audited commit, `pyproject.toml` declares `license = {file = "LICENSE"}`, but the root `LICENSE` file is absent. A fresh editable install can therefore fail while the build metadata is being prepared. Until that repository issue is fixed, install the declared baseline dependencies directly and use `PYTHONPATH=src`:

```bash
python -m pip install \
  numpy==1.25.2 pandas==2.3.3 scipy==1.11.3 \
  anndata==0.9.2 muon==0.1.5 scanpy==1.9.5 \
  h5py==3.16.0 matplotlib==3.10.9 PyYAML==6.0.3 \
  rpy2==3.5.14 anndata2ri==1.3.1 session-info==1.0.1 \
  pytest pytest-cov \
  leidenalg==0.12.0 igraph==1.0.0 pyarrow
```

Once the missing `LICENSE` issue is fixed in the synchronized checkout, the shorter maintained baseline command is:

```bash
python -m pip install -e ".[test,simulation]"
```

The direct fallback makes imports, scripts, and tests work with `PYTHONPATH=src`, but it does not install the `deconvATAC` distribution metadata. Until the editable install is repaired, version/provenance code may report the project as uninstalled or use a fallback version such as `0+unknown`.

Install only the optional method groups needed on the Ubuntu machine. These direct commands also work while the editable project install is blocked:

```bash
# Cell2location
python -m pip install \
  jax==0.4.23 jaxlib==0.4.23 scvi-tools==1.0.3 cell2location==0.1.4

# DestVI
python -m pip install jax==0.4.23 jaxlib==0.4.23 scvi-tools==1.0.3

# Tangram
python -m pip install tangram-sc==1.0.4

# PBMC preparation
python -m pip install celltypist==1.7.1

# ShapeMix
python -m pip install "torch>=2.11,<2.12" "pysam>=0.24,<0.25"

# Notebook and development tools
python -m pip install ipykernel ipython pre-commit "twine>=4.0.2"
```

After the repository metadata issue is fixed, the corresponding extras remain available as `.[cell2location]`, `.[destvi]`, `.[tangram]`, `.[pbmc]`, `.[shapemix]`, `.[doc]`, and `.[dev]`.

For a CUDA machine, select Linux wheels compatible with that machine's NVIDIA driver and CUDA stack. Do not reuse the Mac's PyTorch or JAX installations.

Do **not** use `requirements.txt` verbatim as the migration source. It is a stale local environment snapshot and contains an absolute macOS path such as `deconvATAC @ file:///Users/...`. The dependency declarations in `pyproject.toml`, reflected in the fallback commands above, are the maintained cross-platform source.

## 4. Rebuild `.r-lib/` for local RCTD

This section is needed only when runs use `configs/methods/rctd_local.yaml`. That configuration points to `.r-lib` with a relative path, so launch RCTD from the repository root.

Create the library and install its direct requirements. The following closest-match route requires R 4.6.x and uses Bioconductor 3.23 plus the source machine's pinned `spacexr` commit:

```bash
mkdir -p .r-lib

R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla - <<'RS'
project_lib <- normalizePath(Sys.getenv("R_LIBS_USER"), mustWork = TRUE)
.libPaths(c(project_lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (getRversion() < "4.6" || getRversion() >= "4.7") {
    stop("Closest-version recreation requires R 4.6.x; found ", getRversion())
}

install.packages(c("BiocManager", "remotes"), lib = project_lib)
BiocManager::install(
    version = "3.23",
    lib = project_lib,
    ask = FALSE,
    update = FALSE
)
BiocManager::install(
    c("S4Vectors", "SingleCellExperiment"),
    lib = project_lib,
    ask = FALSE,
    update = FALSE
)
remotes::install_github(
    "dmcable/spacexr@9f5dc33c8060f946c6072a138b70e189636e1435",
    lib = project_lib,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never"
)
RS
```

If Ubuntu must use a different R release, omit both the `getRversion()` guard and the first `BiocManager::install(...)` block that sets `version = "3.23"`; the following package install will use the Bioconductor release compatible with that R version. Keep the `spacexr` Git commit pinned. This produces a functional rebuild but does not preserve the source dependency versions.

The current `.r-lib/` supports RCTD but does **not** contain Giotto. SpatialDWLS uses `configs/methods/spatialdwls.yaml`, whose `r_lib_path` is `null`, and requires a separately managed Giotto plus `data.table` installation. Rebuilding `.r-lib/` as above does not enable SpatialDWLS.

## 5. Handle the cache directories

Neither cache should be copied.

Pip normally uses `~/.cache/pip` on Ubuntu. To retain the source machine's project-local cache layout, set `PIP_CACHE_DIR` as shown before the installation commands in section 3; pip will create `.pip-cache/` as needed.

That setting is optional and affects performance only. The cache is not an environment specification.

Do not create `.pytest_cache/` manually. The first pytest run creates it automatically.

## 6. Validate the rebuilt environments

Verify that Python resolves to a Linux interpreter and that package requirements are consistent:

```bash
readlink -f .venv/bin/python
.venv/bin/python --version
.venv/bin/python -m pip check

PYTHONPATH=src .venv/bin/python - <<'PY'
import platform
import sys

import anndata
import anndata2ri
import deconvatac
import numpy
import rpy2.robjects
import scanpy
import scipy

assert sys.platform.startswith("linux"), sys.platform
print("Python:", sys.version.split()[0])
print("platform:", platform.platform())
print("NumPy:", numpy.__version__)
print("SciPy:", scipy.__version__)
print("AnnData:", anndata.__version__)
print("Scanpy:", scanpy.__version__)
PY
```

The resolved interpreter path must not point under `/Users/` or contain `macos`/`darwin`.

If `.r-lib/` was rebuilt, validate both R and the Python-to-R bridge:

```bash
R_LIBS_USER="$PWD/.r-lib" Rscript --vanilla - <<'RS'
library(S4Vectors)
library(SingleCellExperiment)
library(spacexr)

expected_spacexr_sha <- "9f5dc33c8060f946c6072a138b70e189636e1435"
installed_spacexr_sha <- packageDescription("spacexr")$RemoteSha

stopifnot(
    packageVersion("spacexr") == "2.2.1",
    identical(installed_spacexr_sha, expected_spacexr_sha)
)

cat("R:", R.version.string, "\n")
cat("Bioconductor:", as.character(BiocManager::version()), "\n")
cat("spacexr:", as.character(packageVersion("spacexr")), "\n")
cat("spacexr commit:", installed_spacexr_sha, "\n")
RS

PYTHONPATH=src .venv/bin/python - <<'PY'
from rpy2 import robjects

robjects.r('.libPaths(c(normalizePath(".r-lib"), .libPaths()))')
print(robjects.r('R.version.string')[0])
print("spacexr", robjects.r('as.character(packageVersion("spacexr"))')[0])
PY
```

Finally, run the maintained test suite. The full suite imports PyTorch, so install the ShapeMix optional dependencies from section 3 first even if production runs will not use ShapeMix. The test run also recreates `.pytest_cache/`:

```bash
PYTHONPATH=src .venv/bin/python -m pytest -q
```

## Final checklist

- The Ubuntu checkout is at the intended Git commit.
- `.venv/bin/python` resolves to a Linux Python 3.11 interpreter.
- Python packages were installed from the `pyproject.toml` declarations or the documented direct fallback, not copied from macOS.
- `python -m pip check` reports no broken requirements.
- `.r-lib/` was rebuilt natively if local RCTD is needed.
- `spacexr` reports version 2.2.1 and the install used the pinned Git commit.
- `.pip-cache/` was skipped or regenerated locally.
- `.pytest_cache/` was regenerated by pytest.
- `.git/` was created by the destination clone and was not overwritten.
- The test suite passes on Ubuntu.
