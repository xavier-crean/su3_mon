*cached_data*# Simplicity of confinement in SU(3) Yang-Mills $-$ Analysis Code

## Contents

- Simplicity of confinement in SU(3) Yang-Mills $-$ Analysis Code
  - [Contents](#contents)
  - [Overview](#overview)
  - [Setup](#setup)
    - [Requirements](#requirements)
    - [Installation](#installation)
    - [Extract data](#extract-data)
    - [Binder](#binder)
  - [Processing the raw lattice configuration data (Step 1)](#processing-the-raw-lattice-configuration-data-step-1)
  - [Performing the analysis (Step 2)](#performing-the-analysis-step-2)
      - [Cached data files](#cached-data-files)

## Overview

This repository accompanies the paper [Simplicity of confinement in SU(3) Yang-Mills][paper]. It contains code for

1. [computing betti numbers of Abelian monopole current graphs][step_1] from SU(3) lattice gauge field configurations,
2. [producing plots and results][step_2] included in the [paper][paper].

Alongside this repository, there exists [an accompanying data release][data] where [Step 1][step_1] in the pipeline has been pre-computed from configurations generated using [accompanying Monte-Carlo code][monte_carlo]. Therefore, to reproduce analysis from the paper using the data release, the user is referred to [Step 2][step_2] ([setup instructions](#setup) are included below).

Using Binder, a Docker container has been constructed for reproduction of results from the paper (without the need to install Python or dependencies) and can be accessed online ([see below](#binder)).

## Setup

### Requirements

The code has been tested with Python 3.10. 

Dependencies are documented in `environment.yml` and are most easily managed via a Conda environment. 

The commands below refer to a Linux or compatible environment and need to be run from the root directory of the repository.

### Installation

* Download the repository
* Create a new Conda environment with the necessary dependencies using

      conda env create -f environment.yml
    (Alternatively, create a Python environment with the listed packages.)

### Extract data

* Download the `data.tar.gz` from [the accompanying data release][data] and extract into the root directory of the repository using

        tar -xf data.tar.gz
* To use **pre-computed** [cached data files](#cached-data-files), for fast reproduction of the results in the paper, download `cached_data.tar.gz` from the [data release][data] and extract into the root directory of the repository using

        tar -xf cached_data.tar.gz

---

### Binder

Using [mybinder.org][binder], a Docker container with [necessary dependencies](#requirements) has been constructed for reproduction of the results in the paper. Note that the [data release][data] must be uploaded and extracted into the root directory (see [data release setup](#extract-release)).

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/xavier-crean/su3_mon/HEAD)

## Processing the raw lattice configuration data (Step 1)

The directory `src/` contains a python script and a python module:

* `src/betti_compute.py` for computing the zeroth and first order homology of monopole current graphs from Abelian projected configuration data in effect two U(1) lattice gauge field configurations. It takes six command line arguments: lattice size, $\beta$ (to 4d.p.), number of samples, number of sweeps between measurement, number of sweeps to discard (for thermalisation) and number of parallel computations. E.g.,

      python3 src/betti_compute.py 8.64.64.64 5.8900 400 2000 10000 5
    The output is saved into directory `data/observables/betti/Nt.Nx.Ny.Nz/` in `.csv` format.
* `src/configuration_u1.py` a module used to parse lattice configurations from raw configuration data files.

In the accompanying [data release][data], in `data/cnfgs/8.32.32.32/6.0000/`, a single Abelian projected configuration with $N_t = 8$ and $N_s = 32$ at $\beta = 6.0$, in the format of two U(1) lattice gauge field configurations, outputted by the [Monte Carlo and gauge fixing code][monte_carlo], is included as an example for the user. This data is to be used as the input for `src/betti_compute.py`.

## Performing the analysis (Step 2)

The analysis and plots from the paper are reproduced using jupyter notebooks. The notebooks are of two types:
1. `analysis/Nt=X.ipynb` reproduces the plots from the paper of the observables $\rho_0$, $\rho_1$ and $\lambda$ (and their susceptibilities) from raw per configuration Betti number data. Figures are saved in the directory `reports/Nt=X/figures/`.
2. `analysis/Nt=X_rw.ipynb` performs a multiple histogram reweighting of the observables and a finite-size scaling analysis. It uses parameters set by our analysis that reproduce the results from the paper. Figures are saved in the directory `reports/Nt=X/figures/`.

An accompanying module `reweight.cpp`, written in C++, performs the multiple histogram reweighting part of the analysis (with dependency `ranlxs.cpp`). This has been written with a python wrapper via [pybind11][pybind] as an API for use with the jupyter notebooks. The C++ must be compiled before execution and so requires the installation of a compiler. In a Unix-like environment, the standard [g++ compiler][g++] is appropriate. `setup.py` implements the compilation command via

      python3 setup.py build_ext --inplace

For ease, this is copied as a cell in the relevant jupyter notebooks. However, please note that other non-Unix-like distributions, e.g., MS Windows may require alternative compilation procedure.

### Cached data files

Intermediary files, for caching the results of the histogram reweighting code so that the code may be run more quickly second time round, are stored in the directory `cached_data/`. To use the files in `cached_data.tar.gz` (from the [data release][data]), it is important to ensure the **original** random seeds are used (as found in `analysis/Nt=X_rw.py`). With the cached data files, the code takes a few minutes to run; without them, it takes on the order of hours.

**N.B.** if `cached_data.tar.gz` is extracted into the root, then `analysis/Nt=X_rw.py` will read from this directory. If you are trying to reproduce the analysis using `data.tar.gz`, you will need to delete or rename the directory `cached_data/`.

---

[paper]: https://arxiv.org/abs/2602.10088
[data]: https://doi.org/10.5281/zenodo.18395360
[monte_carlo]: https://doi.org/10.5281/zenodo.18171195
[mc_README]: src/data/README
[mc_AUTHORS]: src/data/AUTHORS
[mc_install]: src/data/INSTALL
[binder]: https://mybinder.org/
[hdf5]: https://www.hdfgroup.org/solutions/hdf5
[step_1]: #code-to-generate-lattice-configuration-data-step-1
[step_2]: #processing-the-raw-lattice-configuration-data-step-2
[pybind]: https://pybind11.readthedocs.io/en/stable/index.html
[g++]: https://gcc.gnu.org/