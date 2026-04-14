# APWP Strat Models

Site-based apparent polar wander path (APWP) workflows for the Keweenawan track, combining paleomagnetic virtual geomagnetic pole (VGP) compilation, stratigraphic height constraints, Bayesian age-depth modeling, and uncertainty propagation through bootstrap APWP construction.

This repository brings together:

- site- and pole-level paleomagnetic compilations
- stratigraphic height calculations for key successions
- Bayesian age-model notebooks built on Chron.jl
- Python utilities for APWP calculation, visualization, and uncertainty propagation

## Overview

The core goal of this project is to build APWPs that are informed not only by paleomagnetic directions and pole statistics, but also by stratigraphic ordering and age-model uncertainty. In practice, the repository supports a workflow that:

1. compiles VGP- and site-level paleomagnetic data for the Keweenawan track
2. assigns or reconstructs stratigraphic heights for samples and formations where possible
3. generates section-specific Bayesian age-depth models in Julia using Chron.jl
4. propagates age and directional uncertainty through bootstrap or running-mean APWP calculations
5. exports tabular APWP products and figures for interpretation and paleogeographic analysis


## Repository Layout

| Path | Purpose |
| --- | --- |
| `code/` | Main notebooks for compilation, APWP construction, and paleogeographic analysis |
| `code/age_models/` | Julia-based age-model notebooks for individual stratigraphic successions |
| `code/strat_height/` | Notebooks and helper code for stratigraphic height calculations |
| `code/vgptools/` | Reusable Python utilities for APWP statistics, compilation, and visualization |
| `data/pmag_compiled/` | Compiled paleomagnetic datasets and study-level source folders |
| `data/age_models_output/` | Posterior age-model outputs exported from Chron.jl workflows |
| `data/models_shapes/` | Rotation files and paleogeographic shape files used for reconstructions |
| `code/apwp_output/` | Exported APWP tables |
| `figure/` | Generated figures |
| `environment_APWP_age_model.yml` | Recommended conda environment specification |
| `environment_APWP_age_model_full.yml` | More tightly pinned environment for reproducibility |

## Notebook Guide

The notebooks are organized roughly in workflow order.

| Notebook | Role |
| --- | --- |
| `code/00_MagIC_files.ipynb` | Prepare or inspect MagIC-format inputs |
| `code/01_VGP_compilation.ipynb` | Compile Keweenawan-track VGP and site-level paleomagnetic data |
| `code/02_APWP_age_model.ipynb` | Construct APWPs while propagating age-model uncertainty |
| `code/03_Keweenawan_VGP_paleogeography.ipynb` | Paleogeographic interpretation and visualization of the resulting track |

Age-model notebooks in `code/age_models/` are organized by stratigraphic succession, including Mamainse Point, Michipicoten, NSVG, Oronto Group, Osler Volcanics, and PLV.

## Data and Methods Notes

Several workflows in this repository extend earlier Keweenawan-track compilation efforts by incorporating site-based information and stratigraphic context. The VGP compilation notebooks explicitly track sample or site heights where available so that posterior age-depth models can be sampled during APWP bootstrap resampling.

The Python APWP utilities in `code/vgptools/` calculate, among other quantities:

- Fisher mean poles and associated statistics
- moving-window or weighted APWP estimates with a triangular kernel
- apparent polar wander rates between successive mean poles
- shape descriptors for pole distributions
- pseudo-VGP or directional resampling products used in uncertainty propagation

The Julia helpers in `code/age_models/MCR_age_models.jl` wrap Chron.jl workflows for building stratigraphic age models, exporting posterior age distributions, and plotting age-depth and accumulation-rate results.

## Environment Setup

This repository includes two conda environment files:

- `environment_APWP_age_model.yml`: recommended environment for most users
- `environment_APWP_age_model_full.yml`: more tightly pinned environment for higher reproducibility on similar platforms

The recommended environment currently specifies Python 3.14 and includes the main scientific stack plus `pmagpy==4.3.11`.

### Create the Python environment

```bash
conda env create -f environment_APWP_age_model.yml
conda activate APWP_age_model
```

### Create from the fully pinned spec

```bash
conda env create -f environment_APWP_age_model_full.yml
conda activate APWP_age_model
```

### Update an existing environment

```bash
conda env update -n APWP_age_model -f environment_APWP_age_model.yml --prune
```

### Launch notebooks

```bash
jupyter lab
```

After launching Jupyter, open notebooks in `code/` and select the `APWP_age_model` kernel.

## Julia Age-Model Environment

The age-model construction notebooks in `code/age_models/` use Julia rather than Python. The Julia project files live in `code/Project.toml` and `code/Manifest.toml` and include Chron, SubsidenceChron, CSV, DataFrames, and related dependencies.

If you want to run the Julia notebooks locally, instantiate the Julia environment from the `code/` directory:

```bash
julia --project=code -e 'using Pkg; Pkg.instantiate()'
```

The checked-in manifest was generated with Julia 1.11.6.

You will also need a Julia notebook kernel available in your Jupyter environment if you want to execute the notebooks interactively.

## Typical Workflow

For a fresh reproduction of the analysis, the practical order is:

1. create the Python environment and verify the Jupyter kernel
2. compile or inspect the paleomagnetic data in `code/01_VGP_compilation.ipynb`
3. build or review section-specific age models in `code/age_models/`
4. confirm posterior age tables in `data/age_models_output/`
5. run the `code/02_APWP_age_model` notebook to generate APWPs
6. inspect exported path tables in `code/apwp_output/` and figures in `figure/`
7. use `code/03_Keweenawan_VGP_paleogeography.ipynb` for paleogeographic visualization

## Outputs

The repository already includes representative derived products, including:

- posterior age tables in `data/age_models_output/`
- compiled paleomagnetic datasets in `data/pmag_compiled/`
- APWP export tables in `code/apwp_output/`

One example output table currently present is `code/apwp_output/Keweenawan_track_VGP_APWP.csv`.

## Reproducibility Notes

- Python and Julia workflows are intentionally separated; both environments are needed for full end-to-end reproduction.
- The repository contains both code and intermediate data products, so some notebooks can be explored without rerunning every upstream step.
- Relative paths in the notebooks assume they are run from within the checked-out repository structure.

## Authors
- [Yiming Zhang](https://github.com/duserzym)
- [Nicholas Swanson-Hysell](https://github.com/Swanson-Hysell)
- [Leandro Gallo](https://github.com/LenGallo)
- [Mathew Domeier](https://github.com/matdomeier)
- [Facundo Sapienza](https://github.com/facusapienza21)
- [Diego Alberto Osorio Afanador](https://github.com/Osorio-AfanadorD)

