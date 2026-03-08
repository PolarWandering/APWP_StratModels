# Site-based APWPs with stratigraphic age models

## Python environment setup

This repository includes two conda environment files:

- `environment_APWP_age_model.yml`: minimal, cross-platform spec (recommended for most users)
- `environment_APWP_age_model_full.yml`: fully pinned spec for high reproducibility on similar platforms

### Create the environment (recommended)

```bash
conda env create -f environment_APWP_age_model.yml
conda activate APWP_age_model
```

### Create from fully pinned spec

Use this if you want package versions as close as possible to the original setup.

```bash
conda env create -f environment_APWP_age_model_full.yml
conda activate APWP_age_model
```

### Update an existing environment

If you already created `APWP_age_model` and want to sync it to a YAML file:

```bash
conda env update -n APWP_age_model -f environment_APWP_age_model.yml --prune
```

### Notebook kernel

After activating the environment, open notebooks in `code/` and select the `APWP_age_model` kernel.

Note: age-model construction notebooks in `code/age_models/` use Julia/Chron.jl workflows and are separate from this Python environment.

## Code

### age_models

The notebooks in this folder utilize a julia kernel to develop Bayesian age models for the stratigraphic successions using the [chron.jl](https://github.com/brenhinkeller/Chron.jl) package. Each succession has an individual .ipynb notebook that brings in functions from the age_models.jl script.


