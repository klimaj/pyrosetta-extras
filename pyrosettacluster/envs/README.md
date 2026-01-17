# Template Environment Configuration Files

This directory contains template configuration files for creating virtual environments using different environment managers. Note that each environment manager is independent; use only the one that best fits your workflow.

## Supported Environment Managers in PyRosettaCluster

### conda / mamba

- `environment.yml`
  - Create a conda environment:
    - `conda env create -f environment.yml`
  - Create a mamba environment:
    - `mamba env create -f environment.yml`

---

### pixi

- `pixi.toml`
  - Create and install the pixi environment:
    - `pixi install`

---

### uv

- `pyproject.toml`
  - Install and sync dependencies for a uv project:
    - `uv sync`
- `requirements.txt`
  - Alternative dependency specification for a uv project:
    - `uv pip install -r requirements.txt`

---

## Notes

- These files are intended as templates and may need customization for your specific project.
- Do not mix tools in the same environment; choose only one approach per project.
- Python version constraints are defined within each configuration file, but you may wish to pin a specific Python version.
