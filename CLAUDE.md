# Nextflow Conventions

- All Python scripts called from Nextflow processes must live in `bin/` — Nextflow automatically adds `bin/` to `$PATH` in every process, so scripts there can be invoked by name without a full path.
- Do NOT place Nextflow-called scripts in `scripts/` or elsewhere; they will not be found at runtime.
- Always run `chmod +x bin/<script>.py` after creating a new script in `bin/`.

# Python Environment

- Use the virtual environment at `.venv/` (activate with `source .venv/bin/activate`)
- Python version: 3.11
- Load modules rust, arrow and bedtools before running python
- Package manager: uv (use `uv add <pkg>` not pip install)
- Run scripts with: `uv run script.py`
- Run tests with: `uv run pytest`