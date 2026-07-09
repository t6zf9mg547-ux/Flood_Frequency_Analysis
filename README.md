# Python Project Template

A reusable starting point for Python projects using [uv](https://docs.astral.sh/uv/) for dependency management, geared toward geospatial/data-analysis scripts (geopandas, pandas, shapely, etc.).

## What's included
Template/
├── Data/            # place input data here (not tracked in git)
├── Module/          # your Python scripts go here
├── Output/          # generated results (not tracked in git)
├── Plot/            # generated figures/plots (not tracked in git)
├── pyproject.toml   # project metadata + dependencies (uv-managed)
└── .gitignore       # excludes venv, cache files, Data/Output/Plot, OS junk, etc.

## How to start a new project from this template

```bash
cp -r Template ~/path/to/NewProjectName
cd ~/path/to/NewProjectName
```

Then edit `pyproject.toml`:
- Update `name` (e.g. `"extract-drainage-system"`)
- Update `description` to describe the new project

Install dependencies and create the virtual environment:
```bash
uv sync
```

Run any script using the project's environment:
```bash
uv run python Module/your_script.py
```

## Default dependencies

- geopandas
- pandas
- fiona
- shapely
- numpy
- pyproj
- pyarrow (needed if reading/writing GeoParquet files)

Add or remove packages as the new project needs:
```bash
uv add <package>
uv remove <package>
```

## Notes

- `tkinter` (used for file/folder selection dialogs in some scripts) is part of the Python standard library and does not need to be added as a dependency — just make sure your Python installation includes it (standard on macOS python.org / Homebrew builds).
- GDAL must be installed at the system level for `fiona`/`geopandas` to work: `brew install gdal`.
- `[tool.uv] package = false` in `pyproject.toml` marks this as a scripts project rather than an installable package — required so `uv sync`/`uv add` don't try (and fail) to build a wheel.
- `Data/`, `Output/`, and `Plot/` are excluded from git by default, since they typically hold large or regenerated files. Adjust `.gitignore` per-project if you want any of them tracked.
