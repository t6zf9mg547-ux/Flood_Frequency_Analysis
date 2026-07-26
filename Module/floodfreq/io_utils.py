"""
Input/output helpers: read an annual-maximum series from CSV or Excel, and
write a FloodFrequencyAnalysis's results back out to a formatted .xlsx
report (so results stay easy to share with people who work in Excel).

Also provides project-layout helpers for the standard folder structure:

    <project_root>/
        Data/<CaseName>.csv     - input data, one file per case
        Module/                 - this package + scripts (Module/floodfreq/...)
        Output/<CaseName>/      - CSV outputs for that case
        Plot/<CaseName>/        - PNG plots for that case

`resolve_case()` derives all of these paths from a case name, given the path
of the calling script (pass `__file__` from a script living in Module/).
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter


@dataclass
class CasePaths:
    case_name: str
    project_root: Path
    data_csv: Path
    output_dir: Path
    plot_dir: Path


def project_root_from(module_file) -> Path:
    """
    Given __file__ of a script living directly under <project_root>/Module/,
    return <project_root>. Walks up until a directory containing both a
    'Data' and a 'Module' subdirectory is found, so it also works for
    scripts nested one level deeper inside Module/.
    """
    p = Path(module_file).resolve().parent
    for candidate in [p] + list(p.parents):
        if (candidate / "Data").is_dir() and (candidate / "Module").is_dir():
            return candidate
    # Fallback: assume the script sits directly in Module/
    return p.parent


def resolve_case(case_name: str, module_file) -> CasePaths:
    """
    Resolve the standard Data/Module/Output/Plot paths for one case, and
    create Output/<CaseName>/ and Plot/<CaseName>/ if they don't exist yet.
    """
    root = project_root_from(module_file)
    data_csv = root / "Data" / f"{case_name}.csv"
    output_dir = root / "Output" / case_name
    plot_dir = root / "Plot" / case_name
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    return CasePaths(case_name, root, data_csv, output_dir, plot_dir)


def read_series(path, value_col=None, year_col=None, sheet_name=0):
    """
    Read an annual-maximum series from a CSV or Excel file.

    If value_col/year_col are not given, the function guesses: the first
    all-numeric column with more than one distinct value is treated as the
    flood series, and a column that looks like a 4-digit year is used as
    the year index if present.
    """
    path = str(path)
    if path.lower().endswith((".csv",)):
        df = pd.read_csv(path)
    else:
        df = pd.read_excel(path, sheet_name=sheet_name)

    if value_col is None:
        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        if not numeric_cols:
            raise ValueError("No numeric column found; pass value_col explicitly.")
        # prefer a column that is NOT a plausible year column
        candidates = [c for c in numeric_cols
                      if not (df[c].between(1800, 2100).mean() > 0.9)]
        value_col = candidates[0] if candidates else numeric_cols[-1]

    if year_col is None:
        for c in df.columns:
            if pd.api.types.is_numeric_dtype(df[c]) and df[c].between(1800, 2100).mean() > 0.9:
                year_col = c
                break

    sub = df[[year_col, value_col]] if year_col else df[[value_col]]
    sub = sub.dropna()
    values = sub[value_col].to_numpy(dtype=float)
    years = sub[year_col].to_numpy() if year_col else None
    return values, years


def write_report(path, ffa, return_periods=(2, 5, 10, 20, 25, 50, 100, 200, 500, 1000)):
    """Write descriptive stats, goodness-of-fit ranking, and a quantile
    table to a formatted .xlsx workbook."""
    wb = Workbook()

    header_font = Font(name="Arial", bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="4472C4")
    body_font = Font(name="Arial")
    title_font = Font(name="Arial", bold=True, size=14)

    def _write_table(ws, df, start_row=1, start_col=1):
        for j, col in enumerate(df.columns):
            cell = ws.cell(row=start_row, column=start_col + j, value=str(col))
            cell.font = header_font
            cell.fill = header_fill
            cell.alignment = Alignment(horizontal="center")
        for i, (_, row) in enumerate(df.iterrows()):
            for j, col in enumerate(df.columns):
                val = row[col]
                if isinstance(val, (np.floating, float)):
                    val = round(float(val), 3)
                elif isinstance(val, (np.integer,)):
                    val = int(val)
                cell = ws.cell(row=start_row + 1 + i, column=start_col + j, value=val)
                cell.font = body_font
        for j, col in enumerate(df.columns):
            ws.column_dimensions[get_column_letter(start_col + j)].width = max(14, len(str(col)) + 2)

    # -- Summary sheet --
    ws = wb.active
    ws.title = "Summary"
    ws["A1"] = f"Flood Frequency Analysis — {ffa.station_id}"
    ws["A1"].font = title_font
    stats = ffa.descriptive_stats()
    rows = [(k, v) for k, v in stats.items()]
    stats_df = pd.DataFrame(rows, columns=["Statistic", "Value"])
    _write_table(ws, stats_df, start_row=3)

    # -- Goodness of fit --
    ws2 = wb.create_sheet("Goodness of fit")
    gof = ffa.goodness_of_fit_table()
    _write_table(ws2, gof, start_row=1)

    # -- Quantile table --
    ws3 = wb.create_sheet("Quantiles")
    q = ffa.quantile_table(return_periods=return_periods)
    _write_table(ws3, q, start_row=1)

    # -- Raw data --
    ws4 = wb.create_sheet("Data")
    data_df = pd.DataFrame({"year": ffa.years, "Q": ffa.data}) if ffa.years is not None \
        else pd.DataFrame({"Q": ffa.data})
    _write_table(ws4, data_df, start_row=1)

    wb.save(path)
    return path
