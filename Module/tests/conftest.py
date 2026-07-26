"""
Shared pytest fixtures. Ensures floodfreq is importable regardless of the
directory pytest is invoked from, and provides the reference validation
dataset embedded directly (NOT read from Data/, which is gitignored and
may not exist on a fresh clone or CI checkout).
"""
import sys
from pathlib import Path
import numpy as np
import pytest

MODULE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(MODULE_DIR))

# The exact station-P1009 series used throughout development to validate
# fitted parameters against the original Excel/VBA workbook (Analyse_
# Frequentielle_V03.xls). Do not change these values -- the regression
# tests in test_distributions.py assert against numbers cross-checked
# against that workbook.
REFERENCE_YEARS = np.array([
    1915, 1916, 1917, 1918, 1919, 1920, 1921, 1922, 1923, 1924, 1925, 1926,
    1927, 1928, 1929, 1930, 1931, 1932, 1933, 1934, 1935, 1936, 1937, 1938,
    1939, 1940, 1941, 1942, 1943, 1944, 1945, 1946, 1947, 1948, 1950, 1952,
    1953, 1954, 1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964,
    1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976,
    1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988,
    1989, 1990, 1991, 1992, 1993, 1994,
])

REFERENCE_DATA = np.array([
    122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205,
    135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240,
    230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 204, 164,
    183, 161, 167, 179, 185, 117, 192, 337, 125, 166, 99.1, 202,
    230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205,
    237, 177, 239, 187, 180, 173, 174, 167, 185, 232, 100, 163,
    203, 219, 182, 184, 118, 155,
], dtype=float)


@pytest.fixture
def reference_data():
    """(values, years) for the validated 78-year reference station."""
    return REFERENCE_DATA.copy(), REFERENCE_YEARS.copy()
