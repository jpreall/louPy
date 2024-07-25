import re
import gzip
from importlib import resources

import numpy as np
import pandas as pd
from anndata import AnnData


def load_valid_10x_barcodes():
    with resources.path("loupy.barcodes", "3M-3pgex-may-2023.txt.gz") as f:
        with gzip.open(f, "r") as fin:
            bcs = list(fin)
    return pd.Index(bcs)


def swap_barcodes(obs_names: pd.Index | list | np.ndarray) -> np.ndarray:
    """
    In the case that one only has access to a raw .mtx file or other strange circumstance,
    we need a mechanism to correct for missing, valid 10x barcodes in adata.obs_names.

    This takes a sequence of observation names (e.g. 0...n, barcode-1...barcode-n, etc.)
    and returns a sequence of valid 10x 3' v4 GEX barcodes.
    """
    n_obs = len(obs_names)
    valid_barcodes = load_valid_10x_barcodes()
    n_valid_barcodes = len(valid_barcodes)
    copy_number = max(n_obs // n_valid_barcodes, 1)

    valid_prefixes = np.repeat(valid_barcodes, copy_number)
    valid_suffixes = np.tile(np.arange(1, copy_number+1, dtype=int), n_valid_barcodes)
    new_obs =  pd.Index(valid_prefixes) + "-" + pd.Index(valid_suffixes).astype(str)

    return new_obs.to_numpy()[:n_obs]


def validate_barcodes(barcodes: np.ndarray | pd.Index | list[str]) -> bool:
    bc_pattern = re.compile("^([ACTG]{6,})(-[0-9]+?)?$")
    return all([bc_pattern.match(bc) for bc in barcodes])


def format_barcodes(adata: AnnData) -> np.ndarray:
    barcodes = adata.obs_names
    if validate_barcodes(barcodes):
        return barcodes.to_numpy()

    bc_pattern = re.compile("^(.*?)(_|-|:)?([ACTG]{6,})(-[0-9]+)?(_|-|:)?(.*?)$")
    expanded = barcodes.str.extractall(bc_pattern).fillna("")
    expanded.columns = ["prefix", "sep1", "barcode", "dashnum", "sep2", "suffix"]
    if expanded.prefix.any():
        expanded.prefix = "-" + expanded.prefix
    if expanded.suffix.any():
        expanded.suffix = "-" + expanded.suffix

    formated = expanded.barcode + expanded.dashnum + expanded.prefix + expanded.suffix
    if len(formated) != len(barcodes):
        return barcodes.to_numpy()

    return formated.to_numpy()