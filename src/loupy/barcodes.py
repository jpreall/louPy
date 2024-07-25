import gzip
from importlib import resources

from pandas import Index


def load_valid_10x_barcodes():
    with resources.path("loupy.barcodes", "3M-3pgex-may-2023.txt.gz") as f:
        with gzip.open(f, "r") as fin:
            bcs = list(fin)
    return Index(bcs)