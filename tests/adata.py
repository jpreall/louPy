from anndata import AnnData
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

from loupy.barcodes import load_valid_10x_barcodes


def create_default_anndata(n_bcs=20, n_features=100):
    adata = create_basic_adata(n_bcs, n_features, valid_barcodes=True)
    add_obsm(adata, "proj1")
    add_cluster(adata, "cluster1")
    return adata


def add_obsm(adata, name):
    n_bcs = adata.shape[0]
    proj = np.random.random_sample(size=(n_bcs, 2))
    adata.obsm[name] = proj


def add_cluster(adata, name):
    clusters = list(map(str, range(adata.shape[0])))
    adata.obs[name] = clusters
    adata.obs[name] = adata.obs[name].astype(dtype="category")


def random_barcode(size=16):
    bc = "".join(np.random.choice(list("ACTG"), size=size, replace=True))
    return f"{bc}-1"


def create_basic_adata(n_rows, n_cols, valid_barcodes=True):
    dense = np.random.randint(0, 101, size=(n_rows, n_cols))
    mat = csc_matrix(dense)

    features = [f"feature-{k}" for k in range(n_cols)]
    barcodes = [
        random_barcode() if valid_barcodes else f"bc-{k}" for k in range(n_rows)
    ]

    return AnnData(
        mat, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=features)
    )