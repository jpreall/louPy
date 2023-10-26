import h5py
import numpy as np
import platform
import pkg_resources
import h5py
import os
import re
import subprocess

"""
Useful hack to write a Loupe file from a Scanpy object using LoupeR (10X)
This module writes a LoupeR-compatible hdf5 file, then passes it to the LoupeR executable

Step 1: Download and install LoupeR from 10X Genomics
https://github.com/10XGenomics/loupeR

Step 2: Run "setup" to download LoupeR executable and agree to EULA

Step 3: Copy LoupeR executable file to the same directory as this .py file

"""

script_dir = os.path.dirname(os.path.abspath(__file__))
louper_path = os.path.join(script_dir, "louper")


def create_hdf5(adata, h5path, scanpy_obj_version, LouPy_style=False):
    if os.path.exists(h5path):
        raise ValueError(f"Cannot create h5 file {h5path}, file already exists.")

    with h5py.File(h5path, 'w') as f:
        write_mat(f, adata)
        write_clusters(f, adata)
        write_projections(f, adata)
        metadata = create_metadata(scanpy_obj_version, LouPy_style=LouPy_style)
        # Assuming you have a function to write metadata to the HDF5 file
        write_metadata(f, metadata)

    return "SUCCESS"

def write_mat(f, adata, layer='counts'):
    """
    Writes data matrix from adata object from the specified layer.

    Currently, 'layer' must be a layer in adata.layers
    and must contain raw integer counts. 

    Currenty, barcodes must be in '10X format': eg. ATAGTGCCTGGTACGA-1

    TODO: write a method to sanitize barcodes to allow weirdly formatted barcodes to be Loupified
    """

    def sanitize_barcodes(barcodes_unmodified):
        """
        TODO: actually write a full-features barcode sanitization function
        """
        barcode_sequence = barcodes_unmodified.str.split('-').str[0]
        barcode_suffix = barcodes_unmodified.str.split('-').str[1].str.zfill(4)
        sanitized = barcode_sequence + '-' + barcode_suffix
        sanitized = np.array(sanitized, dtype=str).astype('S')
        return sanitized

    # Check if provided counts layer is valid
    if layer not in adata.layers:
        raise ValueError(f"{layer} layer not found in adata.layers")
    count_mat = adata.layers[layer]

    if not isIntegers(count_mat):
        raise ValueError(f"The provided layer ({layer}) must contain only integer values")
    
    barcodes_unmodified = adata.obs_names.tolist()
    barcodes_formatted = barcodes_unmodified ## REPLACE THIS WITH SANITIZED BARCODE OUTPUT 
    barcode_count = len(barcodes_unmodified)

    ### Features group
    matrix_group = f.create_group("matrix")

    features = adata.var_names.to_numpy(dtype=str).astype('S')
    feature_count = len(features)
    feature_ids = [f"feature_{i}" for i in range(1, len(features) + 1)]
    # Write Features group
    features_group = matrix_group.create_group("features")
    create_str_dataset(features_group, "name", strs=features)
    create_str_dataset(features_group, "id", strs=feature_ids)
    create_str_dataset(features_group, "feature_type", strs=["Gene Expression"] * len(features))
    create_str_dataset(features_group, "_all_tag_keys", strs=np.array([], dtype='S'))
    
    # Write Matrix Group
    create_str_dataset(matrix_group, "barcodes", strs=barcodes_formatted)
    create_str_dataset(matrix_group, "barcodes_unmodified", strs=barcodes_unmodified)
    matrix_group.create_dataset("data", data=count_mat.data.astype('int32'), compression='gzip')
    matrix_group.create_dataset("indices", data=count_mat.indices.astype('int32'), compression='gzip')
    matrix_group.create_dataset("indptr", data=count_mat.indptr.astype('int32'), compression='gzip')
    matrix_group.create_dataset("shape", data=np.array([feature_count, barcode_count]).astype('int32'), compression='gzip')

def write_clusters(f, adata):
    """
    Searches for categorical observations in adata.obs and writes them all to h5 file.
    """
    clusters_group = f.create_group("clusters")
    categorical_columns = adata.obs.select_dtypes(include=['category'])
    if len(categorical_columns.columns) > 16:
        raise ValueError("Too many categorical columns in adata.obs. Limit is 16.")

    for name, cluster in categorical_columns.items():
        
        group = clusters_group.create_group(name)
    
        group.create_dataset("name", data=name, shape=(1,), dtype=f'S{len(name)}')
        
        # Force data into strings
        group_names = cluster.cat.categories.astype('str').tolist()
        create_str_dataset(group, "group_names", strs=group_names)

        assignments = np.array(cluster.cat.codes, dtype='int32')
        group.create_dataset("assignments", data=assignments, compression="gzip")

        group.create_dataset("score", data=np.array([0.0]), compression="gzip")
        create_str_dataset(group, "clustering_type", strs="unknown")

def write_projections(f, adata):
    projections_group = f.create_group("projections")

    for name, projection in adata.obsm.items():
        name = name.replace('X_','')
        n_dim = projection.shape[1]
        is_umap = bool(re.search("umap", name, re.IGNORECASE))
        is_tsne = bool(re.search("tsne", name, re.IGNORECASE))
        is_tsne_dash = bool(re.search("t-sne", name, re.IGNORECASE))

        if is_umap:
            method = "UMAP"
        elif is_tsne or is_tsne_dash:
            method = "t-SNE"
        else:
            method = name
        
        if n_dim == 2:
            group = projections_group.create_group(name)
            create_str_dataset(group, "name", strs=name)
            create_str_dataset(group, "method", strs=method)
            group.create_dataset("data", data=projection.T, compression='gzip')
            
def create_metadata(scanpy_obj_version=None, LouPy_style=False):
    
    if LouPy_style:
        # Get Python version
        python_version = platform.python_version()
        
        # Create metadata dictionary
        meta = {}
        meta["tool"] = "loupePy"
        #meta["tool_version"] = pkg_resources.get_distribution("loupePy").version if pkg_resources.get_distribution("loupePy") else "n/a"
        meta["tool_version"] = "pre-alpha"
        meta["os"] = platform.system()
        meta["system"] = platform.platform()
        meta["language"] = "Python"
        meta["language_version"] = python_version
        
        # Create extra dictionary
        extra = {}
        extra["loupePy_scanpy_version"] = pkg_resources.get_distribution("scanpy").version if pkg_resources.get_distribution("scanpy") else "n/a"
        extra["loupePy_scanpy_object_version"] = scanpy_obj_version if scanpy_obj_version else "n/a"
        extra["loupePy_hdf5_version"] = h5py.version.hdf5_version
        
        # Add extra to meta
        meta["extra"] = extra
    
    else:
        # Hardcoded metadata from R version
        meta = {}
        meta["tool"] = "loupeR"
        meta["tool_version"] = "1.0.0"
        meta["os"] = "macOS Catalina 10.15.7"
        meta["system"] = 'x86_64-apple-darwin17.0 (64-bit)'
        meta["language"] = "R"
        meta["language_version"] = "4.2.3"

        extra = {}
        extra["loupeR_hdf5_version"] = "1.12.1"
        extra["loupeR_hdf5r_version"] = "1.3.8"
        extra["loupeR_seurat_object_version"] = "5.0.0"
        extra["loupeR_seurat_version"] = "5.0.0"

        meta["extra"] = extra

    return meta

def create_datasets(parent_group, data, groupname):
    group = parent_group.create_group(groupname)

    for name, val in data.items():
        if isinstance(val, dict):
            create_datasets(group, val, name)
        else:
            if isinstance(val, str):
                #dtype = h5py.string_dtype(encoding='utf-8')
                dtype = f'S{len(val)}'
            else:
                dtype = None
            group.create_dataset(name, data=val, shape=(1,), dtype=dtype)

def write_metadata(f, scanpy_obj_version=None):
    metadata = create_metadata(scanpy_obj_version)  # Assume create_metadata function exists

    create_datasets(f, metadata, "metadata")

def create_dataset(obj, key, value, dtype=None):
    """
    not used currently
    """
    obj.create_dataset(key, 
                       data=np.array(value, dtype=dtype), 
                       #data=np.array(strs, dtype=object), 
                       #dtype=dtype
                    )
    obj.close()

def create_str_dataset(obj, key, strs, dtype=None):
    isScalar = isinstance(strs, str)
    compression=None

    if len(strs) == 0:
        max_len = 1
    else:
        max_len = len(strs) if isScalar else max(len(s) for s in strs)
        
    if isScalar:
        shape=(1,)
        data = strs

    else:
        shape=None
        data = np.array(strs, dtype=(f'S{max_len}'))
        
        if len(data) > 2:
            compression='gzip'

    if dtype == None:
        dtype = h5py.string_dtype(encoding='ascii', length=max_len) 
    
    obj.create_dataset(
        key, 
        data=data, 
        dtype=dtype,
        shape=shape,
        compression=compression
        )
def isIntegers(array_like):
    """
    Helper function. Checks if an array contains only integers
    """

    from scipy.sparse import issparse

    elements_to_check = 100000
    if issparse(array_like):
        data = array_like.data[:elements_to_check]
    else:
        rows_to_check = np.ceil(elements_to_check/array_like.shape[1]).astype('int')
        data = array_like[:rows_to_check,:]
    result = np.all(np.mod(data, 1) == 0)
        
    return result

def run_louper(h5_path, cloupe_path, force=False):
    cmd = [louper_path, "create", "--input={}".format(h5_path), "--output={}".format(cloupe_path)]
    if force:
        cmd.append("--force")
    subprocess.run(cmd)

