import os
import re
import numpy as np
import h5py
import platform
import pkg_resources
import subprocess
from datetime import datetime
from scipy.sparse import issparse


"""
Useful hack to write a Loupe file from a Scanpy object using LoupeR (10X Genomics)
This module writes a LoupeR-compatible hdf5 file, then passes it to the LoupeR executable

Step 1: Download and install LoupeR from 10X Genomics
https://github.com/10XGenomics/loupeR

Step 2: Run "setup" to download LoupeR executable and agree to EULA

Step 3: Copy LoupeR executable file to the same directory as this .py file

Usage: 
   
    make_loupe(
        adata, 
        cloupe_path, 
        force=False #Overwrite .cloupe file of the same name, if it exists 
        )

"""

script_dir = os.path.dirname(os.path.abspath(__file__))
louper_path = os.path.join(script_dir, "louper")

def make_loupe(adata, 
               cloupe_path,
               clusters = None, 
               force=False, 
               h5path=None,
               ):
    """
    Purpose: This function orchestrates the process of creating a Loupe file from a Scanpy object. 
    It involves generating an HDF5 file, running the LoupeR tool, and handling the file paths.

    Parameters:
        adata: The Scanpy object containing the data to be converted.
        cloupe_path: The path where the output Loupe file should be saved.
        clusters (optional): Specifies which clusters to include from adata.
        force (optional): If True, overwrites any existing Loupe file with the same name.
        h5path (optional): The path for the intermediate HDF5 file.
    """

    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

    if h5path is None:
        h5file = f'tmp_{timestamp}.h5'
        h5path = os.path.join(os.getcwd(), h5file)

    create_hdf5(adata, 
                h5path, 
                clusters=clusters)

    if cloupe_path is None:
        cloupe_file = f'Converted_{timestamp}.cloupe'
        cloupe_path = os.path.join(os.getcwd(), cloupe_file)

    run_louper(h5path, cloupe_path, force=force)

    os.remove(h5path)

def create_hdf5(adata, 
                h5path, 
                clusters=None):
    """
    Purpose: Generates an HDF5 file to be read by LoupeR binary.

    Parameters:
        adata: The Scanpy object to be converted.
        h5path: Path where the HDF5 file should be saved.
        clusters (optional): Specific clusters to include in the HDF5 file.
        
    """
    if os.path.exists(h5path):
        raise ValueError(f"Cannot create h5 file {h5path}, file already exists.")
    
    with h5py.File(h5path, 'w') as f:
        write_mat(f, adata)
        write_clusters(f, adata, clusters=clusters)
        write_projections(f, adata)
        metadata = create_metadata()
        # Assuming you have a function to write metadata to the HDF5 file
        write_metadata(f)

    return "SUCCESS"

def run_louper(h5_path, cloupe_path, force=False):
    """
    Purpose: Runs the LoupeR tool to convert an HDF5 file into a Loupe file.

    Parameters:
        h5_path: Path to the input HDF5 file.
        cloupe_path: Path where the output Loupe file should be saved.
        force (optional): If True, overwrites an existing Loupe file.
    """

    cmd = [louper_path, "create", "--input={}".format(h5_path), "--output={}".format(cloupe_path)]
    if force:
        cmd.append("--force")
    subprocess.run(cmd)


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
        TODO: actually write a full-featured barcode sanitization function
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

def write_clusters(f, 
                   adata,
                   clusters = None,
                   force=False):
    """
    Writes categorical observations in adata.obs to the h5 file.
    LoupeR calls these 'clusters', but they could be any useful categorical variable.
    eg. Sample, Batch, cell cycle phase, or anything else you have tracked.

    To avoid generating an overly cluttered Loupe file, defaults to limiting the number
    of 'cluster' groupings to 32.

    Parameters:
    clusters (Optional): list of column names in adata.var to include as 'cluster' groups
    force: If you really want more than 32 cluster groups, set force=True.
    """
    max_clusterings = 32
    if isinstance(clusters, str):
        clusters = [clusters]

    categorical_data = adata.obs.select_dtypes(include=['category'])

    if clusters:
        categorical_data = categorical_data.loc[:,clusters]
        if len(clusters) > max_clusterings:
            force = True

    if not force:
        if len(categorical_data.columns) > max_clusterings:
            raise ValueError(f"Too many categorical columns in adata.obs. Limit is {max_clusterings}.")

    # Start writing to the HDF5 file
    clusters_group = f.create_group("clusters")
    for name, cluster in categorical_data.items():
        
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
    """
    Purpose: 
    Writes projection data (like UMAP or t-SNE) from a Scanpy object to an HDF5 file.
    Finds any 2D projection in adata.obsm and writes it to the HDF5
    
    Parameters:
        f: The HDF5 file object.
        adata: The Scanpy object.
    """
    projections_group = f.create_group("projections")

    for name, projection in adata.obsm.items():
        name = name.replace('X_','')
        n_dim = projection.shape[1]
        is_umap = bool(re.search("umap", name, re.IGNORECASE))
        is_tsne = bool(re.search("tsne", name, re.IGNORECASE))
        is_tsne_dash = bool(re.search("t-sne", name, re.IGNORECASE))

        """
        # Deprecated this.  I want to be able to add in any 2D dimensinoality reduction
        # so I just use the name as the "method" variable
        # This allows for projections named "umap" and "umap_harmony", for example
        
        if is_umap:
            method = "UMAP"
        elif is_tsne or is_tsne_dash:
            method = "t-SNE"
        else:
            method = name
        """

        if n_dim == 2:
            group = projections_group.create_group(name)
            create_str_dataset(group, "name", strs=name)
            create_str_dataset(group, "method", strs=name)
            group.create_dataset("data", data=projection.T, compression='gzip')
            
def create_metadata():
    """
    Purpose: 
        Creates metadata for the HDF5 file.
        Returns: A dictionary containing metadata.
    """
    
    # Create metadata dictionary
    meta = {}
    meta["tool"] = "loupePy"
    #meta["tool_version"] = pkg_resources.get_distribution("loupePy").version if pkg_resources.get_distribution("loupePy") else "n/a"
    meta["tool_version"] = "pre-alpha"
    meta["os"] = platform.system()
    meta["system"] = platform.platform()
    meta["language"] = "Python"
    meta["language_version"] = platform.python_version()
    
    # Create extra dictionary
    extra = {}
    extra["loupePy_scanpy_version"] = pkg_resources.get_distribution("scanpy").version if pkg_resources.get_distribution("scanpy") else "n/a"
    #extra["loupePy_scanpy_object_version"] = scanpy_obj_version if scanpy_obj_version else "n/a"
    extra["loupePy_scanpy_object_version"] = "n/a"
    extra["loupePy_hdf5_version"] = h5py.version.hdf5_version
    
    # Add extra to meta
    meta["extra"] = extra

    return meta

def create_datasets(parent_group, data, groupname):
    """
    Purpose: 
        Recursively creates datasets within an HDF5 file from a nested dictionary.
    
    Parameters:
        parent_group: The parent group in the HDF5 file.
        data: The dictionary containing data to be written.
        groupname: The name of the group to be created.
    
    Returns: 
        None. The function writes data directly to the HDF5 file.
    """
    group = parent_group.create_group(groupname)

    for name, val in data.items():
        if isinstance(val, dict):
            create_datasets(group, val, name)
        else:
            if isinstance(val, str):
                dtype = f'S{len(val)}'
            else:
                dtype = None
            group.create_dataset(name, data=val, shape=(1,), dtype=dtype)

def write_metadata(f):
    """
    Purpose: 
        Writes metadata to an HDF5 file.

    Parameters:
        f: The HDF5 file object.
    """
    metadata = create_metadata()

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
    Purpose: 
    Checks if a given array contains only integer values.
    
    Parameters:
        array_like: The array to be checked.
    
    Returns: 
        True if all elements are integers, otherwise False.
    """

    elements_to_check = 100000
    if issparse(array_like):
        data = array_like.data[:elements_to_check]
    else:
        rows_to_check = np.ceil(elements_to_check/array_like.shape[1]).astype('int')
        data = array_like[:rows_to_check,:]
    result = np.all(np.mod(data, 1) == 0)
        
    return result
