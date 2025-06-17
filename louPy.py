import os
import re
import numpy as np
import pandas as pd
import h5py
import platform
import pkg_resources
import subprocess
from scipy.sparse import issparse
import hashlib
import stat
import urllib.request
from datetime import datetime
from pathlib import Path
from platformdirs import user_data_dir

"""
This module facilitates the creation of Loupe files from Scanpy objects using the LoupeR tool (10X Genomics).
It handles the generation of a LoupeR-compatible HDF5 file, invokes the LoupeR executable, and manages file paths.

Setup Instructions:
1. The `setup` function will automatically download the LoupeR executable and ensure the EULA is agreed to.
2. The LoupeR executable will be stored in a user-specific directory unless otherwise specified.
3. Ensure you have an active internet connection for the initial setup.

Scanpy object requirements:
- The Scanpy object (`adata`) should contain the necessary data in its layers, obs, and obsm attributes.
- Integer counts should be present in the specified layer (default is 'counts').

Usage:
    make_loupe(
        adata, 
        cloupe_path, 
        overwrite=False,  # Overwrite an existing .cloupe file with the same name, if True
    )
"""

######## Setup LoupeR executable ########
# artifacts mapping (R version hard-codes v1.1.4); 
# TODO: dynamically fetch latest release
_artifacts = {
    'linux': {
        'url':    "https://github.com/10XGenomics/loupeR/releases/download/v1.1.4/louper-linux-x64",
        'md5':    "b3fd93fd88a43fbcf3f6e40af3186eaa"
    },
    'mac': {
        'url':    "https://github.com/10XGenomics/loupeR/releases/download/v1.1.4/louper-macos-x64",
        'md5':    "ea65a2ec372d623c54d45c51793014e2"
    },
    'windows': {
        'url':    "https://github.com/10XGenomics/loupeR/releases/download/v1.1.4/louper-windows-x64.exe",
        'md5':    "f5d1e99138e840169a19191d10bb25ab"
    },
}

def executable_basename():
    return "louper.exe" if platform.system().lower().startswith("win") else "louper"

def default_executable_path():
    try:    # Check if the executable can be written to the current directory
        d = os.path.dirname(os.path.abspath(__file__))
        if os.access(d, os.W_OK):
            return os.path.join(d, executable_basename()) 
        else:   # Otherwise, use a user-specific data directory
            d = user_data_dir("louPy", "10XGenomics")
            Path(d).mkdir(parents=True, exist_ok=True)
            return os.path.join(d, executable_basename())
    except Exception as e:
        raise RuntimeError(f"Error determining the default executable path: {e}")

def bundled_executable_path():
    """
    Path to the bundled executable, if it exists.
    Currently, the executable is not bundled with the package.
    """
    return os.path.join(os.path.dirname(__file__), "bin", executable_basename())

def get_artifact():
    system = platform.system().lower()
    if system == "darwin":
        key = "mac"
    elif system in ("windows", "cygwin", "msys"):
        key = "windows"
    else:
        key = "linux"
    return _artifacts[key]

def verify_executable(path):
    if not os.path.isfile(path):
        return False, "not found"
    data = open(path, "rb").read()
    if hashlib.md5(data).hexdigest() != get_artifact()["md5"]:
        return False, "md5 mismatch"
    return True, ""

def find_executable():
    for p in (default_executable_path(), bundled_executable_path()):
        ok, _ = verify_executable(p)
        if ok:
            return p
    return None

def install_executable(force=False):
    dest = default_executable_path()
    if os.path.exists(dest) and not force:
        return True
    art = get_artifact()
    req = urllib.request.Request(art["url"])
    token = os.getenv("GITHUB_TOKEN") or os.getenv("GITHUB_PAT")
    if token:
        req.add_header("Authorization", f"token {token}")
        req.add_header("Accept", "application/octet-stream")
    with urllib.request.urlopen(req) as r, open(dest, "wb") as out:
        out.write(r.read())
    os.chmod(dest, stat.S_IRUSR | stat.S_IXUSR)
    ok, msg = verify_executable(dest)
    return ok

# EULA acceptance marker
_eula_file = os.path.join(os.path.dirname(__file__), "eula_accepted")

def eula_have_agreed():
    return os.path.isfile(_eula_file)

def agree_eula():
    prompt = (
        "Please read and accept the LoupeR EULA at\n"
        "https://github.com/10XGenomics/loupeR/blob/main/EULA.md\n"
        "Type 'yes' to accept: "
    )
    if input(prompt).strip().lower() == "yes":
        #Path(_eula_file).write_text(datetime.utcnow().isoformat())
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
        Path(_eula_file).write_text(f'EULA agreed on: {timestamp}')
        return True
    return False

def setup(executable_path=None):
    # install if missing
    if executable_path is None and find_executable() is None:
        install_executable()
    # require EULA
    if not eula_have_agreed():
        agree_eula()

def needs_setup(executable_path=None):
    has_exec = (executable_path or find_executable()) is not None
    has_eula = eula_have_agreed()
    if not (has_exec and has_eula):
        print('has_exec',has_exec)
        print('has_eula',has_eula)
        raise RuntimeError("Call setup() to install LoupeR/LouPy and agree to the EULA")
    return True

#script_dir = os.path.dirname(os.path.abspath(__file__))
#louper_path = os.path.join(script_dir, "louper")

######## LouPy functions ########
def make_loupe(adata, 
               cloupe_path,
               layer='counts',
               clusters = None, 
               overwrite=False, 
               h5path=None,
               ):
    """
    Purpose: This function orchestrates the process of creating a Loupe file from a Scanpy object. 
    It involves generating an HDF5 file, running the LoupeR tool, and handling the file paths.

    Parameters:
        adata: The Scanpy object containing the data to be converted.
        cloupe_path: The path where the output Loupe file should be saved.
        clusters (optional): Specifies which clusters to include from adata.
        overwrite (optional): If True, overwrites any existing Loupe file with the same name.
        h5path (optional): The path for the intermediate HDF5 file.
    """
    # Ensure LoupeR is installed and EULA is agreed to
    setup()         # install & EULA on first run
    needs_setup()   # guard before calling the binary

    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

    if h5path is None:
        h5file = f'tmp_{timestamp}.h5'
        h5path = os.path.join(os.getcwd(), h5file)

    # Initialize a dictionary for logging purposes
    global _log 
    _log = {'n_cells': adata.shape[0], 
            'n_features': adata.shape[1]}

    create_hdf5(adata, 
                h5path,
                layer=layer, 
                clusters=clusters)

    if cloupe_path is None:
        cloupe_file = f'Converted_{timestamp}.cloupe'
        cloupe_path = os.path.join(os.getcwd(), cloupe_file)

    if os.path.exists(cloupe_path) and not overwrite:
        raise FileExistsError(f"Loupe file {cloupe_path} already exists. Use overwrite=True to replace it.")
    
    #print log:
    log_string = f"Writing Loupe file with shape: {adata.shape[0]:,} cells X {adata.shape[1]:,} features"
    print()
    if 'n_clusters' in _log:
        log_string += f" and {_log['n_clusters']} annotations"
    log_string += '...'
    print(log_string)

    run_louper(h5path, cloupe_path, overwrite=overwrite)
    
    # Clean up temporary HDF5 file
    os.remove(h5path)

def create_hdf5(adata, 
                h5path,
                layer='counts', 
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
        write_mat(f, adata, layer=layer)
        write_clusters(f, adata, clusters=clusters)
        write_projections(f, adata)
        metadata = create_metadata()
        # Assuming you have a function to write metadata to the HDF5 file
        write_metadata(f)

    return "SUCCESS"

def run_louper(h5_path, cloupe_path, overwrite=False):
    """
    Purpose: Runs the LoupeR tool to convert an HDF5 file into a Loupe file.

    Parameters:
        h5_path: Path to the input HDF5 file.
        cloupe_path: Path where the output Loupe file should be saved.
        force (optional): If True, overwrites an existing Loupe file.
    """
    louper_path = find_executable()
    cmd = [louper_path, 
           "create", 
           "--input={}".format(h5_path), 
           "--output={}".format(cloupe_path)]
    if overwrite:
        cmd.append("--force")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        # Print path and file size of the generated Loupe file
        file_size = os.path.getsize(cloupe_path) / (1024**2)  # Convert to MB
        print(f"Loupe file successfully created  {cloupe_path}: {file_size:,.2f} MB")
        
    except subprocess.CalledProcessError as e:
        print(f"Error running LoupeR: {e.stderr}")
        raise RuntimeError(f"LoupeR failed with exit code {e.returncode}")


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

    # Prepare data matrix
    if layer == 'X':
        count_mat = adata.X
    elif layer == 'raw':
        if adata.raw is None:
            raise ValueError("adata.raw is None, cannot use layer='raw'")
        count_mat = adata.raw.X
    else:
        if layer not in adata.layers:
            raise ValueError(f"{layer} layer not found in adata.layers")
        count_mat = adata.layers[layer]
    if not issparse(count_mat):
        raise TypeError(f"The provided layer ({layer}) is not a sparse matrix")
    if not hasattr(count_mat, 'shape') or not hasattr(count_mat, '__getitem__'):
        raise TypeError(f"The provided layer ({layer}) is not a valid array-like object")
    if not isIntegers(count_mat):
        raise ValueError(f"The provided layer ({layer}) must contain only integer values")
    
    barcodes_unmodified = adata.obs_names.tolist()
    barcodes_formatted = barcodes_unmodified ## REPLACE THIS WITH SANITIZED BARCODE OUTPUT 
    barcode_count = len(barcodes_unmodified)

    ### Features group
    matrix_group = f.create_group("matrix")

    # Use adata.var_names as features
    features = adata.var_names.to_numpy(dtype=str).astype('S')
    feature_count = len(features)

    # Define feature IDs
    ensembl_pattern = re.compile(r'\bENS[A-Z]{1,5}\d{11}\b')
    def contains_ensembl_ids(column: pd.Series) -> bool:
        """
        Checks if a given AnnData .var columns contains stable Ensembl IDs.
        Allow for up to 5% of the features to be non-ensembl IDs to account for 
        added transgenes or other non-ensembl features.
        """
        # Regular expression to match stable Ensembl IDs (e.g., ENSG00000123456)
        matches = column.astype(str).str.contains(ensembl_pattern)
        ensembl_matches = matches.sum()
        return ensembl_matches / len(column) >= 0.95
    
    ensembl_cols = [name for name in adata.var.columns 
                    if contains_ensembl_ids(adata.var[name])]

    if ensembl_cols:
        if len(ensembl_cols) > 1:
            print(f"Warning: Multiple columns in adata.var match Ensembl ID format: {ensembl_cols}. Using the first one.")
            ensembl_col = ensembl_cols[0]
        else:
            ensembl_col = ensembl_cols[0]
        feature_ids = adata.var[ensembl_col].to_numpy(dtype=str).astype('S')
        _log.update({'ensembl_ids': True})
        _log.update({'ensembl_id_column': ensembl_col})
    else:
        # Otherwise, create generic numeric feature IDs
        _log.update({'ensembl_ids': False})
        _log.update({'ensembl_id_column': 'None_Generic_IDs_Created'})
        feature_ids = [f"feature_{i}" for i in range(1, len(features) + 1)]
    
    # Check for feature_type column (eg. for 10X feature barcoding):
    if 'feature_type' in adata.var.columns:
        feature_types = adata.var['feature_type'].to_numpy(dtype=str).astype('S')
    else:
        feature_types = np.array(['Gene Expression'] * len(features), dtype='S')
        
    # Write Features group
    features_group = matrix_group.create_group("features")
    create_str_dataset(features_group, "name", strs=features)
    create_str_dataset(features_group, "id", strs=feature_ids)
    create_str_dataset(features_group, "feature_type", strs=feature_types)
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
    if clusters is None:
        clusters = categorical_data.columns.tolist()

    if len(clusters) > max_clusterings:
        raise ValueError(f"Number of specified clusters ({len(clusters)}) exceeds the maximum allowed ({max_clusterings}).")
    if len(clusters) > max_clusterings:
        force = True

    if not force:
        if len(categorical_data.columns) > max_clusterings:
            raise ValueError(f"Too many categorical columns in adata.obs. Limit is {max_clusterings}.")
    
    _log.update({'n_clusters': len(clusters)})

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
        
        """
        # Deprecated this.  I want to be able to add in any 2D dimensinoality reduction
        # so I just use the name as the "method" variable
        # This allows for projections named "umap" and "umap_harmony", for example
        
        is_umap = bool(re.search("umap", name, re.IGNORECASE))
        is_tsne = bool(re.search("tsne", name, re.IGNORECASE))
        is_tsne_dash = bool(re.search("t-sne", name, re.IGNORECASE))
        
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
    meta["tool"] = "louPy"
    
    try:
        meta["tool_version"] = pkg_resources.get_distribution("louPy").version
    except pkg_resources.DistributionNotFound:
        meta["tool_version"] = "unknown"
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
