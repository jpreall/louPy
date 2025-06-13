# louPy
Loupe File Generator for Scanpy

### Overview
This Python script converts a Scanpy AnnData object into a Loupe-compatible HDF5 file and then runs the LoupeR executable to generate a .cloupe file. For best results, store your original integer UMI counts in the layer `adata.layers['counts']`. 

### Requirements
Python 3.x   
h5py>=3.8.0  
scanpy>=1.5.0  
numpy>=1.20.3  
pandas>=1.0.0  
scipy>=1.5.3  
platformdirs>=2.5.0  

### Installation
Clone this repository. When first running `make_loupe`, setup will prompt you to agree to the 10X Genomics EULA and download the necessary Rust executable to create the Loupe file.

### Usage
```
from louPY impmort make_loupe

make_loupe(
    adata, 
    cloupe_path, 
    )
```

#### Optional parameters:  
`force (Bool)`: If True, overwrite .cloupe file of the same name, if it exists. Default = `False`  
`clusters (list)`: list of categorical columns in `adata.obs` to include as 'cluster' groups in Loupe file.  Defaults to use all available catagorical columns in `adata`, if fewer than 16.


For example:
```
make_loupe(
    adata, 
    cloupe_path = '/path/to/write/myloupe.cloupe',
    clusters = ['Clusters_res0.3', 'Clusters_res0.5','Batch','Donor'], 
    force = True
    )
```

### Functions
`make_loupe()`: Main function to create a .cloupe file.  

*helper functions*  
`create_hdf5()`: Writes the HDF5 file from the AnnData object.  
`run_louper()`: Runs the LoupeR executable to convert the HDF5 to .cloupe.  
`write_mat()`: Writes the data matrix to the HDF5 file.  
`write_clusters()`: Writes cluster information to the HDF5 file.  
`write_projections()`: Writes dimensionality reduction data to the HDF5 file.  
`create_metadata()`: Creates metadata for the HDF5 file.  

### Limitations
The `write_mat()` function currently requires the data matrix to contain only integer values.  
The `write_clusters()` function limits the number of cluster groups to 16 unless force=True.  
