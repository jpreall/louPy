# louPy
Loupe File Generator for Scanpy

### Overview
This Python script converts a Scanpy AnnData object into a Loupe-compatible HDF5 file and then runs the LoupeR executable to generate a .cloupe file.

### Prerequisites
Python 3.x
h5py
numpy
Scanpy
LoupeR executable

### Installation
Download and install LoupeR from 10X Genomics.
Place the LoupeR executable in the same directory as this script.

### Usage
```
make_loupe(
    adata, 
    cloupe_path, 
    force=False # Overwrite .cloupe file of the same name, if it exists 
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
