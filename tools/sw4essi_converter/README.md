# SW4 ESSI converter
This code converts the SW4 SSIoutput (w/ ZFP compression) to be used in SW4 with DRM.
It requires scipy, numpy, h5py, mpi4py, hdf5plugin (for ZFP), matplotlib python packages.

Supported command line arguments:
```
-v or --verbose: increase output verbosity
-p or --plotonly: only generate plots of the input nodes
-s or --savepath: full path for saving the result files, default=""
-d or --drm: full path to the OpenSees DRM template file with node coordinates, default=""
-e or --essi: full path to the SW4 ESSI output file (required), default=""
```

Example usage:
```
mpirun -np 1 python convert_sw4drm.py -e /resnick/groups/enceladus/fxia/2025_lbnl/06_basin_drm/drm_files/ssi.ssi -d /resnick/groups/enceladus/fxia/2025_lbnl/06_basin_drm/drm_config.csv -s /resnick/groups/enceladus/fxia/2025_lbnl/06_basin_drm/drm_inputs -v
```

