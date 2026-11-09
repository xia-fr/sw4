# SW4 ESSI converter

`convert_sw4drm.py` extracts DRM boundary motions from an SW4 `ssioutput` (ESSI)
HDF5 file, optionally ZFP-compressed, and writes the HDF5 file that SW4 reads
with its `DRM` command.

It requires the Python packages numpy, scipy, pandas, h5py, mpi4py, hdf5plugin
(for ZFP) and matplotlib.

## Inputs

1. **The ESSI file** (`-e`), written by SW4's `ssioutput` command. It holds the
   velocity time series `vel_0/1/2 ijk layout` on a regular grid.
2. **A DRM configuration CSV** (`-d`) that defines the DRM box, for example:

   ```
   drm_origin_x,drm_origin_y,x_length,y_length,depth,h
   2700,2700,2100,2100,1000,50
   ```

   | Column | Meaning |
   |---|---|
   | `drm_origin_x`, `drm_origin_y` | DRM box origin in the original SW4 coordinates (the box starts at the surface, z = 0) |
   | `x_length`, `y_length`, `depth` | Box size |
   | `h` | DRM node spacing, and the grid spacing of the SW4 DRM run. Equal to the ESSI spacing is ideal; see *Spatial interpolation*. Each box length must be a whole multiple of `h`, otherwise the script stops with an error |

## Command-line options

```
-v or --verbose:   increase output verbosity
-p or --plotonly:  only generate plots of the input nodes
-s or --savepath:  full path for saving the result files, default=""
-d or --drm:       full path to the DRM configuration CSV file, default=""
-e or --essi:      full path to the SW4 ESSI output file (required), default=""
-t or --tstep:     keep every n-th ESSI time step in the output, default=1 (all steps)
-c or --chunktime: time samples per HDF5 chunk of node_motion_data, default=512
-k or --essidumpinterval: legacy ESSI files only, the ssioutput dumpInterval of the source run, default=0 (every step)
```

## What the script does

1. **Generates the DRM nodes.** Four layers of nodes on the four sides and the
   bottom of the box, at spacing `h`.
2. **Checks the nodes** lie inside the ESSI volume, and plots them (rank 0, at
   most about 50,000 nodes, to `input_coordsXYZ/XZ/XY.png`).
3. **Splits the work** over MPI ranks by ESSI chunk, and each rank reads only the
   ESSI chunks its nodes need.
4. **Interpolates in space if needed** (see below), and **downsamples in time**
   if `-t` is greater than 1.
5. **Writes `DRMNodeMotion.h5`** in the save path, one rank at a time.

## Output file

| Dataset | Content |
|---|---|
| `node_motion_data` | velocity, 3N x nt (float32); rows 3i, 3i+1, 3i+2 are ux, uy, uz of node i |
| `xyz` | node coordinates relative to the DRM origin, N x 3 |
| `dt` | time step of the series (ESSI dt x tstep) |
| `h` | DRM node spacing (`h` from the CSV); SW4 builds the DRM grid with it |
| `essi_h` | ESSI grid spacing, for reference only (SW4 does not read it) |
| `t_start`, `t_end` | time range |
| `grid_lon_lat`, `grid_az` | SW4 grid geolocation, copied from the ESSI file |
| `drm_origin`, `drm_shape` | DRM box origin and size |

## Spatial interpolation

DRM nodes that fall between ESSI grid points (for example a DRM `h` finer than
the ESSI spacing, or a DRM origin off the ESSI grid) are interpolated trilinearly
from the 8 surrounding grid points, and the script prints
`Use trilinear interpolation (k of N nodes are off the ESSI grid)`. Otherwise it
prints `No interpolation is needed` and copies each node's series. A node counts
as on the grid when its grid index, relative to the ESSI origin, is within 1e-6
of a whole number.

Interpolation adds no information beyond the ESSI resolution, so a DRM `h` equal
to the ESSI spacing, with the DRM origin on the ESSI grid, is still the most
accurate choice. Halving `h` gives about 4 times as many nodes, and the output
file, memory per rank and run time grow by the same factor. The SW4 DRM run uses
this `h` as its grid spacing, so its CFL time step shrinks with it (see *Using
the file in SW4*).

## Downsampling in time (`-t n`)

The output keeps ESSI samples n-1, 2n-1, 3n-1, ..., and its `dt` is n times the
ESSI sample spacing. These are the samples at t = dt, 2 dt, 3 dt, ..., which is
where SW4 places the output columns (see *Time convention*). Versions before
this fix kept samples 0, n, 2n, ..., which put every boundary motion (n-1) ESSI
steps late; regenerate files made with `-t` > 1. SW4 reads `dt` from the file
and spline-interpolates each node's series onto its own time step, so nothing
else needs changing. No anti-alias filter is
applied: the ESSI motions must have no content above the new Nyquist frequency,
1/(2 n dt_ESSI), which the script prints. A safe choice keeps it at least 2 to 3
times the simulation's maximum resolved frequency. The series are downsampled
right after reading, so memory per rank and the file size drop by about n.

## Time convention

SW4 stores ESSI sample i at t = (i+1) x `timestep` (the first sample is the
first stored step, not t = 0), and its DRM code places output column m at
t = (m+1) x `dt`, after a zero sample at t = 0. The script relies on the ESSI
`time start` being one `timestep`, i.e. a source run that starts at t = 0, and
warns if it is not.

**`ssioutput dumpInterval`.** An SW4 with the ESSI fix writes `timestep` as the
spacing of the stored samples (solver dt x dumpInterval), `time start` as the
time of the first stored sample, and the interval itself as `dump interval`.
The script reads these and needs no option. Older SW4 versions wrote the solver
dt as `timestep` whatever the dumpInterval, and could overwrite the last stored
sample with the solution at the final step. For an ESSI file from an older SW4
whose source run used `dumpInterval=k` with k > 1, pass `-k k`: the script then
scales the sample spacing by k and drops the last sample. Without `-k` such a
file is read as if every step were stored, which it cannot detect. Files from
runs without `dumpInterval` (or with 1) need nothing.

## Chunk layout (`-c`)

`node_motion_data` is written in chunks of 192 rows (64 nodes x 3 components) by
`--chunktime` time samples, both shrunk if the dataset is smaller. SW4 reads the
file one time window at a time (`window=` on the `DRM` line), so keep
`--chunktime` equal to the SW4 window or a divisor of it. A mismatch only costs
I/O speed, never accuracy.

To check or change the layout of an existing file:

```
h5dump -pH DRMNodeMotion.h5 | grep -A3 node_motion_data
h5repack -l node_motion_data:CHUNK=192x512 DRMNodeMotion.h5 DRMNodeMotion_chunked.h5
```

## Using the file in SW4

```
DRM file=/path/to/DRMNodeMotion.h5 window=512 overlap=16
```

`window=` is the number of time intervals each DRM source keeps in memory
(default 512, `0` = the whole series), and `overlap=` is the number of extra
samples fitted on each side of a window (default 16). At startup SW4 prints the
number of DRM nodes each rank owns and the estimated memory on the busiest rank.

**Time step and run length.** SW4 computes the DRM run's time step from the CFL
condition, exactly as for a regular run (grid spacing `h`, material, attenuation),
and interpolates the node motions in time; the file's `dt` is only the sample
spacing. It prints both. The node motions cover t = 0 to nt x `dt`:

- `time t=...` is honoured, but capped at the end of the node motions (with a
  warning). To run the whole record, give the record length or anything larger.
- `time steps=...` is kept as given; SW4 warns if the steps run past the end of
  the node motions, where the boundary motion is held constant.

## Example usage

```
mpirun -np 1 python convert_sw4drm.py -e /resnick/groups/enceladus/fxia/2025_lbnl/06_basin_drm/drm_files/ssi.ssi -d /resnick/groups/enceladus/fxia/2025_lbnl/06_basin_drm/drm_config.csv -s /resnick/groups/enceladus/fxia/2025_lbnl/06_basin_drm/drm_inputs -v
```

On Perlmutter, with the conda environment that has mpi4py:

```
srun -n 77 python convert_sw4drm.py -e $PSCRATCH/run/ssi.ssi -d drm_config.csv -s $PSCRATCH/run/drm -t 1 -c 512 -v
```

The script prints the recommended number of ranks, which is the number of ESSI
chunks the DRM nodes touch.
