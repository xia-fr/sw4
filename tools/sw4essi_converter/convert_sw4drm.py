#!/usr/bin/env python3
"""
Modified version of the original convert.py, specifically for use
with SW4 -> SW4 DRM.
"""

# from genericpath import exists
import os
# import sys
import argparse
import h5py
import math
import scipy
from scipy import integrate
import numpy as np
import pandas as pd
import datetime
import time
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.dpi'] = 150
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from mpi4py import MPI
import functools
print = functools.partial(print, flush=True) # Don't buffer print
import hdf5plugin # Use this when SW4 output uses ZFP compression, can be installed with "pip install hdf5plugin"

def plot_cube(save_path, cube_definition, x, y, z, view):
    """
    Plots a 3D cube and the grid points of interest specified by 
    the x, y, z arrays.

    Parameters
    ----------
    save_path : str
        Path to where png images of the generated plots will be saved.
    cube_definition : list of tuple
        A list containing four tuples of (float, float, float).
        Each tuple represents:
            - sw4_y_coord (float)
            - sw4_x_coord (float) 
            - sw4_z_coord (float)
        Combined, the four tuples define the 3D cube in the SW4 coordinate frame.
    x, y, z : (N,) np.ndarray
        Contains the (x,y,z) coordinates for the N grid points of interest to plot.
    view : str
        Defines the desired view for the plot. Options are: 'XYZ', 'XZ', or 'XY'.
    """
    cube_definition_array = [
        np.array(list(item))
        for item in cube_definition
    ]

    points = []
    points += cube_definition_array
    vectors = [
        cube_definition_array[1] - cube_definition_array[0],
        cube_definition_array[2] - cube_definition_array[0],
        cube_definition_array[3] - cube_definition_array[0]
    ]

    points += [cube_definition_array[0] + vectors[0] + vectors[1]]
    points += [cube_definition_array[0] + vectors[0] + vectors[2]]
    points += [cube_definition_array[0] + vectors[1] + vectors[2]]
    points += [cube_definition_array[0] + vectors[0] + vectors[1] + vectors[2]]

    points = np.array(points)

    edges = [
        [points[0], points[3], points[5], points[1]],
        [points[1], points[5], points[7], points[4]],
        [points[4], points[2], points[6], points[7]],
        [points[2], points[6], points[3], points[0]],
        [points[0], points[2], points[4], points[1]],
        [points[3], points[6], points[7], points[5]]
    ]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    faces = Poly3DCollection(edges, linewidths=1, edgecolors='xkcd:grey')
    faces.set_facecolor((0,0,1,0.05))
    ax.add_collection3d(faces)

    # Plot the points themselves to force the scaling of the axes
    ax.scatter(points[:,0], points[:,1], points[:,2], s=0)
    #[(1500, 500, 0), (1700, 500, 0), (1500, 650, 0), (1500, 500, 200)]
    x_min = cube_definition[0][0]
    x_max = cube_definition[1][0]
    y_min = cube_definition[0][1]
    y_max = cube_definition[2][1]
    z_min = cube_definition[0][2]
    z_max = cube_definition[3][2]
    x_len = x_max - x_min
    y_len = y_max - y_min
    z_len = z_max - z_min

    x_plot_min = x_min
    x_plot_max = x_max
    y_plot_min = y_min
    y_plot_max = y_max
    z_plot_min = z_min
    z_plot_max = z_max
    
    #print('plot min/max %.2f %.2f %.2f %.2f %.2f %.2f' %( x_plot_min, x_plot_max, y_plot_min, y_plot_max, z_plot_min, z_plot_max))

    x_y_len_diff = abs(x_len-y_len)
    if x_len < y_len:
        x_plot_min = x_min - x_y_len_diff/2
        x_plot_max = x_max + x_y_len_diff/2
    elif x_len > y_len:
        y_plot_min = y_min - x_y_len_diff/2
        y_plot_max = y_max + x_y_len_diff/2
    else:
        tmp0 = 0.95
        tmp1 = 1+(1-tmp0)
        x_plot_min *= tmp0
        x_plot_max *= tmp1
        y_plot_min *= tmp0
        y_plot_max *= tmp1
        
    #print('plot min/max', x_plot_min, x_plot_max, y_plot_min, y_plot_max, z_plot_min, z_plot_max)
       
    ax.set_xlabel('Y(SW4)')
    ax.set_ylabel('X(SW4)')
    ax.set_zlabel('Z(SW4)')
    ax.set_xlim(x_plot_min, x_plot_max) 
    ax.set_ylim(y_plot_min, y_plot_max) 
    ax.set_zlim(z_plot_max, 0) 
    
    lblsize = 5
    ax.zaxis.set_tick_params(labelsize=lblsize)
    ax.yaxis.set_tick_params(labelsize=lblsize)
    ax.xaxis.set_tick_params(labelsize=lblsize)
    
    # ax.dist = 12
    #ax.set_aspect('equal')
    ax.set_box_aspect(None, zoom=0.5)
    
    ax.text(cube_definition[2][0], cube_definition[2][1], cube_definition[2][2]-z_len*.05, 'SW4-ESSI domain', fontsize=7)

    xcolor = 'xkcd:azure'
    ycolor = 'xkcd:green'
    zcolor = 'xkcd:goldenrod'
    xyzmarker = 'x'
    xyzalpha = 0.1
    markersize=2

    # xs = x + cube_definition[0][1]
    # ys = y + cube_definition[0][0]
    # zs = z + cube_definition[0][2]
    xs = x
    ys = y
    zs = z
    
    #print(xs)
    #print(ys)
    #print(zs)
    
    ax.scatter(ys, xs, zs, c='r', marker='.')
    ax.plot(ys, zs, linestyle = 'None', marker=xyzmarker, markersize=markersize, color=ycolor, zdir='y', zs=y_plot_max, alpha=xyzalpha)
    ax.plot(xs, zs, linestyle = 'None', marker=xyzmarker, markersize=markersize,  color=xcolor, zdir='x', zs=x_plot_min, alpha=xyzalpha)
    ax.plot(ys, xs, linestyle = 'None', marker=xyzmarker, markersize=markersize,  color=zcolor, zdir='z', zs=z_plot_max, alpha=xyzalpha)
    
    if view == 'XZ':
        ax.view_init(azim=0, elev=0)    # XZ
        ax.set_proj_type('ortho')
    elif view == 'XY':
        ax.view_init(azim=0, elev=90)   # XY
        ax.set_proj_type('ortho')
    #ax.view_init(azim=0, elev=-90)   # XZ
    fname = save_path + '/input_coords' + view + '.png'
    plt.savefig(fname)

def plot_coords(
    sbdmn_x0, sbdmn_y0, sbdmn_z0, sbdmn_h, 
    sbdmn_nx, sbdmn_ny, sbdmn_nz, 
    user_sbdmn_x, user_sbdmn_y, user_sbdmn_z, 
    save_path='./'
) -> None:
    """
    Plot user specified grid points along with the DRM domain, 
    and its relative location in the SW4 domain.

    Parameters
    ----------
    sbdmn_x0, sbdmn_y0, sbdmn_z0 : float
        Coordinates that define the origin of the DRM subdomain 
        relative to the SW4 origin.
    sbdmn_h : float
        Grid spacing of the DRM subdomain.
    sbdmn_nx, sbdmn_ny, sbdmn_nz : float
        Defines the number of grid points in the DRM subdomain in the
        x, y, and z directions.
    user_sbdmn_x, user_sbdmn_y, user_sbdmn_z : (N,) np.ndarray
        Contains the (x,y,z) coordinates for the N grid points of interest to plot.
    save_path : str
        The path to save the generated png plots to.
    """
    sw4_start_x = sbdmn_x0
    sw4_end_x   = sbdmn_x0 + (sbdmn_nx-1)*sbdmn_h
    sw4_start_y = sbdmn_y0
    sw4_end_y   = sbdmn_y0 + (sbdmn_ny-1)*sbdmn_h
    sw4_start_z = sbdmn_z0
    sw4_end_z   = sbdmn_z0 + (sbdmn_nz-1)*sbdmn_h

    cube_definition = [ (sw4_start_y,sw4_start_x,sw4_start_z), 
                        (sw4_end_y,sw4_start_x,sw4_start_z), 
                        (sw4_start_y,sw4_end_x,sw4_start_z), 
                        (sw4_start_y,sw4_start_x,sw4_end_z)   ]

    # print(cube_definition)
    plot_cube(save_path, cube_definition, user_sbdmn_x, user_sbdmn_y, user_sbdmn_z, 'XYZ')
    plot_cube(save_path, cube_definition, user_sbdmn_x, user_sbdmn_y, user_sbdmn_z, 'XZ')
    plot_cube(save_path, cube_definition, user_sbdmn_x, user_sbdmn_y, user_sbdmn_z, 'XY')

def convert_to_essi_coord(from_x, from_y, from_z, ref_essi_xyz, ref_sbdmn):
    """
    Converts coordinates from the full SW4 domain coordinate frame
    to the SW4 ssioutput ESSI subdomain coordinate frame (not always the same).
    Carried out effectively as a coordinate offset.
    
    Parameters
    ----------
    from_x, from_y, from_z : np.ndarray
        The (x,y,z) coordinate arrays in the original SW4 coordinate frame to convert from.
    ref_essi_xyz : (3,) np.ndarray
        Contains the coordinates of the ESSI subdomain origin relative to the SW4 coord. frame.
    
    Returns
    -------
    Arrays of (x,y,z) grid point coordinates relative to the ESSI subdomain coordinate frame,
    rounded to 8 decimal points.
    """
    user_essi_x = from_x + ref_sbdmn[0] - ref_essi_xyz[0]
    user_essi_y = from_y + ref_sbdmn[1] - ref_essi_xyz[1]
    user_essi_z = from_z + ref_sbdmn[2] - ref_essi_xyz[2]
    
    return np.round(user_essi_x, decimals=8), np.round(user_essi_y, decimals=8), np.round(user_essi_z, decimals=8)

def get_csv_meta(csv_fname):
    """
    Gets the parameter values from the CSV settings file.

    Parameters
    ----------
    csv_fname : str
        The name/path of the csv file to read the parameters from.

    Returns
    -------
    ref_coord : (3,) np.ndarray
        The subdomain origin in the SW4 coordinate frame.
    start_t : float
        Start time.
    end_t : float
        End time.
    tstep : float
        Time step.
    """
    df = pd.read_csv(csv_fname)

    # reference point, which is the ESSI or OpenSees origin in the SW4 coordinate system
    ref_coord = np.zeros(3)
    ref_coord[0] = df['sbdmnXstart'][0]
    ref_coord[1] = df['sbdmnYstart'][0]
    ref_coord[2] = df['sbdmnZstart'][0]

    # start time and end time for truncation
    start_t = df['startTime'][0]
    end_t = df['endTime'][0]
    tstep = int(df['tstep'][0])

    # print('In csv file: ref_coord, start_t, end_t, rotate_angle:', ref_coord, start_t, end_t, rotate_angle)
    return ref_coord, start_t, end_t, tstep

def get_essi_meta(ssi_fname, verbose):
    """
    Get parameter values from the SW4 "ssioutput" HDF5 data file.

    Returns
    -------
    x0, y0, z0 : float
        Origin of the extracted volumetric data region relative to the
        original full SW4 domain.
    h : float
        Grid spacing of the volumetric data.
    nx, ny, nz : int
        Number of grid points in the x, y, z direction for the extracted data.
    nt : int
        Number of time steps.
    dt : float
        Timestep.
    timeseq : np.ndarray
        The sequence of time from the simulation.
    """
    essiout = h5py.File(ssi_fname, 'r')
    h  = essiout['ESSI xyz grid spacing'][0]
    x0 = essiout['ESSI xyz origin'][0]
    y0 = essiout['ESSI xyz origin'][1]
    z0 = essiout['ESSI xyz origin'][2]
    t0 = essiout['time start'][0]
    dt = essiout['timestep'][0]
    nt = essiout['vel_0 ijk layout'].shape[0]
    nx = essiout['vel_0 ijk layout'].shape[1]
    ny = essiout['vel_0 ijk layout'].shape[2]
    nz = essiout['vel_0 ijk layout'].shape[3]
    ll = essiout['Grid lon-lat origin'][:]
    az = essiout['Grid azimuth'][:]
    t1 = t0 + dt*(nt-1)
    timeseq = np.linspace(t0, t1, nt)
    # print('dt, t0, t1, timeseq =', dt, t0, t1, timeseq)
    essiout.close()
    
    return x0, y0, z0, h, nx, ny, nz, nt, dt, timeseq, ll, az

def write_to_hdf5_range(h5_fname, gname, dname, data, mystart, myend):
    """
    Writes data to the hdf5 file within a given range.

    Parameters
    ----------
    h5_fname : str
        The filename for the hdf5 file to write into.
    gname : str
        The HDF5 file group to write the data into. If
        writing into the base file, gname = '/'.
    dname : str
        Name of the data being written into the file.
    data : np.ndarray
        The data to be written.
    mystart, myend : int
        The indices of the dataset within the hdf5 file to
        which the data will be written.
    """
    h5file = h5py.File(h5_fname, 'r+')
    if gname == '/':
        dset = h5file[dname]
    else:
        grp = h5file[gname]
        dset = grp[dname]

    #print('write_to_hdf5_range, data shape:', data.shape, 'dset shape:', dset.shape)
    #print('mystart=%d, myend=%d' %(mystart, myend))
    dset[mystart:myend, :] = data[:]
    h5file.close()

def write_to_hdf5_range_2d(h5_fname, gname, dname, data, mystart, myend):
    """
    Writes 2D data to the hdf5 file within a given range.
    """
    h5file = h5py.File(h5_fname, 'r+')
    if gname == '/':
        dset = h5file[dname]
    else:
        grp = h5file[gname]
        dset = grp[dname]

    #print('mystart=%d, myend=%d' %(mystart, myend))
    dset[mystart:myend,:] = data[:,:]
    h5file.close()

def create_hdf5_csv(
    h5_fname, ncoord, tsteprange, 
    essi_dt, essi_h, 
    gen_vel, gen_acc, gen_dis,
    lonlat, az, sbdmn_0, sbdmn_len
) -> None:
    """
    Sets up the HDF5 file to save node motions into.
    """
    h5file = h5py.File(h5_fname, 'w')
    tstart, tend, tstep = (tsteprange[0]+1)*essi_dt, (tsteprange[-1]+1)*essi_dt, tsteprange.step
    nstep = len(tsteprange)
    dt = tstep*essi_dt
    
    print('Create HDF5 file with ', ncoord, ' coordinates, ', nstep, ' steps')
    dset = h5file.create_dataset('node_motion_data', (ncoord*3, nstep), dtype='f4')
    dset = h5file.create_dataset('xyz', (ncoord, 3), dtype='f4')
    
    # use array form [dt] instead of scalar dt, so that it has a shape and facilitate postprocessing
    dset = h5file.create_dataset('dt', data=[dt], dtype='f8')
    dset = h5file.create_dataset('h', data=[essi_h], dtype='f8')

    dset = h5file.create_dataset('t_start', data=[dt], dtype='f8')
    dset = h5file.create_dataset('t_end', data=[tend - tstart + dt], dtype='f8')
    dset = h5file.create_dataset('grid_lon_lat', data=lonlat, dtype='f8')
    dset = h5file.create_dataset('grid_az', data=az, dtype='f8')
    dset = h5file.create_dataset('drm_origin', data=sbdmn_0, dtype='f8')
    dset = h5file.create_dataset('drm_shape', data=sbdmn_len, dtype='f8')
    
    h5file.close()

def coord_to_chunkid(x, y, z, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z):
    val = int(np.floor(x/chk_x)*nchk_y*nchk_z + np.floor(y/chk_y)*nchk_z + np.floor(z/chk_z))
    #print('coord_to_chunkid:', x, y, z, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z, val)
    return val

def chunkid_to_start(cid, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z):
    #print('cid2:', cid, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z)
    x = math.floor(cid / (nchk_z * nchk_y))
    y = math.floor((cid - x*nchk_z*nchk_y) / nchk_z)
    z = cid - y*nchk_z - x*nchk_z*nchk_y
    return int(x*chk_x), int(y*chk_y), int(z*chk_z)

def get_chunk_size(ssi_fname):
    fid = h5py.File(ssi_fname, 'r')
    dims = fid['vel_0 ijk layout'].chunks
    if not dims:
        dims = fid['vel_0 ijk layout'].shape
    fid.close()
    #print('Chunk size:', dims)
    return int(dims[0]), int(dims[1]), int(dims[2]), int(dims[3])

def get_nchunk_from_coords(x, y, z, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z):
    if len(x) != len(y) or len(y) != len(z):
        print('Not equal sizes of the x,y,z coordinates array')
    chk_ids = {}
    cnt = 0
    for i in range(0, len(x)):
        cid = coord_to_chunkid(x[i], y[i], z[i], chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z)
        if cid not in chk_ids:
            chk_ids[cid] = cnt
            cnt += 1
    return len(chk_ids), chk_ids

def coord_to_str_3d(x, y, z):
    return str(x)+','+str(y)+','+str(z)

def str_to_coord_3d(s):
    val = s.split(',')
    return int(val[0]), int(val[1]), int(val[2])
    
def allocate_neighbor_coords_8(data_dict, x, y, z, n, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z):
    nadd = 0
    add_cids_dict = {}
    neighbour = 2
    for i0 in range(0, neighbour):
        for i1 in range(0, neighbour):
            for i2 in range(0, neighbour):
                intx, inty, intz = int(x+i0), int(y+i1), int(z+i2)
                coord_str = coord_to_str_3d(intx, inty, intz)
                if coord_str not in data_dict:
                    data_dict[coord_str] = np.zeros(n)
                    nadd += 1

                    cid = coord_to_chunkid(intx, inty, intz, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z)
                    if cid in add_cids_dict:
                        add_cids_dict[cid].add(coord_str)
                    else:
                        add_cids_dict[cid] = {coord_str}

                    #print(coord_str)
                #else:
                #    print(coord_str, 'alread in dict')
                
    return nadd, add_cids_dict

def read_hdf5_by_chunk(ssi_fname, data_dict, comp, cids_dict, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z, chk_t, mpi_rank, verbose):
    fid = h5py.File(ssi_fname, 'r')
    dset_name = 'vel_' + str(int(comp)) + ' ijk layout'
    for cids_iter in cids_dict:
        # Read chunk
        nread = math.ceil(fid[dset_name].shape[0] / chk_t)
        for start_t in range(0, nread): 
            start_x, start_y, start_z = chunkid_to_start(cids_iter, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z)
            # print('Read chunk cid =', cids_iter, start_x, chk_x, start_y, chk_y, start_z, chk_z)
            starttime = time.time()
            chk_data = fid[dset_name][int(chk_t*start_t):int(chk_t*(start_t+1)), int(start_x):int(start_x+chk_x), int(start_y):int(start_y+chk_y), int(start_z):int(start_z+chk_z)]
            endtime = time.time()
            if verbose: 
                # print('Rank', mpi_rank, 'read: cid', cids_iter, dset_name, ',time sub chunk', start_t+1, '/', nread, 'time:', endtime-starttime)
                print('Rank {:3d} read: chunk cid {:4d} {}, time slice {:3d}/{} took {:.3f}s'.format(mpi_rank, cids_iter, dset_name, start_t+1, nread, endtime-starttime))
                #sys.stdout.flush()

            starttime = time.time()
            for coord_str in cids_dict[cids_iter]:
                x, y, z = str_to_coord_3d(coord_str)
                # assign values from chunk to data_dict[coord_str][0:3]
                # print('==assign values for', x, y, z, '->', x%chk_x, y%chk_y, z%chk_z, 'cid', cids_iter, 'is in ', cids_iter, 'timestep', chk_t*start_t)
                # print('shape is:', data_dict[coord_str].shape, chk_data.shape)
                data_dict[coord_str][chk_t*start_t:chk_t*(start_t+1)] = chk_data[:,x%chk_x,y%chk_y,z%chk_z]
            endtime = time.time()
            #print('assign value time', endtime-starttime)
    fid.close()

def generate_acc_dis_time(ssi_fname,
    subdomain_origin, subdomain_xyz_len,
    user_x0, user_y0, user_z0, n_coord,
    start_t, end_t, tstep,
    gen_vel, gen_acc, gen_dis,
    verbose, plot_only, output_fname,
    mpi_rank, mpi_size, output_format,
) -> None:
    """
    The function which does the extraction of DRM boundary data and
    saves the extracted data into a hdf5 file for use with DRM.

    Parameters
    ----------
    ssi_fname : str
        Filename of SW4 ssioutput file to extract data from.
    subdomain_origin : list of float
        The origin for the subdomain reference frame that we want to extract
        motions from. This origin coordinate is in the coordinate frame of
        the original SW4 simulation domain.
    user_x0, user_y0, user_z0 : np.ndarray
        Arrays containing the node coords (in the subdomain ref. frame) that
        we want to extract motions from.
    """
    # Read ESSI metadata
    essi_x0, essi_y0, essi_z0, essi_h,  essi_nx, essi_ny, essi_nz,\
        essi_nt, essi_dt, essi_timeseq, essi_lonlat, essi_az = get_essi_meta(ssi_fname, verbose)
    essi_x_len_max = (essi_nx-1) * essi_h
    essi_y_len_max = (essi_ny-1) * essi_h
    essi_z_len_max = (essi_nz-1) * essi_h
    
    # Start and end time step
    if start_t > -1e-6 and end_t > -1e-6:
        start_ts = int(abs(start_t)/essi_dt)
        end_ts   = int(abs(end_t)/essi_dt)
        # if start and end time step equals, we are likely to want all following steps till the end 
        if end_ts > essi_nt or start_ts == end_ts:
            end_ts = int(essi_nt)
        if end_ts <= start_ts:
            print('End time step {} <= start time step {}, no need to extract motions, exit...'.format(end_ts, start_ts))
            exit(0)
    else:
        print('Error getting start and end time step: start_t, end_t, essi_dt =', start_t, end_t, essi_dt)
        exit(0)

    # Time step information
    tsteprange = range(start_ts, end_ts, tstep)
    nsteps = len(tsteprange)
    dt = tstep * essi_dt

    # Print metadata for ssioutput file
    if verbose and mpi_rank == 0:
        print('\nESSI origin x0, y0, z0, h: ', essi_x0, essi_y0, essi_z0, essi_h)
        print('ESSI origin nx, ny, nz, nt, dt: ', essi_nx, essi_ny, essi_nz, essi_nt, essi_dt)
        print('ESSI max len x, y, z: ', essi_x_len_max, essi_y_len_max, essi_z_len_max)
        print('ESSI max x, y, z: ', essi_x0+essi_x_len_max, essi_y0+essi_y_len_max, essi_z0+essi_z_len_max)
        print(' ')
        print('Generate output file with timesteps between', start_ts, 'and', end_ts, 'with step interval', tstep, 'in', output_format, 'format')

    # Convert user coordinate to sw4 coordinate, relative to ESSI domain (potentially a subset of SW4 domain)
    # user_x,user_y,user_z := coordinates in ESSI subdomain
    # user_x0,user_y0,user_z0 := coordinates in DRM boundary subdomain
    ref_coord = np.array([essi_x0, essi_y0, essi_z0])
    user_x, user_y, user_z = user_x0, user_y0, user_z0
    user_essi_x, user_essi_y, user_essi_z = convert_to_essi_coord(user_x, user_y, user_z, ref_coord, subdomain_origin)

    # debug print
    nprint = 0
    for i in range(0, nprint):
        if i == 0:
            print('converted essi coordinate:')
        print('(%d, %d, %d)' % (user_essi_x[i], user_essi_y[i], user_essi_z[i]))

    # Plot
    if mpi_rank == 0:
        plot_coords(essi_x0, essi_y0, essi_z0, essi_h, essi_nx, essi_ny, essi_nz, user_essi_x, user_essi_y, user_essi_z, save_path)

    if plot_only:
        if mpi_rank == 0:
            print('Only generate the plots of input nodes')
        exit(0)

    # Check if all node coordinates are within the sw4 domain
    if np.min(user_essi_x) < essi_x0 or np.max(user_essi_x) > essi_x0+essi_x_len_max or \
       np.min(user_essi_y) < essi_y0 or np.max(user_essi_y) > essi_y0+essi_y_len_max or \
       np.min(user_essi_z) < essi_z0 or np.max(user_essi_z) > essi_z0+essi_z_len_max:
        if mpi_rank == 0:
            print('Error: all node coordinates should be within the sw4 domain for extracting the motion')
            print('while:')
            print('\t','Min/Max SW4 x:',essi_x0,essi_x0+essi_x_len_max,'Min/Max user x:',np.min(user_essi_x),np.max(user_essi_x))
            print('\t','Min/Max SW4 y:',essi_y0,essi_y0+essi_y_len_max,'Min/Max user y:',np.min(user_essi_y),np.max(user_essi_y))
            print('\t','Min/Max SW4 z:',essi_z0,essi_z0+essi_z_len_max,'Min/Max user z:',np.min(user_essi_z),np.max(user_essi_z))
            
            debugfile = save_path + '/user_essi_xyz.npy'
            print('\tcheck user_essi_xyz in file \'{}\''.format(debugfile))
            np.save(debugfile, np.c_[user_essi_x, user_essi_y, user_essi_z])
        exit(0)
    
    # Convert to array location (spacing is 1), floating-point
    coord_x = (user_essi_x - essi_x0) / essi_h
    coord_y = (user_essi_y - essi_y0) / essi_h
    coord_z = (user_essi_z - essi_z0) / essi_h  
    
    # Check if we actually need interpolation
    # ghost_cell = 0
    # do_interp = True
    do_interp = False
    for nid in range(0, n_coord):
        if user_essi_x[nid] % essi_h != 0 or user_essi_y[nid] % essi_h != 0 or user_essi_z[nid] % essi_h != 0:
            do_interp = True
            # ghost_cell = 1
            break    
    if mpi_rank == 0:
      if do_interp:
        print('Use spline interpolation.')
      else:
        print('No spline interpolation is needed.')
    
    #for i in range(0, len(user_essi_x)):
    #    print('(%.2f, %.2f, %.2f)' % (coord_x[i], coord_y[i], coord_z[i]))
    
    chk_t, chk_x, chk_y, chk_z = get_chunk_size(ssi_fname)
    if chk_t <= 0 or chk_x <= 0 or chk_y <= 0 or chk_z <= 0:
        print('Error getting chunk size from essi file', chk_t, chk_x, chk_y, chk_z)
        exit(0)
    
    nchk_x = int(np.ceil(essi_nx/chk_x))
    nchk_y = int(np.ceil(essi_ny/chk_y))
    nchk_z = int(np.ceil(essi_nz/chk_z))
    if nchk_x <= 0 or nchk_y <= 0 or nchk_z <= 0:
        print('Error getting number of chunks', nchk_x, nchk_y, nchk_z)
        exit(0) 
    
    if verbose and mpi_rank == 0:
        print('Essi file: chk_t, chk_x, chk_y, chk_z =', chk_t, chk_x, chk_y, chk_z, ', nchk_x, nchk_y, nchk_z =', nchk_x, nchk_y, nchk_z)
    
    ntry = 0
    ntry_max = 1
    nchk = 0
    # Try to reduce the chunk size if the number of chunks is less than half the number of ranks
    while(nchk < 0.5*mpi_size):
        if ntry > 0:
            if ntry % 3 == 1 and chk_x % 2 == 0:
                # chk_x /= 2
                chk_x = int(chk_x/2)
            elif ntry % 3 == 2 and chk_y % 2 == 0:
                # chk_y /= 2     
                chk_y = int(chk_y/2)    
            elif ntry % 3 == 0 and chk_z % 2 == 0:
                # chk_z /= 2  
                chk_z = int(chk_z/2)
        
        # Find chunks where all the user input coordinates are (not including adjacent chunks for interpolation yet)
        # cids_dict format: {cid1:index1_in_all_cids,}
        nchk, cids_dict = get_nchunk_from_coords(coord_x, coord_y, coord_z, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z)
        # if mpi_rank == 0:
        #   print('ntry, nchk, mpi_size, cids_dict, chk_x, chk_y, chk_z = ', ntry, nchk, mpi_size, cids_dict, chk_x, chk_y, chk_z)
        
        if ntry == 0 and mpi_rank == 0 and nchk != mpi_size:
            print('\nRecommend using', nchk, 'MPI rank(s)', 'instead of currently used', mpi_size, '\n')
        
        # Don't try too many times
        ntry += 1
        if ntry > ntry_max:
            break
    
    if verbose and mpi_rank == 0:
        print(nchk, 'total chunks to read/distribute', 'using chunk size (', chk_x, chk_y, chk_z, ')')
        print('All needed chuck ids and their order: cids_dict =', cids_dict)
    
    # Get the coordinates assigned to this rank
    read_coords_vel_0 = {}
    read_coords_vel_1 = {}
    read_coords_vel_2 = {}
    # coords_str_dict = {}
    my_ncoord = np.zeros(1, dtype='int')
    my_user_coordinates = np.zeros((n_coord,3), dtype='f8')
    my_converted_coordinates = np.zeros((n_coord,3), dtype='f8')
    my_cids_dict = {} # format: {cid1:{coord_str1,},}, includes all the chunks for interpolation in this rank
    
    # Iterate over all of the coordinates
    for i in range(0, n_coord):
        cid = coord_to_chunkid(coord_x[i], coord_y[i], coord_z[i], chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z)
        if cid < 0:
            print('Error with coord_to_chunkid', coord_x[i], coord_y[i], coord_z[i], cid)
            exit(0)
        # # Debug
        # if mpi_rank == 0:
        #     tmp0, tmp1, tmp2 = chunkid_to_start(cid, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z)
        #     print('cid', cid, coord_x[i], coord_y[i], coord_z[i], 'reverse:', tmp0, tmp1, tmp2)
            
        # cids_dict stores the actual unique ids of chunks that contain input coordinates
        if cids_dict[cid] % mpi_size == mpi_rank:
            # if verbose:
            #    print(i, coord_x[i], coord_y[i], coord_z[i], 'goes to chunk', cid, 'and rank', mpi_rank)
            my_user_coordinates[my_ncoord[0], 0] = user_x0[i]
            my_user_coordinates[my_ncoord[0], 1] = user_y0[i]
            my_user_coordinates[my_ncoord[0], 2] = user_z0[i]
            
            my_converted_coordinates[my_ncoord[0], 0] = coord_x[i]
            my_converted_coordinates[my_ncoord[0], 1] = coord_y[i]
            my_converted_coordinates[my_ncoord[0], 2] = coord_z[i]
            
            coord_str = coord_to_str_3d(int(coord_x[i]), int(coord_y[i]), int(coord_z[i]))
            # coords_str_dict[coord_str] = 1
                    
            # if coord_x[i] % 1 != 0 or coord_y[i] % 1 != 0 or coord_z[i] % 1 != 0:
            if do_interp:
                # Linear interpolation requires 8 neighbours' data
                nadded, add_cids_dict = allocate_neighbor_coords_8(read_coords_vel_0, coord_x[i], coord_y[i], coord_z[i], essi_nt, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z)
                nadded, add_cids_dict = allocate_neighbor_coords_8(read_coords_vel_1, coord_x[i], coord_y[i], coord_z[i], essi_nt, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z)
                nadded, add_cids_dict = allocate_neighbor_coords_8(read_coords_vel_2, coord_x[i], coord_y[i], coord_z[i], essi_nt, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z)

                # print('Rank', mpi_rank, ': add_cids_dict =', add_cids_dict)

                for iadd in add_cids_dict:
                    if iadd in my_cids_dict:
                        my_cids_dict[iadd] |= add_cids_dict[iadd]
                    else:
                        my_cids_dict[iadd] = add_cids_dict[iadd]

                #print(int(coord_x[i]), int(coord_y[i]), int(coord_z[i]), 'added', nadded, 'nodes /', len(read_coords_vel_0))
            else:
                if coord_str not in read_coords_vel_0:
                    read_coords_vel_0[coord_str] = np.zeros(essi_nt)
                    read_coords_vel_1[coord_str] = np.zeros(essi_nt)
                    read_coords_vel_2[coord_str] = np.zeros(essi_nt)

                if cid in my_cids_dict:
                    my_cids_dict[cid].add(coord_str)
                else:
                    my_cids_dict[cid] = {coord_str}
                    
            my_ncoord[0] += 1        
        #end if assigned to my rank
    #end for i in all coordinates

    if verbose:
        print('Rank', mpi_rank, 'has my_cids_dict.keys() =', my_cids_dict.keys())

    # Allocated more than needed previously, adjust
    my_user_coordinates.resize(my_ncoord[0], 3)
    my_converted_coordinates.resize(my_ncoord[0], 3)
    
    # if mpi_rank == 0:
    #     # print('read_coords_vel_0 =', read_coords_vel_0)
    #     print('Rank', mpi_rank, ': my_converted_coordinates =', my_converted_coordinates)

    comm = MPI.COMM_WORLD
    all_ncoord = np.empty(mpi_size, dtype='int')
    comm.Allgather([my_ncoord, MPI.INT], [all_ncoord, MPI.INT])
    
    my_nchk = len(my_cids_dict)
    if verbose:
        print('Rank', mpi_rank, ': assigned', my_ncoord, 'nodes, need to read', len(read_coords_vel_0), 'nodes, in', my_nchk, 'chunk')

    if my_ncoord[0] > 0:
        # Read data by chunk and assign to read_coords_vel_012
        read_hdf5_by_chunk(ssi_fname, read_coords_vel_0, 0, my_cids_dict, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z, chk_t, mpi_rank, verbose)
        read_hdf5_by_chunk(ssi_fname, read_coords_vel_1, 1, my_cids_dict, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z, chk_t, mpi_rank, verbose)
        read_hdf5_by_chunk(ssi_fname, read_coords_vel_2, 2, my_cids_dict, chk_x, chk_y, chk_z, nchk_x, nchk_y, nchk_z, chk_t, mpi_rank, verbose)

    # Calculate the offset from the global array
    my_offset = 0
    for i in range(0, mpi_rank):
        my_offset += all_ncoord[i]

    output_acc_all = np.zeros((my_ncoord[0]*3, essi_nt), dtype='f4')
    output_dis_all = np.zeros((my_ncoord[0]*3, essi_nt), dtype='f4')   
    output_vel_all = np.zeros((my_ncoord[0]*3, essi_nt), dtype='f4')    
    
    # Iterate over all coordinates, all the vel data (vel_0 to 2) in read_coords_vel_012 dict for this rank
    # no interpolation needed, just go through all coordinates' data and convert to acc and dis
    iter_count = 0 
    print('Rank', mpi_rank, 'size of read_coords_vel_0:', len(read_coords_vel_0))
    for iter_count in range(0, my_ncoord[0]):
        x = my_converted_coordinates[iter_count, 0]
        y = my_converted_coordinates[iter_count, 1]
        z = my_converted_coordinates[iter_count, 2]
        coord_str = coord_to_str_3d(int(x), int(y), int(z))

        if gen_vel:
            output_vel_all[iter_count*3+0, :] = read_coords_vel_0[coord_str][:]
            output_vel_all[iter_count*3+1, :] = read_coords_vel_1[coord_str][:]
            output_vel_all[iter_count*3+2, :] = read_coords_vel_2[coord_str][:]
            # debug
            #if iter_count == 0:
            #    print('vel_0 for', vel_iter, 'is:', read_coords_vel_0[vel_iter])
        # iter_count += 1
    #end for
    print('Rank', mpi_rank, 'has written', iter_count+1, 'coordinates')
    # end no interpolation

    # Write coordinates and boundary nodes (file created previously), in serial with baton passing
    comm.Barrier()

    if output_format == "csv" or output_format == "DRM":
        if mpi_rank == 0:
            create_hdf5_csv(output_fname, n_coord, tsteprange, 
                essi_dt, essi_h, gen_vel, gen_acc, gen_dis, 
                essi_lonlat, essi_az, subdomain_origin, subdomain_xyz_len)
            if my_ncoord[0] > 0:
                write_to_hdf5_range_2d(output_fname, '/', 'xyz', my_user_coordinates, my_offset, (my_offset+my_ncoord[0]))
                write_to_hdf5_range(output_fname, '/', 'node_motion_data', output_vel_all[:,tsteprange], my_offset*3, (my_offset+my_ncoord[0])*3)
            if mpi_size > 1:
                comm.send(my_ncoord, dest=1, tag=11)
        else:
            data = comm.recv(source=mpi_rank-1, tag=11)
            if my_ncoord[0] > 0:
                if verbose:
                    print('Rank', mpi_rank, 'start to write data')
                write_to_hdf5_range_2d(output_fname, '/', 'xyz', my_user_coordinates, my_offset, (my_offset+my_ncoord[0]))
                write_to_hdf5_range(output_fname, '/', 'node_motion_data', output_vel_all[:,tsteprange], my_offset*3, (my_offset+my_ncoord[0])*3)
            if mpi_rank != mpi_size-1:
                comm.send(my_ncoord, dest=mpi_rank+1, tag=11) 
    else:
        if mpi_rank == 0:
            print('Invalid output format', output_format)
    
    comm.Barrier()
    if mpi_rank == 0:
        print('Rank', mpi_rank, 'Finished writing data')    
    return

def dframeToDict(dFrame):
    dFrame = list(dFrame.iterrows())
    return {i[1].to_list()[0] : i[1].to_list()[1] for i in dFrame}

def convert_csv(
    csv_fname, ssi_fname, save_path, 
    ref_coord, start_t, end_t, tstep, 
    plot_only, mpi_rank, mpi_size, verbose
) -> None:
    """
    Docstring
    """
    # TODO: THIS FUNCTION
    print('Warning: CSV settings file not implemented yet for convert_sw4drm.py...')
    return

def generate_drm_planes(xyz_len, h):
    """
    plane format = [x0 x1 y0 y1 z0 z1]
    if *0 = *1, that is the plane,
    e.g. if x0 = x1, then we are in plane x=x0=x1

    Parameters
    ----------
    xyz_len : list of float
        Defines the length of the DRM boundary in each axis.
    h : float
        Grid spacing.

    Returns
    -------
    A list of planes that make up the DRM boundary, with no overlaps.
    """
    numLayers = 4
    
    # Compute boundary min/max
    x_min = 0.0
    x_max = xyz_len[0]
    y_min = 0.0
    y_max = xyz_len[1]
    z_min = 0.0
    z_max = xyz_len[2]
    
    plane_data = []
    
    # Planes 1
    z0 = 0
    z1 = z_max
    y0 = y_min
    y1 = y_max - numLayers*h
    for xx in np.arange(x_min, x_min+numLayers*h, h):
        plane_data.append([xx, xx, y0, y1, z0, z1])
    
    # Planes 2
    x0 = x_min
    x1 = x_max - numLayers*h
    for yy in np.arange(y_max-numLayers*h+h, y_max+h, h):
        plane_data.append([x0, x1, yy, yy, z0, z1])
    
    # Planes 3
    y0 = y_min+numLayers*h
    y1 = y_max
    for xx in np.arange(x_max-numLayers*h+h, x_max+h, h):
        plane_data.append([xx, xx, y0, y1, z0, z1])
    
    # Planes 4
    x0 = x_min+numLayers*h
    x1 = x_max
    for yy in np.arange(y_min, y_min+numLayers*h, h):
        plane_data.append([x0, x1, yy, yy, z0, z1])
    
    # Bottom plane
    x0 = x_min+numLayers*h
    x1 = x_max-numLayers*h
    y0 = y_min+numLayers*h
    y1 = y_max-numLayers*h
    for zz in np.arange(z_max-numLayers*h+h, z_max+h, h):
        plane_data.append([x0, x1, y0, y1, zz, zz])
    
    return np.array(plane_data)

def generate_drm_nodes(drm_planes, h):
    """
    Takes a list of DRM planes and turns it into arrays of node coords.
    """
    user_x = []
    user_y = []
    user_z = []

    for plane in drm_planes:
        for xx in np.arange(plane[0], plane[1]+h, h):
            for yy in np.arange(plane[2], plane[3]+h, h):
                for zz in np.arange(plane[4], plane[5]+h, h):
                    user_x.append(xx)
                    user_y.append(yy)
                    user_z.append(zz)
    
    return np.array(user_x), np.array(user_y), np.array(user_z)

def convert_drm(
    csv_fname, ssi_fname, save_path, 
    ref_coord, start_t, end_t, tstep, 
    plot_only, mpi_rank, mpi_size, verbose
) -> None:
    """
    Reads the DRM boundary configuration from the CSV file, and calls 
    the function `generate_acc_dis_time` to extract and save the data.
    
    Parameters
    ----------
    csv_fname : str
        Filename of the CSV file containing the DRM configuration info.
    ssi_fname : str
        Filename of the SW4 ssioutput file containing the simulation data.
    save_path : str
        Filepath to save the hdf5 containing the extracted motions to.
    ref_coord : (3,) np.ndarray
        The DRM subdomain origin in the SW4 coordinate frame.
    """
    if mpi_rank == 0:
        print('Input  CSV [%s]' %csv_fname)
        print('Input ESSI [%s]' %ssi_fname)
    
    # By default, we assume the input and output are in velocity for now
    gen_vel = True
    gen_dis = False
    gen_acc = False
    
    # Load DRM boundary information (in coordinate frame of original SW4 simulation)
    df = np.loadtxt(csv_fname, delimiter=',', dtype=float, skiprows=1)
    drm_xyz_0 = [df[0], df[1], 0.0] # by default we assume drm boundary origin lies on surface
    drm_xyz_len = [df[2], df[3], df[4]]
    drm_h = df[5]
    
    # Generate the list of node coordinates given boundary information
    # - In SW4, the DRM boundary requires four layers of node motions.
    # - These node coordinates are in the coordinate frame of the
    #   DRM boundary (so the origin is 0,0,0).
    drm_planes = generate_drm_planes(drm_xyz_len, drm_h)
    user_x, user_y, user_z = generate_drm_nodes(drm_planes, drm_h)
    n_coord = user_x.shape[0]
    
    if mpi_rank == 0:
        print(f'Generating motions for {n_coord} nodes...')
    
    output_format = 'DRM'
    output_fname = save_path + '/' + output_format + 'NodeMotion.h5'
    
    generate_acc_dis_time(ssi_fname,
        drm_xyz_0, drm_xyz_len,
        user_x, user_y, user_z, n_coord,
        start_t, end_t, tstep,
        gen_vel, gen_acc, gen_dis,
        verbose, plot_only, output_fname,
        mpi_rank, mpi_size, output_format)
    
    return

if __name__ == "__main__":
    # Define defaults
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    os.environ['PYTHONUNBUFFERED'] = 'TRUE'
    verbose=False
    plotonly=False
    use_drm=False
    use_csv=False
    ssi_fname=''
    drm_fname=''
    csv_fname=''
    save_path='./'
    ref_coord=np.zeros(3)
    start_t=0
    end_t=0
    tstep = 1
    
    # Parse arguments
    parser=argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("-p", "--plotonly", help="only generate plots of the input nodes", action="store_true")
    parser.add_argument("-s", "--savepath", help="full path for saving the result files", default="")
    parser.add_argument("-d", "--drm", help="full path to the CSV file with DRM boundary information", default="")
    parser.add_argument("-e", "--essi", help="full path to the SW4 ESSI output file (required)", default="")
    args = parser.parse_args()
    
    if args.verbose:
        verbose=True
    if args.plotonly:
        plotonly=True
    if args.savepath:
        save_path = args.savepath
    if args.drm:
        drm_fname=args.drm
        use_drm=True
    if args.essi:
        ssi_fname=args.essi
    
    comm = MPI.COMM_WORLD
    mpi_size = comm.Get_size()
    mpi_rank = comm.Get_rank()
    
    if mpi_rank == 0:
        print('Running with ', mpi_size, 'MPI processes')
        os.makedirs(save_path, exist_ok=True)
        
        startTime = time.time()
        startTimeStr = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime(startTime))
        print('Start time:', startTimeStr)
    
    # Check that all required data has been provided, and no extra
    if csv_fname == '' and drm_fname == '':
        print('Error: no node coordinate input file is provided, exit...')
        exit(0)
    if csv_fname != '' and drm_fname != '':
        print('Error: both a csv node coordinate file and DRM configuration file were provided, please only provide one, exit...')
        exit(0)
    if ssi_fname == '':
        print('Error: no SW4 ESSI output file is provided, exit...')
        exit(0) 
    
    # EXTRACT THE DATA AND SAVE
    if use_drm:
        convert_drm(drm_fname, ssi_fname, save_path, ref_coord, start_t, end_t, tstep, plotonly, mpi_rank, mpi_size, verbose)
    # elif use_csv:
    #     convert_csv(csv_fname, ssi_fname, save_path, ref_coord, start_t, end_t, tstep, plotonly, mpi_rank, mpi_size, verbose)
    
    if mpi_rank == 0:
        endTime = time.time()
        endTimeStr = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime(endTime))
        print(f'End time: {endTimeStr} (time spent: {endTime-startTime:.2f} s)')