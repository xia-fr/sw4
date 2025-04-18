//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// # Produced at the Lawrence Livermore National Laboratory. 
// # 
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// # 
// # LLNL-CODE-643337 
// # 
// # All rights reserved. 
// # 
// # This file is part of SW4, Version: 1.0
// # 
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991. 
// # 
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details. 
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 

#include "mpi.h"
#include "EW.h"
#include "Require.h"
#include "nearlyEqual.h"

#ifdef USE_HDF5
  #include "hdf5.h"
#endif

#include <cstring>
#include <iostream>
#include <fstream>

#define SQR(x) ((x)*(x))

int computeEndGP( float_sw4 maxval, float_sw4 dh )
{
  // We round up one, so that the end point
  // specified by the user is always included
  // in the domain.  i.e. if z was specified
  // as 15.0, and dh was computed to be 3.33,
  // the number of grid points would be 15.0/3.33 + 1
  // or 5, giving us a new max z = 16.64.
  int pts = 0;
  float_sw4 x = 0.0;
  while (x < maxval && !dbc::nearlyEqual(x, maxval) ){
    x += dh;
    pts++;
  }
  // 1 based indexing
  pts++;
  return pts;
}

// ======================================================================
// [1] DRM loading and setup functions
// ======================================================================

// Processes the DRM command from an input file
//-----------------------------------------------------------------------
void EW::processDRM(char* buffer)
{
#ifdef USE_HDF5
  m_use_DRM = true;
  char dfile[1000];
  bool dfile_set = false;
  bool lon_p_set = false, lat_p_set = false, datum_set = false;
  bool proj_set = false, ellps_set = false;
  bool use_geoprojection = false;
  stringstream proj0;

  double stime, etime;
  stime = MPI_Wtime();

  char* token = strtok(buffer, " \t");
  REQUIRE2(strcmp("DRM", token) == 0, "ERROR: not a DRM...: " << token);
  token = strtok(NULL, " \t");

  if (proc_zero_evzero() )
    cout << endl << "* Processing the DRM command..." << endl;

  // Assume supergrid command has already been processed.
  if( m_sg_damping_order == 6 )
  {
     m_ghost_points = 3;
     m_ppadding = 3;
  }

  while (token != NULL)
  {
    // while there are tokens in the string still
    if (startswith("#", token) || startswith(" ", buffer))
      // Ignore comments
      break;
    if (startswith("file=",token)){
      token += 5;
      strncpy(dfile, token, 1000);
      dfile_set = true;
    } 
    else if( startswith("proj=",token))
    {
      token += 5;
      // accumulate new style string
      proj0 << "+proj=" << token;
      use_geoprojection = true;
      proj_set = true;
    }
    else if( startswith("datum=",token))
    {
      token += 6;
      proj0 << " +datum=" << token;
      use_geoprojection = true;
      datum_set = true;
    }
    else if( startswith("ellps=",token))
    {
      token +=6;
      // accumulate new style string
      proj0 << " +ellps=" << token;
      use_geoprojection = true;
      ellps_set = true;
    }
    else if( startswith("lon_p=",token))
    {
      token += 6;
      proj0 << " +lon_0=" << atof(token);
      use_geoprojection = true;
      lon_p_set=true;
    }
    else if( startswith("lat_p=",token))
    {
      token += 6;
      proj0 << " +lat_0=" << atof(token);
      use_geoprojection = true;
      lat_p_set=true;
    }
    else if( startswith("scale=",token))
    {
      token += 6;
      proj0 << " +scale=" << atof(token);
      use_geoprojection = true;
    }
    else
    {
      badOption("DRM", token);
    }
    token = strtok(NULL, " \t");
  }

  // From the hdf5 file containing the DRM info, we read in:
  // - grid size
  // - original SW4 simulation azimuth
  // - original SW4 simulation origin lon/lat
  // - size of the DRM domain
  // This information in combination with the supergrid is used
  // to construct a grid model around the DRM domain.
  // Grid's (x, y, z, h, lat, lon, az) all read/computed from DRM file.
  float_sw4 x = 0.0, y = 0.0, z = 0.0, h = 0.0;
  float_sw4 x0 = 0.0, y0 = 0.0;
  double lat = 37.0, lon = -118.0; // Default is NTS
  mGeoAz = 0; // default azimuth
  m_DRM_dt = 1.0; // default dt
  if ( dfile_set ){
    // Read in the grid information from the hdf5 file
    readDRMHDF5_grid( dfile, x, y, z, h, x0, y0, lat, lon );

    // Compute the correct latitude and longitude of the DRM simulation origin based on
    // the lat/lon of the original SW4 simulation origin,the location of the DRM boundary origin,
    // and the size of the supergrid
    if (use_geoprojection) { // default arguments for proj4 projection
      if (!proj_set){ // Default projection: Universal Transverse Mercator (UTM)
        proj0 << " +proj=utm";
      }
      if (!ellps_set && !datum_set){ // default ellipse
        proj0 << " +ellps=WGS84";
      }
      if (!lon_p_set){ // if lon_p not given, use lon
        proj0 << " +lon_0=" << lon;
      }
      if (!lat_p_set){ // if lat_p not given, use lat
        proj0 << " +lat_0=" << lat;
      }
      // hard code units to be in meters
      proj0 << " +units=m";
    }
    #if !defined(ENABLE_PROJ)
      CHECK_INPUT( !use_geoprojection, "ERROR: need to configure SW4 with proj=yes to use projections "
                  "from the PROJ library (version 6 or later)");
    #endif
    if( use_geoprojection ){
      // 1. Create temporary GeographicProjection using latp and lonp of original full grid
      GeographicProjection* temp_geoproj = new GeographicProjection( lon, lat, proj0.str(), mGeoAz );

      // 2. Compute the DRM simulation grid origin lat and lon using the temporary geoproj object
      //    and set the computed values to mLonOrigin and mLatOrigin.
      //    - The origin of the DRM simulation grid origin is: (x0-6h-(Nsg-1)h, y0-6h-(Nsg-1)h, 0)
      float_sw4 x0_drm, y0_drm;
      x0_drm = x0 - 6*h - m_supergrid_width;
      y0_drm = y0 - 6*h - m_supergrid_width;
      if ( x0_drm < 0.0 || y0_drm < 0.0){
        cerr << "Warning: x0_drm=" << x0_drm << ", y0_drm=" << y0_drm << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      temp_geoproj->computeGeographicCoord(x0_drm, y0_drm, mLonOrigin, mLatOrigin);

      // 3. Create the global geoprojection object for this simulation
      m_geoproj = new GeographicProjection( mLonOrigin, mLatOrigin, proj0.str(), mGeoAz );
    } else {
      m_geoproj = static_cast<GeographicProjection*>(0);
      mLatOrigin = lat;
      mLonOrigin = lon;
    }
    if (proc_zero())
        printf("> Computed DRM grid origin lat = %.8f, lon = %.8f\n", mLatOrigin, mLonOrigin);

    if( mVerbose >= 1 && proc_zero() )
    {
      cout << "\n=== Setup DRM domain boundaries: ===" << endl;
      cout << "x: " << m_DRM_xMin << " to " << m_DRM_xMax << endl;
      cout << "i: " << m_DRM_iMin << " to " << m_DRM_iMax << endl;
      cout << "====================================" << endl;
      cout << "y: " << m_DRM_yMin << " to " << m_DRM_yMax << endl;
      cout << "j: " << m_DRM_jMin << " to " << m_DRM_jMax << endl;
      cout << "====================================" << endl;
      cout << "z: Surface to " << m_DRM_zMax << endl;
      cout << "k: 0 to " << m_DRM_kMax << endl;
      cout << "====================================\n" << endl;
    }

    // Now we do the other parts of constructing the simulation grid
    constructDRMgrid( x, y, z, h, lat, lon );

    // // From the hdf5 file containing the DRM info, we read in:
    // // - 3 component time series for each grid point on the DRM boundary
    // // - the corresponding xyz coordinate of each g.p. on the DRM boundary,
    // //   relative to the DRM boundary origin
    // // Notes:
    // // - relative to this DRM simulation's grid origin, the DRM boundary origin 
    // //   will have the xyz coordinates: ((Nsg-1)*h + 5*h, (Nsg-1)*h + 5*h, 0.0)
    // readDRMHDF5_data( drm_filename );
  } else {
    cerr << "Error: DRM data filepath not set." << endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  etime = MPI_Wtime();
  if (proc_zero())
    cout << "* Processed DRM command, took " << etime-stime << "seconds.\n" << endl;
  fflush(stdout);
  // end processing DRM input command
#else // if SW4 is not compiled with HDF5
  if (proc_zero())
    cout << "Loading DRM data which requires HDF5, but sw4 is not compiled with HDF5!"<< endl;
#endif
}

// 
//-----------------------------------------------------------------------
void EW::readDRMHDF5_grid(
  char *fname,
  float_sw4 &x, float_sw4 &y, float_sw4 &z, float_sw4 &h,
  float_sw4 &x0, float_sw4 &y0,
  double &lat, double &lon )
{
#ifdef USE_HDF5
  int ret;
  double stime, etime;
  stime = MPI_Wtime();

  // Setup for parallel n readers
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int nreader = m_nwriters;
  if (nreader <= 0) 
      nreader = 1;
  if (nreader > world_size)
      nreader = world_size;
  
  if ( world_size % nreader != 0){
    printf("Warning: %i processes not evenly divided by %i readers, may run into issues", 
      world_size, nreader);
  }
  if ( world_rank == 0 )
      printf("> Using %i readers to access data\n", nreader);

  int read_color = world_rank % (world_size / nreader) == 0 ? 0 : 1;
  int node_color = world_rank / (world_size / nreader);
  int read_rank, read_size;
  MPI_Comm read_comm, node_comm;
  MPI_Comm_split(MPI_COMM_WORLD, read_color, world_rank, &read_comm);
  MPI_Comm_split(MPI_COMM_WORLD, node_color, world_rank, &node_comm);
  MPI_Comm_rank(MPI_COMM_WORLD, &read_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &read_size);
  // end setup for parallel n readers

  // hid_t := HDF5 library signed integer handle for HDF5 objects
  hid_t fid, fapl, dset;
  string gname;
  double temp_shape[3];

  // Only color 0 reads data, then broadcast to all other processes
  if (read_color == 0) {
    // Open file 'fname'
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, read_comm, MPI_INFO_NULL);
    fid = H5Fopen(fname, H5F_ACC_RDONLY, fapl);

    // Printing diagnostics
    if (fid <= 0) {
      cerr << "DRM HDF5 file " << fname << " not found" << endl;
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    if ( world_rank == 0 )
      printf("> Opened DRM file '%s'\n", fname);

    // Open groups and read data
    string groupArray[] = {"h", "dt", "grid_az"};
    readDRMHDF5_double( fid, groupArray[0].c_str(), h );
    if ( world_rank == 0 ) 
      printf("> Read %s = %.2f\n", groupArray[0].c_str(), h);

    readDRMHDF5_double( fid, groupArray[1].c_str(), m_DRM_dt );
    if ( world_rank == 0 ) 
      printf("> Read %s = %.8f\n", groupArray[1].c_str(), m_DRM_dt);

    readDRMHDF5_double( fid, groupArray[2].c_str(), mGeoAz );
    if ( world_rank == 0 ) 
      printf("> Read %s = %.8f\n", groupArray[2].c_str(), mGeoAz);

    // Read original SW4 simulation grid origin lon/lat
    double templl[2];
    gname = "grid_lon_lat";
    dset = H5Dopen(fid, gname.c_str(), H5P_DEFAULT);
    if (dset < 0) {
      printf("Error with DRM HDF5 file, no '%s' dataset found!\n", gname);
      MPI_Abort( MPI_COMM_WORLD, 1 );
    } else {
      ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &templl);
      ASSERT(ret >= 0);
      H5Dclose(dset);
    }
    lon = templl[0];
    lat = templl[1];
    if ( world_rank == 0 )
        printf("> Read original grid lat = %.8f, lon = %.8f\n", lat, lon);

    // Read DRM boundary origin relative to original SW4 simulation grid
    double temp_orig[3];
    gname = "drm_origin";
    dset = H5Dopen(fid, gname.c_str(), H5P_DEFAULT);
    if (dset < 0) {
      printf("Error with DRM HDF5 file, no '%s' dataset found!\n", gname);
      MPI_Abort( MPI_COMM_WORLD, 1 );
    } else {
      ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &temp_orig);
      ASSERT(ret >= 0);
      H5Dclose(dset);
    }
    x0 = temp_orig[0];
    y0 = temp_orig[1];
    if ( world_rank == 0 )
        printf("> Read DRM origin (x,y,z) = %.2f %.2f %.2f\n", x0, y0, temp_orig[2]);

    // Read DRM boundary shape (ostensibly, in meters with grid spacing h)
    gname = "drm_shape";
    dset = H5Dopen(fid, gname.c_str(), H5P_DEFAULT);
    if (dset < 0) {
      printf("Error with DRM HDF5 file, no '%s' dataset found!\n", gname);
      MPI_Abort( MPI_COMM_WORLD, 1 );
    } else {
      ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &temp_shape);
      ASSERT(ret >= 0);
      H5Dclose(dset);
    }
    if ( world_rank == 0 )
        printf("> Read DRM domain shape (x,y,z) = %.2f %.2f %.2f\n", temp_shape[0], temp_shape[1], temp_shape[2]);

    // Compute x,y,z for DRM simulation grid based on DRM boundary shape and supergrid size
    if ( !m_use_sg_width ){
      m_supergrid_width = (m_sg_gp_thickness-1)*h;
      m_use_sg_width = true;
    }
    x = temp_shape[0] + 2*m_supergrid_width + 12*h;
    y = temp_shape[1] + 2*m_supergrid_width + 12*h;
    z = temp_shape[2] + m_supergrid_width + 6*h;
    if ( world_rank == 0 ) 
        printf("> Computed total simulation grid shape (x,y,z) = %.2f %.2f %.2f\n", x, y, z);
    
    H5Fclose(fid);
  }
  // Broadcast to other processors
  MPI_Bcast(&m_DRM_dt, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&mGeoAz, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&m_supergrid_width, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&m_sg_gp_thickness, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&m_use_sg_width, 1, MPI_CXX_BOOL, 0, node_comm);

  MPI_Bcast(&h, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&lat, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&lon, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&x0, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&y0, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&x, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&y, 1, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&z, 1, MPI_DOUBLE, 0, node_comm);

  MPI_Bcast(&temp_shape, 3, MPI_DOUBLE, 0, node_comm);

  MPI_Comm_free(&node_comm);
  MPI_Comm_free(&read_comm);

  // Save DRM min/max values
  m_DRM_xMin = m_supergrid_width + 6*h;
  m_DRM_yMin = m_supergrid_width + 6*h;
  m_DRM_xMax = m_DRM_xMin + temp_shape[0];
  m_DRM_yMax = m_DRM_yMin + temp_shape[1];
  m_DRM_zMax = temp_shape[2];

  // Compute corresponding index values
  m_DRM_iMin = static_cast<int>( floor( m_DRM_xMin/h+1) );
  m_DRM_iMax = static_cast<int>( floor( m_DRM_xMax/h+1) );
  m_DRM_jMin = static_cast<int>( floor( m_DRM_yMin/h+1) );
  m_DRM_jMax = static_cast<int>( floor( m_DRM_yMax/h+1) );
  m_DRM_kMax = static_cast<int>( floor( m_DRM_zMax/h+1) );
  
#endif // USE_HDF5
}

// By default, float_sw4 is of the type 'double', which is what we want
// to be using for the DRM data...
//-----------------------------------------------------------------------
void EW::readDRMHDF5_double( hid_t fid, const char *dname, float_sw4 &data )
{
#ifdef USE_HDF5
  hid_t dset;
  int ret;
  // Open group and read data
  dset = H5Dopen(fid, dname, H5P_DEFAULT);
  if (dset < 0) {
    printf("Error with DRM HDF5 file, no '%s' dataset found!\n", dname);
    MPI_Abort( MPI_COMM_WORLD, 1 );
  } else {
    ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    ASSERT(ret >= 0);
    H5Dclose(dset);
  }
#endif // USE_HDF5
}

// Carries out the grid setup in EW::processGrid(char* buffer) using grid
// parameters loaded previously from DRM hdf5 file.
// By construction, the information we have for specifying the grid are:
// - grid spacing h
// - the total width of the grid in 3 axes, x, y, z
//-----------------------------------------------------------------------
void EW::constructDRMgrid( float_sw4 &x, float_sw4 &y, float_sw4 &z,
  float_sw4 &h, float_sw4 &lat, float_sw4 &lon )
{
  int nx = 0, ny = 0, nz = 0;
  int nxprime, nyprime, nzprime;
  float_sw4 xprime, yprime, zprime;
  nxprime = computeEndGP(x, h);
  nyprime = computeEndGP(y, h);
  nzprime = computeEndGP(z, h);

  if (proc_zero_evzero() && mVerbose >=3)
      printf("**** Setting up the grid for a non-periodic problem\n");
  if (nxprime != nx && proc_zero_evzero())
    cout << "> Setting nx to " << nxprime << " to be consistent with h=" << h << endl;
  if (nyprime != ny && proc_zero_evzero())
    cout << "> Setting ny to " << nyprime << " to be consistent with h=" << h << endl;
  if (nzprime != nz && proc_zero_evzero())
    cout << "> Setting nz to " << nzprime << " to be consistent with h=" << h << endl;

  // -------------------------------------------------------------
  // Now we adjust the geometry bounds based on the actual 
  // number of grid points used in each dimension.
  // -------------------------------------------------------------
  xprime = (nxprime-1)*h;
  zprime = (nzprime-1)*h;
  yprime = (nyprime-1)*h;

  float_sw4 eps = 1.e-9*sqrt(SQR(xprime)+SQR(yprime)+SQR(zprime));
  if( sizeof(float_sw4)==4)
      eps=eps*1e4;

  if (fabs(xprime-x) > eps && proc_zero_evzero())
    cout << "* Changing x from " << x << " to " << xprime << " to be consistent with h=" << h << endl;
  if (fabs(zprime-z) > eps && proc_zero_evzero())
    cout << "* Changing z from " << z << " to " << zprime << " to be consistent with h=" << h << endl;
  if (fabs(yprime-y) > eps && proc_zero_evzero())
    cout << "* Changing y from " << y << " to " << yprime << " to be consistent with h=" << h << endl;

  // Finally we save these values
  m_nx_base = nxprime;
  m_ny_base = nyprime;
  m_nz_base = nzprime;
  m_h_base = h;
  m_global_xmax = xprime;
  m_global_ymax = yprime;
  m_global_zmax = zprime;
}

// Reads in the displacement DRM data and creates the corresponding
// DRMSource objects.
// Note: all XYZ coordinates in this will be offset, because the coordinates
//       in the hdf5 file are relative to the DRM boundary origin
//-----------------------------------------------------------------------
void EW::readDRMHDF5_data( char* buffer )
{
  // We grab the name of the hdf5 file again
  char dfile[1000];
  bool dfile_set = false;
  double stime, etime1, etime2;

  char* token = strtok(buffer, " \t");
  REQUIRE2(strcmp("DRM", token) == 0, "ERROR: not a DRM...: " << token);
  token = strtok(NULL, " \t");

  while (token != NULL)
  {
    // while there are tokens in the string still
    if (startswith("#", token) || startswith(" ", buffer))
      // Ignore comments
      break;
    if (startswith("file=",token)){
      token += 5;
      strncpy(dfile, token, 1000);
      dfile_set = true;
    }
    token = strtok(NULL, " \t");
  }

  if ( !dfile_set ){
    cerr << "Error: DRM data filepath not set." << endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  if ( proc_zero_evzero() )
    cout << endl << "* Loading DRM node motions..." << endl;
  stime = MPI_Wtime();

#ifdef USE_HDF5
  string gname;
  int status;
  bool myPoint;

  // hid_t := HDF5 library signed integer handle for HDF5 objects
  hid_t fid, fapl, dset;

  int ndims = 2;
  float_sw4 *xyz_data, *node_data;
  hid_t xyz_set, xyz_space, xyz_memspace;
  hsize_t xyz_shape[ndims];
  hsize_t xyz_offset[ndims], xyz_count[ndims];

  hid_t node_set, node_space, node_memspace;
  hsize_t node_shape[ndims];
  hsize_t node_offset[ndims], node_count[ndims];

  // Setup for parallel n readers
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int nreader = m_nwriters;
  if (nreader <= 0) 
      nreader = 1;
  if (nreader > world_size)
      nreader = world_size;
  
  if ( world_size % nreader != 0){
    printf("Warning: %i processes not evenly divided by %i readers", 
      world_size, nreader);
  }
  if ( world_rank == 0 )
      printf("> Using %i readers to access data\n", nreader);

  int read_color = world_rank % (world_size / nreader) == 0 ? 0 : 1;
  int node_color = world_rank / (world_size / nreader);
  int read_rank, read_size;
  MPI_Comm read_comm, node_comm;
  MPI_Comm_split(MPI_COMM_WORLD, read_color, world_rank, &read_comm);
  MPI_Comm_split(MPI_COMM_WORLD, node_color, world_rank, &node_comm);
  MPI_Comm_rank(MPI_COMM_WORLD, &read_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &read_size);
  // end setup for parallel n readers
  
  // Only color 0 reads data, then broadcast to all other processes
  if (read_color == 0) {
    // Open file
    stime = MPI_Wtime();
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, read_comm, MPI_INFO_NULL);
    fid = H5Fopen(dfile, H5F_ACC_RDONLY, fapl);

    // Printing diagnostics
    if (fid <= 0) {
      cerr << "DRM HDF5 file " << dfile << " not found" << endl;
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    if ( proc_zero() )
      printf("> Opened DRM file '%s'\n", dfile);

    // Grab the dataset containing the node motion data
    gname = "node_motion_data";
    node_set = H5Dopen(fid, gname.c_str(), H5P_DEFAULT);

    // Get the dataspace and shape of the data
    node_space = H5Dget_space(node_set); 
    status = H5Sget_simple_extent_dims(node_space, node_shape, NULL);
    if ( proc_zero() ) 
      printf("> Group 'node_motion_data': dimensions %d x %d \n", node_shape[0], node_shape[1]);
    node_count[0] = node_shape[0]; // for the 3 components of motion
    node_count[1] = node_shape[1];
    node_offset[0] = 0;
    node_offset[1] = 0;

    // Define the memory dataspace (same size as selection)
    node_memspace = H5Screate_simple(ndims, node_count, NULL);

    // Grab the dataset containing the xyz coordinate data
    gname = "xyz";
    xyz_set = H5Dopen(fid, gname.c_str(), H5P_DEFAULT);
    
    // Get the dataspace and shape of the data
    xyz_space = H5Dget_space(xyz_set); 
    status = H5Sget_simple_extent_dims(xyz_space, xyz_shape, NULL);
    if ( proc_zero() ) 
      printf("> Group 'xyz': dimensions %d x %d \n", xyz_shape[0], xyz_shape[1]);
    fflush(stdout);
    xyz_count[0] = xyz_shape[0];
    xyz_count[1] = xyz_shape[1];
    xyz_offset[0] = 0;
    xyz_offset[1] = 0;

    // Define the memory dataspace (same size as selection)
    xyz_memspace = H5Screate_simple(ndims, xyz_count, NULL);

    // Allocate temporary arrays to hold the data once it has been read
    // - index using row major syntax: data[i][j] = data[i*data_shape[1] + j]
    xyz_data = (float_sw4 *)malloc(xyz_shape[0]*xyz_shape[1]*sizeof(float_sw4));
    node_data = (float_sw4 *)malloc(node_shape[0]*node_shape[1]*sizeof(float_sw4));

    // Select the hyperslab and read data from the selection
    status = H5Sselect_hyperslab(xyz_space, H5S_SELECT_SET, xyz_offset, NULL, xyz_count, NULL);
    ASSERT(status >= 0);
    status = H5Dread(xyz_set, H5T_NATIVE_DOUBLE, xyz_memspace, xyz_space, H5P_DEFAULT, xyz_data);
    ASSERT(status >= 0);

    // Read the node motion data
    status = H5Sselect_hyperslab(node_space, H5S_SELECT_SET, node_offset, NULL, node_count, NULL);
    ASSERT(status >= 0);
    status = H5Dread(node_set, H5T_NATIVE_DOUBLE, node_memspace, node_space, H5P_DEFAULT, node_data);
    ASSERT(status >= 0);

    // Close resources
    H5Sclose(xyz_memspace);
    H5Sclose(xyz_space);
    H5Dclose(xyz_set);
    H5Sclose(node_memspace);
    H5Sclose(node_space);
    H5Dclose(node_set);
    H5Fclose(fid);
  }
  // Broadcast data to the other processors
  MPI_Bcast(&node_shape, ndims, MPI_DOUBLE, 0, node_comm);
  MPI_Bcast(&xyz_shape, ndims, MPI_DOUBLE, 0, node_comm);

  // Allocate data for the processors that weren't reading
  if (read_color != 0) {
    xyz_data = (float_sw4 *)malloc(xyz_shape[0]*xyz_shape[1]*sizeof(float_sw4));
    node_data = (float_sw4 *)malloc(node_shape[0]*node_shape[1]*sizeof(float_sw4));
  }

  // Broadcast node motion data
  MPI_Bcast(xyz_data, xyz_shape[0]*xyz_shape[1], m_mpifloat, 0, node_comm);
  MPI_Bcast(node_data, node_shape[0]*node_shape[1], m_mpifloat, 0, node_comm);

  MPI_Comm_free(&node_comm);
  MPI_Comm_free(&read_comm);

  etime1 = MPI_Wtime();
  if ( proc_zero() ){
    printf("* Loaded %i DRM node motions, took %.8f seconds.\n", xyz_shape[0], etime1-stime);
    printf("> Now generating DRM source objects...\n");
    fflush(stdout);
  }

  // We set the number of time steps here while we have the information
  int e = 0;
  if (mTimeIsSet[e])
    mNumberOfTimeSteps[e] = node_shape[1];

  // if ( m_myRank == 1 ){
  //   int i0=1500;
  //   printf("[x0]: %.4f\n", xyz_data[(i0+0)*xyz_shape[1] + 0]);
  //   printf("[y0]: %.4f\n", xyz_data[(i0+0)*xyz_shape[1] + 1]);
  //   printf("[z0]: %.4f\n", xyz_data[(i0+0)*xyz_shape[1] + 2]);
  //   for(int ii=0;ii<39;ii++){
  //     printf("[%i]: %.8e\n", ii+1, node_data[(i0*3+1)*node_shape[1] + ii]);
  //   }
  // }
  
  float_sw4 x = 0.0, y = 0.0, z = 0.0;
  DRMSource* drmsrcPtr;
  float_sw4* par;
  for ( int ix=0; ix<xyz_shape[0]; ix++ ){
    // The xyz coordinates of the drm nodes are saved with respect to the DRM boundary origin,
    // so we have to offset them to be in the proper coordinate frame.
    x = xyz_data[ix*xyz_shape[1] + 0] + m_supergrid_width + 6*mGridSize[0];
    y = xyz_data[ix*xyz_shape[1] + 1] + m_supergrid_width + 6*mGridSize[0];
    z = xyz_data[ix*xyz_shape[1] + 2];
  
    // Create new DRM source and add it to list of source terms
    for ( int comp=0; comp<3; comp++){
      par = &node_data[(ix*3+comp)*node_shape[1]];
      drmsrcPtr = new DRMSource(this, x, y, z, comp+1,
        par, node_shape[1], m_DRM_dt);
      if (!drmsrcPtr->myPoint()){
        delete drmsrcPtr;
      } else { // We only save this source object if it belongs to the proc
        m_DRMSources.push_back(drmsrcPtr);
      }
    }
  }

  delete[] xyz_data;
  delete[] node_data;

  MPI_Barrier(MPI_COMM_WORLD);
  etime2 = MPI_Wtime();
  if ( proc_zero() ){
    printf("* Generated DRM source objects, took %.8f seconds.\n", xyz_shape[0], etime2-etime1);
    printf("* Total time to process DRM node motions: %.8f seconds.\n\n", xyz_shape[0], etime2-stime);
    fflush(stdout);
  }
#endif // USE_HDF5
}

// ======================================================================
// [2] DRM computation functions
// ======================================================================

// Loads all of the DRM sources from the vector into displacement field U
// TODO: optimize?
// also TODO: remove the "inner/outer" differentiation for this
//-----------------------------------------------------------------------
void EW::DRMSourcesToUField( vector<Sarray>& a_U, float_sw4 t, bool inner )
{
  // Presently, we assume all the DRM happens in the top-most grid
  int g = 0;

  // We iterate over all the source objects in m_DRMSources
  float_sw4 value;
  for (int i = 0; i < m_DRMSources.size(); i++){
    // If this source belongs to the processor
    if (m_DRMSources[i] -> myPoint()){
      // Save into corresponding location in displacement array
      int comp = m_DRMSources[i]->getComp();
      int i0 = m_DRMSources[i]->m_i0;
      int j0 = m_DRMSources[i]->m_j0;
      int k0 = m_DRMSources[i]->m_k0;

      // If we're loading displacements that are on or inside of boundary
      if (inner && point_inside_DRM_boundary(i0, j0, k0, g)){
        // Compute the source value for time t
        a_U[g](comp, i0, j0, k0) = m_DRMSources[i] -> getU(t, mDt);
      }
      // If we're loading displacements on the outside of boundary
      else if (!inner && point_outside_DRM_boundary(i0, j0, k0, g))
        a_U[g](comp, i0, j0, k0) = 0.0; 
      // NOTE: SARRAY COMPONENTS INDEX FROM ONE !!!!!!!!!!!!!
      // 1: X, 2: Y, 3: Z !!!!!

      // // DEBUG:
      // cout << "Processor " << m_myRank << ": ";
      // cout << "Loaded value " << srcVal << " into ";
      // cout << "x" << m_DRMSources[i]->getX0() << " y" << m_DRMSources[i]->getY0() 
      //      << " z" << m_DRMSources[i]->getZ0() << " at time " << t;
      // cout << " for component=" << comp << " i0=" << i0 << " j0=" << j0 << " k0=" << k0 << endl;
    }
  }
  communicate_array( a_U[g], g );
}

// Computes the force for DRM given displacement fields
// TODO: optimize?
//---------------------------------------------------------------------------
void EW::evalForce(vector<Sarray> & a_Up, vector<Sarray> & a_U, vector<Sarray> & a_Um,
  vector<Sarray> & a_Rho, vector<Sarray> & a_Lu, vector<Sarray> & a_F, bool inner )
{
  int g = 0;
  int ifirst, ilast, jfirst, jlast, kfirst, klast;
  ifirst = m_iStart[g];
  ilast  = m_iEnd[g];
  jfirst = m_jStart[g];
  jlast  = m_jEnd[g];
  kfirst = m_kStart[g];
  klast  = m_kEnd[g];
  float_sw4 dt2 = mDt*mDt;

  a_F[g].set_to_zero();

  for (int ii = ifirst; ii <= ilast; ii++)
  for (int jj = jfirst; jj <= jlast; jj++)
  for (int kk = kfirst; kk <= klast; kk++)
  for (int cc = 1; cc <= 3; cc++){
    if ( interior_point_in_proc(ii, jj, g) ){
      float_sw4 rhodt2 = a_Rho[0](ii, jj, kk)/dt2;
      a_F[g](cc, ii,jj,kk) = rhodt2*(a_Up[g](cc,ii,jj,kk) - 2*a_U[g](cc,ii,jj,kk) + a_Um[g](cc,ii,jj,kk)) 
                           - a_Lu[g](cc,ii,jj,kk);
    }
  }
  communicate_array( a_F[g], g );
}

// Modifies force array to remove forcing terms outside the DRM boundary
// TODO: optimize this by instead removing the relevant GridPointSource objects
//       instead of zero-ing out the force array
//---------------------------------------------------------------------------
void EW::maskForce( vector<Sarray> & a_F )
{
  int g = 0;
  for (int ii = m_iStart[0]; ii <= m_iEnd[0]; ii++)
  for (int jj = m_jStart[0]; jj <= m_jEnd[0]; jj++)
  for (int kk = m_kStart[0]; kk <= m_kEnd[0]; kk++)
    if ( interior_point_in_proc(ii, jj, g) ){
      if ( point_on_DRM_boundary(ii,jj,kk,0,true) || point_outside_DRM_boundary(ii,jj,kk,0) )
        for (int cc = 1; cc <= 3; cc++){
          a_F[g](cc,ii,jj,kk) = 0.0;
        }
    }
}

// Checks if the point (i,j,k) lies on a DRM boundary that we want to
// compute forces for (aka the defined boundary in DRMinfo)
//-----------------------------------------------------------------------
bool EW::point_on_DRM_boundary(int a_i, int a_j, int a_k, int a_g, bool inner)
{
  bool retval = true;
  if (inner){ // Checks if the point is on the inner boundary
    // It's easier to check if the point ISNT on the DRM boundary
    if (a_i < m_DRM_iMin || a_i > m_DRM_iMax ||
        a_j < m_DRM_jMin || a_j > m_DRM_jMax ||
        a_k > m_DRM_kMax)
          retval = false;
    if (a_i > m_DRM_iMin+1 && a_i < m_DRM_iMax-1 &&
        a_j > m_DRM_jMin+1 && a_j < m_DRM_jMax-1 &&
        a_k < m_DRM_kMax-1)
          retval = false;
  } else { // Checks if the point is on the outer boundary
    if (a_i < m_DRM_iMin-2 || a_i > m_DRM_iMax+2 ||
        a_j < m_DRM_jMin-2 || a_j > m_DRM_jMax+2 ||
        a_k > m_DRM_kMax+2)
          retval = false;
    if (a_i > m_DRM_iMin-1 && a_i < m_DRM_iMax+1 &&
        a_j > m_DRM_jMin-1 && a_j < m_DRM_jMax+1 &&
        a_k < m_DRM_kMax+1)
          retval = false;
  }
  return retval;
}

// Checks if the point (i,j,k) lies inside the DRM boundary (inclusive)
//-----------------------------------------------------------------------
bool EW::point_inside_DRM_boundary(int a_i, int a_j, int a_k, int a_g)
{
  bool retval = false;
  if (a_i >= m_DRM_iMin && a_i <= m_DRM_iMax &&
      a_j >= m_DRM_jMin && a_j <= m_DRM_jMax &&
      a_k <= m_DRM_kMax)
    retval = true;
  return retval;
}

// Checks if the point (i,j,k) lies outside the DRM boundary (exclusive)
//-----------------------------------------------------------------------
bool EW::point_outside_DRM_boundary(int a_i, int a_j, int a_k, int a_g)
{
  bool retval = false;
  if (a_i < m_DRM_iMin || a_i > m_DRM_iMax ||
      a_j < m_DRM_jMin || a_j > m_DRM_jMax ||
      a_k > m_DRM_kMax)
    retval = true;
  return retval;
}

//------------------------------------------------------------------------
void EW::debugPrinter_drm( Sarray & a_Var )
{
  int i_dbg = 47;
  int j_dbg = 52;
  int k_dbg = 21;
  cout << a_Var(1,i_dbg,j_dbg,k_dbg-4) << " ";
  cout << a_Var(1,i_dbg,j_dbg,k_dbg-3) << " ";
  cout << a_Var(1,i_dbg,j_dbg,k_dbg-2) << " ";
  cout << a_Var(1,i_dbg,j_dbg,k_dbg-1) << " ";
  cout << a_Var(1,i_dbg,j_dbg,k_dbg-0) << " ";
  cout << a_Var(1,i_dbg,j_dbg,k_dbg+1) << endl;
}