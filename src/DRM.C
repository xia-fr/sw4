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
// [DRM-9] added for std::vector / std::min used in readDRMHDF5_data
#include <vector>
#include <algorithm>

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
    // [SW-4] NEW DRM options: sliding time window for the node motions.
    //   window=N  : number of time intervals kept in memory per source (0 = whole series)
    //   overlap=M : extra samples fitted on each side of a window
    else if( startswith("window=",token))
    {
      token += 7;
      m_DRM_window = atoi(token);
      CHECK_INPUT( m_DRM_window == 0 || m_DRM_window >= 8,
                   "ERROR: DRM window must be 0 (whole series) or >= 8, not " << m_DRM_window );
    }
    else if( startswith("overlap=",token))
    {
      token += 8;
      m_DRM_overlap = atoi(token);
      CHECK_INPUT( m_DRM_overlap >= 0, "ERROR: DRM overlap must be >= 0, not " << m_DRM_overlap );
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
    if (proc_zero())  // [SW-4]
        printf("> DRM time window = %d intervals (0 = whole series), overlap = %d samples\n",
               m_DRM_window, m_DRM_overlap);

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
      printf("Error with DRM HDF5 file, no '%s' dataset found!\n", gname.c_str()); // [SW-9] BUG FIX: was passing std::string to %s
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
      printf("Error with DRM HDF5 file, no '%s' dataset found!\n", gname.c_str()); // [SW-9] BUG FIX: was passing std::string to %s
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
      printf("Error with DRM HDF5 file, no '%s' dataset found!\n", gname.c_str()); // [SW-9] BUG FIX: was passing std::string to %s
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
  // ======================================================================
  // [DRM-1] REWRITTEN: per-rank, batched reading of node_motion_data.
  //   OLD: reader ranks malloc'd the ENTIRE node_motion_data dataset
  //        (3929724 x 19403 doubles = ~610 GB), read it all, MPI_Bcast it to
  //        the other ranks (non-reader ranks also malloc'd the full size),
  //        and EVERY rank then constructed a DRMSource for EVERY node and
  //        threw away the ones it did not own.
  //   NEW: every rank
  //        (a) reads only the small xyz table (~31 MB),
  //        (b) works out which DRM nodes lie in its own patch   [DRM-4],
  //        (c) builds DRMSource objects only for those nodes [SW-5],
  //        (d) loads their time series one window at a time [SW-6]
  //            (the earlier full-series batched read [DRM-6] is replaced).
  //   There is no broadcast and no reader/node communicator any more.
  //   REQUIRES [DRMS-1] in DRMSource.C (removal of the per-source
  //   MPI_Allreduce); otherwise ranks would no longer make matching
  //   collective calls and the run would hang.
  // ======================================================================
  const int ndims = 2;
  int status;
  hid_t fid, fapl;
  hid_t xyz_set, xyz_space;
  hid_t node_set, node_space;
  hsize_t xyz_shape[2], node_shape[2];
  float_sw4 *xyz_data;

  // [DRM-2] float_sw4 is normally double; pick the matching HDF5 memory type
  // instead of hard-coding NATIVE_DOUBLE (HDF5 converts from the on-disk type).
  const hid_t h5type = (sizeof(float_sw4) == 8) ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;

  // [DRM-2] Every rank opens the file (collective open over MPI_COMM_WORLD)
  // and then does its own independent reads (default transfer property).
  // The access property list is now closed (it was leaked before).
  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  fid = H5Fopen(dfile, H5F_ACC_RDONLY, fapl);
  H5Pclose(fapl);
  if (fid < 0) {
    cerr << "DRM HDF5 file " << dfile << " not found" << endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  if ( proc_zero() )
    printf("> Opened DRM file '%s'\n", dfile);

  // [DRM-3] Dataset-open and shape checks (the old code did not check
  // H5Dopen for either dataset).
  node_set = H5Dopen(fid, "node_motion_data", H5P_DEFAULT);
  xyz_set  = H5Dopen(fid, "xyz", H5P_DEFAULT);
  if (node_set < 0 || xyz_set < 0) {
    cerr << "Error: DRM HDF5 file is missing 'node_motion_data' and/or 'xyz'" << endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  node_space = H5Dget_space(node_set);
  xyz_space  = H5Dget_space(xyz_set);
  status = H5Sget_simple_extent_dims(node_space, node_shape, NULL);
  ASSERT(status == ndims);
  status = H5Sget_simple_extent_dims(xyz_space, xyz_shape, NULL);
  ASSERT(status == ndims);
  if ( proc_zero() ) {
    // [DRM-8] hsize_t printed with %llu (was %d)
    printf("> Group 'node_motion_data': dimensions %llu x %llu \n",
           (unsigned long long)node_shape[0], (unsigned long long)node_shape[1]);
    printf("> Group 'xyz': dimensions %llu x %llu \n",
           (unsigned long long)xyz_shape[0], (unsigned long long)xyz_shape[1]);
    fflush(stdout);
  }
  if (xyz_shape[1] != 3 || node_shape[0] != 3*xyz_shape[0]) {
    cerr << "Error: DRM HDF5 shape mismatch: expected xyz = N x 3 and "
         << "node_motion_data = 3N x nt" << endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  const hsize_t nnodes = xyz_shape[0];   // number of DRM nodes
  const hsize_t nt     = node_shape[1];  // number of time steps

  // Every rank reads the whole xyz table (small).
  xyz_data = (float_sw4 *)malloc(nnodes*xyz_shape[1]*sizeof(float_sw4));
  status = H5Dread(xyz_set, h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz_data);
  ASSERT(status >= 0);
  H5Dclose(xyz_set);
  H5Sclose(xyz_space);

  // The xyz coordinates of the drm nodes are saved with respect to the DRM
  // boundary origin, so we offset them into this simulation's coordinate frame
  // (same offset as the original code applied to x and y).
  const float_sw4 xyoff = m_supergrid_width + 6*mGridSize[0];

  // ----------------------------------------------------------------------
  // [DRM-4] Ownership scan: which DRM nodes lie in this rank's patch?
  //   Uses only positions (no time series). The test is the same one
  //   DRMSource::compute_grid_point uses to set myPoint().
  //   ASSUMPTION: computeNearestGridPoint2 is purely local (no MPI calls),
  //   which matches how DRMSource::compute_grid_point uses it.
  // ----------------------------------------------------------------------
  std::vector<hsize_t> mine;
  for ( hsize_t ix = 0; ix < nnodes; ix++ ){
    float_sw4 xs = xyz_data[ix*xyz_shape[1] + 0] + xyoff;
    float_sw4 ys = xyz_data[ix*xyz_shape[1] + 1] + xyoff;
    float_sw4 zs = xyz_data[ix*xyz_shape[1] + 2];
    int ii, jj, kk, gg;
    int success = computeNearestGridPoint2( ii, jj, kk, gg, xs, ys, zs );
    if ( success && interior_point_in_proc(ii, jj, gg) )
      mine.push_back(ix);
  }
  const size_t nloc = mine.size();

  // [DRM-5] Diagnostics: ownership sanity check and memory estimate.
  //   Each DRMSource stores 6*nt+1 spline coefficients per component, i.e.
  //   about 6x the raw time series, which is the dominant memory cost.
  {
    unsigned long long nloc_ull = nloc, nsum = 0, nmax = 0;
    MPI_Allreduce(&nloc_ull, &nsum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nloc_ull, &nmax, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    if ( proc_zero() ){
      // [SW-5] estimate uses the window length, not the whole series
      double nint_kept = (m_DRM_window > 0 && (hsize_t)m_DRM_window < nt) ? m_DRM_window : (double)nt;
      double bytes_per_node_src = 3.0*(6.0*nint_kept)*sizeof(float_sw4);
      printf("> DRM nodes owned: %llu of %llu total; busiest rank owns %llu\n",
             nsum, (unsigned long long)nnodes, nmax);
      printf("> Estimated DRM source storage: %.1f GB total, %.1f GB on the busiest rank\n",
             nsum*bytes_per_node_src/1.0e9, nmax*bytes_per_node_src/1.0e9);
      if ( nsum != (unsigned long long)nnodes )
        printf("Warning: %llu DRM nodes are owned by no rank (outside the domain?) "
               "or by more than one rank\n", (unsigned long long)nnodes - nsum);
      fflush(stdout);
    }
  }
  etime1 = MPI_Wtime();
  if ( proc_zero() ){
    printf("* Located DRM nodes, took %.8f seconds.\n", etime1-stime);
    printf("> Now generating DRM source objects and loading the first time window...\n");
    fflush(stdout);
  }

  // We set the number of time steps here while we have the information
  int e = 0;
  if (mTimeIsSet[e])
    mNumberOfTimeSteps[e] = nt;

  // ----------------------------------------------------------------------
  // [SW-5] REPLACES the [DRM-6] full-series batched read.
  //   Sources are created here WITHOUT data; this rank only records which file
  //   nodes it owns. The time series are then loaded one window at a time by
  //   refreshDRMWindow(), starting with the window that contains t = 0.
  //   Node ix occupies rows 3*ix .. 3*ix+2 of node_motion_data, and the sources
  //   are pushed in the same order (3 per node), which refreshDRMWindow relies on.
  // ----------------------------------------------------------------------
  if ( nt < 6 ){
    cerr << "Error: DRM node_motion_data has only " << nt << " samples (need >= 6)" << endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  m_DRM_nt   = (int)nt;
  m_DRM_file = dfile;
  m_DRM_myNodes.assign(mine.begin(), mine.end());
  m_DRM_winStart = m_DRM_winEnd = 0;   // nothing loaded yet

  m_DRMSources.reserve(m_DRMSources.size() + 3*nloc);
  for ( size_t kk2 = 0; kk2 < nloc; kk2++ ){
    hsize_t ix = mine[kk2];
    float_sw4 x = xyz_data[ix*xyz_shape[1] + 0] + xyoff;
    float_sw4 y = xyz_data[ix*xyz_shape[1] + 1] + xyoff;
    float_sw4 z = xyz_data[ix*xyz_shape[1] + 2];
    for ( int comp = 0; comp < 3; comp++ ){
      DRMSource* drmsrcPtr = new DRMSource(this, x, y, z, comp+1, (int)nt, m_DRM_dt);
      if (!drmsrcPtr->myPoint()){
        // cannot happen after [DRM-4]; abort rather than break the 3-per-node ordering
        cerr << "Error: DRM node " << (unsigned long long)ix << " selected but not owned" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
      }
      m_DRMSources.push_back(drmsrcPtr);
    }
  }

  // [DRM-7] The old code's MPI_Bcast calls are gone. They had three problems:
  //   * node_shape/xyz_shape (hsize_t) were broadcast as MPI_DOUBLE;
  //   * the node_motion_data broadcast count node_shape[0]*node_shape[1]
  //     (~7.6e10) does not fit in the int count argument of MPI_Bcast;
  //   * MPI_Comm_rank/size for read_rank/read_size were called on
  //     MPI_COMM_WORLD instead of read_comm (and the values were unused).
  //   Buffers are now released with free() to match malloc() (the old code
  //   used delete[] on malloc'd memory).
  free(xyz_data);
  H5Sclose(node_space);
  H5Dclose(node_set);
  H5Fclose(fid);   // collective (mpio) handle used for setup only

  // [SW-5] load the first window (independent, per-rank reads; no collectives)
  refreshDRMWindow(0);

  MPI_Barrier(MPI_COMM_WORLD);
  etime2 = MPI_Wtime();
  if ( proc_zero() ){
    // [DRM-8] removed stray extra printf arguments (xyz_shape[0] was passed
    // with no matching format specifier)
    printf("* Generated DRM source objects, took %.8f seconds.\n", etime2-etime1);
    printf("* Total time to process DRM node motions: %.8f seconds.\n\n", etime2-stime);
    fflush(stdout);
  }
#endif // USE_HDF5
}

// ======================================================================
// [SW-6] NEW: sliding time window for the DRM node motions
// ======================================================================

// Loads, for every DRM source on this rank, the spline coefficients of the
// window that contains global interval a_kNeed.
//
// Windows are aligned to multiples of W = m_DRM_window intervals (so that a
// time-chunked file is read chunk by chunk) and extended by 2 intervals on each
// side, so that intervals a_kNeed-1 .. a_kNeed+1 are always inside the window
// that is loaded. The spline is fitted on the window's samples plus
// m_DRM_overlap extra samples on each side.
//
// Global sample convention (unchanged from the original DRMSource): g[0] = 0
// at t = 0, and g[m] = node_motion_data column m-1 for m = 1..nt; interval i
// lies between samples i and i+1, so there are nt intervals.
//
// Purely local: no MPI collectives, so it is safe for each rank to call it
// whenever it needs to. Ranks without DRM sources never open the file.
//-----------------------------------------------------------------------
void EW::refreshDRMWindow( int a_kNeed )
{
#ifdef USE_HDF5
  if ( m_DRMSources.empty() )
    return;

  const int nInt = m_DRM_nt;            // global number of intervals
  const int W = ( m_DRM_window > 0 && m_DRM_window < nInt ) ? m_DRM_window : nInt;
  if ( a_kNeed < 0 )      a_kNeed = 0;
  if ( a_kNeed > nInt-1 ) a_kNeed = nInt-1;

  const int block    = (a_kNeed / W) * W;
  const int intStart = std::max(0,    block - 2);
  const int intEnd   = std::min(nInt, block + W + 2);

  // global samples to fit: [s0, s1]
  int s0 = std::max(0,        intStart - m_DRM_overlap);
  int s1 = std::min(m_DRM_nt, intEnd   + m_DRM_overlap);
  while ( s1 - s0 + 1 < 7 ){             // Qspline needs >= 7 points
    if ( s0 > 0 ) s0--; else if ( s1 < m_DRM_nt ) s1++; else break;
  }
  const int nsamp = s1 - s0 + 1;
  const int m0    = std::max(s0, 1);     // first sample that comes from the file
  const int lead  = m0 - s0;             // 1 if the synthetic g[0] = 0 is included
  const hsize_t c0   = (hsize_t)(m0 - 1);
  const hsize_t ncol = (hsize_t)(s1 - m0 + 1);

  // open lazily (per-rank, read-only, default driver: independent I/O)
  if ( m_DRM_fid < 0 ){
    m_DRM_fid = H5Fopen(m_DRM_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( m_DRM_fid < 0 ){
      cerr << "Error: rank " << m_myRank << " could not reopen DRM file " << m_DRM_file << endl;
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    m_DRM_dset = H5Dopen(m_DRM_fid, "node_motion_data", H5P_DEFAULT);
    if ( m_DRM_dset < 0 ){
      cerr << "Error: no 'node_motion_data' in DRM file " << m_DRM_file << endl;
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }
  }
  const hid_t h5type = (sizeof(float_sw4) == 8) ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;
  hid_t fspace = H5Dget_space(m_DRM_dset);

  // read in batches of nodes (<= ~256 MB), merging consecutive node indices
  const size_t nloc = m_DRM_myNodes.size();
  const size_t bytes_per_node = 3*(size_t)ncol*sizeof(float_sw4);
  size_t batch_nodes = (256ULL*1024ULL*1024ULL) / bytes_per_node;
  if ( batch_nodes < 1 )    batch_nodes = 1;
  if ( batch_nodes > nloc ) batch_nodes = nloc;

  std::vector<float_sw4> buf(batch_nodes*3*(size_t)ncol);
  std::vector<float_sw4> samples(nsamp, 0.0);

  for ( size_t b0 = 0; b0 < nloc; b0 += batch_nodes ){
    size_t b1 = std::min(nloc, b0 + batch_nodes);
    size_t k = b0;
    while ( k < b1 ){
      size_t run = 1;
      while ( k+run < b1 && m_DRM_myNodes[k+run] == m_DRM_myNodes[k] + run )
        run++;
      hsize_t start[2] = { 3*m_DRM_myNodes[k], c0 };
      hsize_t count[2] = { 3*(hsize_t)run, ncol };
      herr_t status = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
      ASSERT(status >= 0);
      hid_t mspace = H5Screate_simple(2, count, NULL);
      status = H5Dread(m_DRM_dset, h5type, mspace, fspace, H5P_DEFAULT,
                       &buf[(k-b0)*3*(size_t)ncol]);
      ASSERT(status >= 0);
      H5Sclose(mspace);
      k += run;
    }
    // fit each source's window; sources 3*k..3*k+2 belong to owned node k
    for ( size_t kk = b0; kk < b1; kk++ )
      for ( int comp = 0; comp < 3; comp++ ){
        const float_sw4* row = &buf[((kk-b0)*3 + comp)*(size_t)ncol];
        for ( int s = 0; s < lead; s++ ) samples[s] = 0.0;   // g[0] = 0
        for ( hsize_t c = 0; c < ncol; c++ ) samples[lead + c] = row[c];
        m_DRMSources[3*kk + comp]->setWindow( &samples[0], nsamp, s0, intStart, intEnd );
      }
  }
  H5Sclose(fspace);

  m_DRM_winStart = intStart;
  m_DRM_winEnd   = intEnd;
  if ( mVerbose >= 3 )
    printf("> rank %d: DRM window intervals [%d,%d) loaded (samples %d..%d)\n",
           m_myRank, intStart, intEnd, s0, s1);
#endif // USE_HDF5
}

// Closes this rank's handle on the DRM file. refreshDRMWindow reopens it if a
// later solve (e.g. another EQL iteration) needs a different window.
//-----------------------------------------------------------------------
void EW::closeDRMWindowFile()
{
#ifdef USE_HDF5
  if ( m_DRM_dset >= 0 ) { H5Dclose(m_DRM_dset); m_DRM_dset = -1; }
  if ( m_DRM_fid  >= 0 ) { H5Fclose(m_DRM_fid);  m_DRM_fid  = -1; }
#endif
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

  // [SW-7] NEW: make sure the loaded time window covers t. getU(t) uses interval
  //   k = floor(t/dt); we require k-1..k+1 to be loaded (slack for rounding), and
  //   otherwise load the window containing k. This also handles jumps backwards
  //   in time (a restart, or the next EQL iteration starting again at t = 0).
  if ( !m_DRMSources.empty() ){
    int kNeed = static_cast<int>( floor( t/m_DRM_dt ) );
    if ( kNeed < 0 )            kNeed = 0;
    if ( kNeed > m_DRM_nt-1 )   kNeed = m_DRM_nt-1;
    int kLo = std::max(0, kNeed-1);
    int kHi = std::min(m_DRM_nt-1, kNeed+1);
    if ( kLo < m_DRM_winStart || kHi >= m_DRM_winEnd )
      refreshDRMWindow( kNeed );
  }

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