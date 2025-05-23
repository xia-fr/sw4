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
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License"
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
#include "ESSI3DHDF5.h"
#include "EW.h"
#include "mpi.h"

#ifdef USE_ZFP
#include "H5Zzfp_lib.h"
#include "H5Zzfp_props.h"
#endif

#ifdef USE_SZ
#include "H5Z_SZ.h"
#endif

using std::cerr;
using std::cout;

ESSI3DHDF5* ESSI3DHDF5::nil = static_cast<ESSI3DHDF5*>(0);

//-----------------------------------------------------------------------
ESSI3DHDF5::ESSI3DHDF5(const std::string& filename, int (&global)[3],
                       int (&window)[6], bool ihavearray, int precision)
    : m_end_cycle(-1),
      m_filename(filename),
      m_ihavearray(ihavearray),
      m_precision(precision) {
#ifdef USE_HDF5
  bool debug = false;
  /* debug = true; */
  for (int d = 0; d < 3; d++) {
    m_global[d] = global[d];
    m_window[2 * d] = window[2 * d];          // lo, relative to global
    m_window[2 * d + 1] = window[2 * d + 1];  // hi
  }
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if (debug) {
    char msg[1000];
    sprintf(msg, "Rank %d global size = %d %d %d\n", myRank, m_global[0],
            m_global[1], m_global[2]);
    cerr << msg;
    sprintf(msg, "Rank %d ihavearray = %s\n", myRank,
            (m_ihavearray) ? "yes" : "no");
    cerr << msg;
    sprintf(msg, "Rank %d index range = (%d:%d %d:%d %d:%d)\n", myRank,
            m_window[0], m_window[1], m_window[2], m_window[3], m_window[4],
            m_window[5]);
    cerr << msg;
  }

  // Time step is 1st comp
  m_slice_dims[0] = 1;   // write one time step
  m_global_dims[0] = 1;  // Extendible dimension for cycle
  m_window_dims[0] = 1;  // write one time step
  m_cycle_dims[0] = 1;   // this will become cycle-1

  for (int d = 1; d < 4; d++) {
    m_window_dims[d] = m_window[2 * (d - 1) + 1] - m_window[2 * (d - 1)] + 1;
    m_global_dims[d] = m_global[d - 1];
    m_cycle_dims[d] = m_global[d - 1];
  }

  // 1 i,j column per slice
  m_slice_dims[1] = 1;
  m_slice_dims[2] = 1;
  m_slice_dims[3] = m_global[2];

  m_file_id = 0;
  m_es_id = 0;
#endif
}

//-----------------------------------------------------------------------
ESSI3DHDF5::~ESSI3DHDF5() {}

//-----------------------------------------------------------------------
void ESSI3DHDF5::create_file(bool is_restart, bool is_root) {
#ifdef USE_HDF5
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_libver_bounds(fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

  int m_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
  if (is_restart) {
    if (!is_root) {
      H5Pset_fapl_mpio(fapl, comm, info);
      H5Pset_coll_metadata_write(fapl, 1);
      H5Pset_all_coll_metadata_ops(fapl, 1);
    }

#ifdef USE_HDF5_ASYNC
    size_t num_in_progress;
    hbool_t op_failed;

    if (m_es_id == 0)
      m_es_id = H5EScreate();
    else
      H5ESwait(m_es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);

    m_file_id = H5Fopen_async(const_cast<char*>(m_filename.c_str()),
                              H5F_ACC_RDWR, fapl, m_es_id);
#else
    m_file_id =
        H5Fopen(const_cast<char*>(m_filename.c_str()), H5F_ACC_RDWR, fapl);
#endif
  } 
  else {
    m_file_id = H5Fcreate(const_cast<char*>(m_filename.c_str()), H5F_ACC_TRUNC,
                          H5P_DEFAULT, fapl);
  }

  if (m_file_id <= 0) {
    cerr << "Could not open hdf5 file: " << m_filename << endl;
    MPI_Abort(comm, m_file_id);
  }
  H5Pclose(fapl);
#endif
  return;
}

void ESSI3DHDF5::write_header(double h, double (&lonlat_origin)[2], double az,
                              double (&origin)[3], int cycle, double t,
                              double dt) {
#ifdef USE_HDF5
  bool debug = false;
  /* debug=true; */
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  // Add all the metadata to the file
  if (debug && (myRank == 0))
    cerr << "Writing hdf5 metadata for file: " << m_filename << endl;

  // Global grid lat,lon origin
  hsize_t dim = 2;
  hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
  hid_t dataset_id =
      H5Dcreate2(m_file_id, "Grid lon-lat origin", H5T_NATIVE_DOUBLE,
                 dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t ierr = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                         H5P_DEFAULT, &lonlat_origin);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // Grid azimuth
  dim = 1;
  dataspace_id = H5Screate_simple(1, &dim, NULL);
  dataset_id = H5Dcreate2(m_file_id, "Grid azimuth", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ierr = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &az);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // grid spacing
  dim = 1;
  dataspace_id = H5Screate_simple(1, &dim, NULL);
  dataset_id = H5Dcreate2(m_file_id, "ESSI xyz grid spacing", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ierr = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &h);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // x,y origin
  dim = 3;
  dataspace_id = H5Screate_simple(1, &dim, NULL);
  dataset_id = H5Dcreate2(m_file_id, "ESSI xyz origin", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ierr = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &origin);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // t start
  dim = 1;
  dataspace_id = H5Screate_simple(1, &dim, NULL);
  dataset_id = H5Dcreate2(m_file_id, "time start", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ierr = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &t);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // dt
  dim = 1;
  dataspace_id = H5Screate_simple(1, &dim, NULL);
  dataset_id = H5Dcreate2(m_file_id, "timestep", H5T_NATIVE_DOUBLE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ierr = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &dt);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  if (debug && (myRank == 0)) cerr << "Writing hdf5 metadata done" << endl;
#endif
  return;
}

void ESSI3DHDF5::write_topo(void* window_array) {
#ifdef USE_HDF5
  bool debug = false;
  /* debug=true; */
  MPI_Comm comm = MPI_COMM_WORLD;
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  hid_t dtype = H5T_NATIVE_DOUBLE;
  if (m_precision == 4) dtype = H5T_NATIVE_FLOAT;

  if (debug && (myRank == 0))
    cerr << "Writing hdf5 z coordinate: " << m_filename << endl;

  // Create the data space and chunk for z coordinate
  char msg[1000];
  if (debug) {
    sprintf(msg, "Rank %d creating z dataspace size = %d %d %d\n", myRank,
            (int)m_global_dims[0], (int)m_global_dims[1],
            (int)m_global_dims[2]);
    cerr << msg;
  }

  hsize_t z_dims = 3;
  // Modify dataset creation properties to enable chunking
  hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
  herr_t ierr;
  /* ierr = H5Pset_chunk (prop_id, z_dims, m_slice_dims); */

  // Create dataset for z field
  hid_t dataspace_id = H5Screate_simple(z_dims, &m_global_dims[1], NULL);
  hid_t dataset_id = H5Dcreate2(m_file_id, "z coordinates", dtype, dataspace_id,
                                H5P_DEFAULT, prop_id, H5P_DEFAULT);
  ierr = H5Pclose(prop_id);

  if (debug) {
    sprintf(msg, "Rank %d creating chunk/slice size = %d %d %d\n", myRank,
            (int)m_slice_dims[1], (int)m_slice_dims[2], (int)m_slice_dims[3]);
    cerr << msg;
  }

  /*
  if (debug)
     cout << "Rank " << myRank <<
       " in write_topo_real: starting slab loop" << endl;
  sprintf(msg, "Rank %d looping over hyperslab = (%d:%d %d:%d %d:%d)\n",
      myRank, m_window[0], m_window[1], m_window[2], m_window[3],
      m_window[4], m_window[5]);
  cout << msg;
  */

  hsize_t start[3] = {0, 0, 0};
  start[0] = m_window[0];  // i window offset
  start[1] = m_window[2];  // j window offset
  start[2] = m_window[4];  // k local index lo relative to global
  // We should write data
  if (debug) {
    sprintf(msg, "Rank %d selecting z hyperslab = %d %d %d\n", myRank,
            (int)start[0], (int)start[1], (int)start[2]);
    cerr << msg;
  }
  hid_t window_id = H5Screate_simple(z_dims, &m_window_dims[1], NULL);
  ierr = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL,
                             &m_window_dims[1], NULL);
  if (ierr < 0) {
    cerr << "Error from z H5Sselect_hyperslab" << endl;
    MPI_Abort(comm, ierr);
  }
  if (debug) {
    char msg[1000];
    sprintf(msg, "Writing z array Rank %d\n", myRank);
    cerr << msg;
  }

  ierr = H5Dwrite(dataset_id, dtype, window_id, dataspace_id, H5P_DEFAULT,
                  window_array);
  if (ierr < 0) {
    cerr << "Error from z H5Dwrite " << endl;
    MPI_Abort(comm, ierr);
  }

  ierr = H5Sclose(dataspace_id);
  ierr = H5Dclose(dataset_id);

  H5Fflush(m_file_id, H5F_SCOPE_LOCAL);

  if (debug && (myRank == 0))
    cerr << "Done writing hdf5 z coordinate: " << m_filename << endl;
#endif
  return;
}
void ESSI3DHDF5::init_write_vel(bool isRestart, int ntimestep,
                                int compressionMode, double compressionPar,
                                int bufferInterval) {
  bool debug = false;
  /* debug=true; */

#ifdef USE_HDF5
  int myRank, nProc;

  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);

  hid_t dtype = H5T_NATIVE_DOUBLE;
  if (m_precision == 4) dtype = H5T_NATIVE_FLOAT;

  // Create the extendible velocity data space and chunk
  hsize_t num_dims = 4;
  // m_xvel_dataspace_id = H5Screate_simple(num_dims, m_global_dims, dims);
  hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_alloc_time(prop_id, H5D_ALLOC_TIME_LATE);
  H5Pset_fill_time(prop_id, H5D_FILL_TIME_NEVER);

  if (ntimestep > 0)
    m_cycle_dims[0] = ntimestep;
  else
    printf("Error with m_ntimestep=%d!\n", ntimestep);

  hsize_t my_chunk[4] = {0, 0, 0, 0};

  if (compressionMode > 0) {
    MPI_Allreduce(&m_window_dims[1], &my_chunk[1], 3, MPI_UNSIGNED_LONG_LONG,
                  MPI_MAX, MPI_COMM_WORLD);

    // Keep last dimension (z)
    for (int i = 1; i < 3; i++) {
      if (my_chunk[i] > 4) my_chunk[i] = ((hsize_t)(my_chunk[i] / 4)) * 4;
    }

    my_chunk[0] = bufferInterval;

    char* env_char = NULL;
    env_char = getenv("SSI_CHUNK_X");
    if (env_char != NULL) my_chunk[1] = atoi(env_char);

    env_char = getenv("SSI_CHUNK_Y");
    if (env_char != NULL) my_chunk[2] = atoi(env_char);

    env_char = getenv("SSI_CHUNK_Z");
    if (env_char != NULL) my_chunk[3] = atoi(env_char);

    hsize_t total_chunk_size = m_precision;
    for (int i = 0; i < num_dims; i++) total_chunk_size *= my_chunk[i];

    // Make sure chunk size is less than 4GB, which is the HDF5 chunk limit
    while (total_chunk_size >= 4294967295llu) {
      if (my_chunk[1] > my_chunk[2])
        my_chunk[1] /= 2;
      else
        my_chunk[2] /= 2;
      total_chunk_size /= 2;
    }

    H5Pset_chunk(prop_id, num_dims, my_chunk);

    if (myRank == 0) {
      /* if (debug && myRank == 0) { */
      printf("SSI output chunk sizes:");
      for (int i = 0; i < num_dims; i++) printf("%llu  ", my_chunk[i]);
      printf("\n");
      fflush(stdout);
    }

    if (compressionMode == SW4_SZIP) {
      H5Pset_szip(prop_id, H5_SZIP_NN_OPTION_MASK, 32);
    } else if (compressionMode == SW4_ZLIB) {
      H5Pset_deflate(prop_id, (int)compressionPar);
    }
#ifdef USE_ZFP
    else if (compressionMode == SW4_ZFP_MODE_RATE) {
      H5Pset_zfp_rate(prop_id, compressionPar);
    } else if (compressionMode == SW4_ZFP_MODE_PRECISION) {
      H5Pset_zfp_precision(prop_id, (unsigned int)compressionPar);
    } else if (compressionMode == SW4_ZFP_MODE_ACCURACY) {
      H5Pset_zfp_accuracy(prop_id, compressionPar);
    } else if (compressionMode == SW4_ZFP_MODE_REVERSIBLE) {
      H5Pset_zfp_reversible(prop_id);
    }
#endif
#ifdef USE_SZ
    else if (compressionMode == SW4_SZ) {
      size_t cd_nelmts;
      unsigned int* cd_values = NULL;
      int dataType = SZ_DOUBLE;
      if (m_precision == 4) dataType = SZ_FLOAT;
      SZ_metaDataToCdArray(&cd_nelmts, &cd_values, dataType, 0, m_cycle_dims[3],
                           m_cycle_dims[2], m_cycle_dims[1], m_cycle_dims[0]);
      H5Pset_filter(prop_id, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts,
                    cd_values);
    }
#endif
  }

  if (debug && myRank == 0) {
    char msg[1000];
    sprintf(msg, "Rank %d creating vel dataspaces size = %d %d %d %d\n", myRank,
            (int)m_cycle_dims[0], (int)m_cycle_dims[1], (int)m_cycle_dims[2],
            (int)m_cycle_dims[3]);
    cerr << msg;
    fflush(stdout);
  }

  hid_t dset, dspace;
  if (myRank == 0) {
    for (int c = 0; c < 3; c++) {
      dspace = H5Screate_simple(num_dims, m_cycle_dims, m_cycle_dims);
      char var[100];
      sprintf(var, "vel_%d ijk layout", c);

      if (isRestart) {
          dset = H5Dopen(m_file_id, var, H5P_DEFAULT);
      }
      else {
          dset = H5Dcreate2(m_file_id, var, dtype, dspace, H5P_DEFAULT, prop_id,
                            H5P_DEFAULT);
      }
      if (dset < 0) {
        cerr << "Error with H5Dcreate/open!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
      H5Dclose(dset);
    }
    H5Fclose(m_file_id);
    m_file_id = 0;
  }
  H5Pclose(prop_id);
#endif
  return;
}

void ESSI3DHDF5::write_vel(void* window_array, int comp, int cycle, int nstep) {
  bool enable_timing = true;
  bool debug = false;
  /* debug=true; */
#ifdef USE_HDF5

  herr_t ierr;
  double write_time_start, write_time;
  int myRank;
  int write_size = m_precision;
  m_end_cycle = cycle;  // save for header for later when we close the file
  time_t now;

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  hid_t dtype = H5T_NATIVE_DOUBLE;
  if (m_precision == 4) dtype = H5T_NATIVE_FLOAT;

  for (int i = 0; i < 4; i++) write_size *= m_window_dims[i];

  if (enable_timing) write_time_start = MPI_Wtime();

  hsize_t vel_dims = 4;
  hsize_t buf_window_dims[4];
  buf_window_dims[0] = nstep;
  buf_window_dims[1] = m_window_dims[1];
  buf_window_dims[2] = m_window_dims[2];
  buf_window_dims[3] = m_window_dims[3];

  hid_t window_id = H5Screate_simple(vel_dims, buf_window_dims, NULL);

  hsize_t start[4] = {0, 0, 0, 0};
  start[0] = cycle - nstep;
  start[1] = m_window[0];  // local index lo relative to global
  start[2] = m_window[2];  // local index lo relative to global
  start[3] = m_window[4];  // local index lo relative to global

  if (debug && comp == 2) {
    time(&now);
    printf(
        "Rank %d: cycle=%d, nstep=%d, writing vel start = %d %d %d %d, size = "
        "%d %d %d %d, time: %s",
        myRank, cycle, nstep, (int)start[0], (int)start[1], (int)start[2],
        (int)start[3], (int)buf_window_dims[0], (int)buf_window_dims[1],
        (int)buf_window_dims[2], (int)buf_window_dims[3], ctime(&now));
    fflush(stdout);
  }

  hid_t dset, dspace;
  hsize_t my_size = 1;
  hsize_t num_dims = 4;
  for (int i = 0; i < vel_dims; i++) my_size *= m_window_dims[i];

  dspace = H5Screate_simple(num_dims, m_cycle_dims, m_cycle_dims);

  if (my_size == 0) {
    H5Sselect_none(dspace);
    H5Sselect_none(window_id);
  } else {
    ierr = H5Sselect_hyperslab(dspace, H5S_SELECT_SET, start, NULL,
                               buf_window_dims, NULL);
    if (ierr < 0) {
      cerr << "Error from vel H5Sselect_hyperslab" << endl;
      MPI_Abort(comm, ierr);
    }
  }

  char var[100];
  sprintf(var, "vel_%d ijk layout", comp);
#ifdef USE_HDF5_ASYNC
  dset = H5Dopen_async(m_file_id, var, H5P_DEFAULT, m_es_id);
#else
  dset = H5Dopen(m_file_id, var, H5P_DEFAULT);
#endif
  if (dset < 0) {
    cerr << "Error from vel H5Dopen " << myRank << endl;
    MPI_Abort(comm, ierr);
  }

#ifdef USE_HDF5_ASYNC
  ierr = H5Dwrite_async(dset, dtype, window_id, dspace, dxpl, window_array,
                        m_es_id);
#else
  ierr = H5Dwrite(dset, dtype, window_id, dspace, dxpl, window_array);
#endif
  if (ierr < 0) {
    cerr << "Error from vel H5Dwrite " << endl;
    MPI_Abort(comm, ierr);
  }

#ifdef USE_HDF5_ASYNC
  H5Dclose_async(dset, m_es_id);
#else
  H5Dclose(dset);
#endif

  if (enable_timing) {
    time(&now);
    write_time = MPI_Wtime() - write_time_start;
    if (myRank == 0 && buf_window_dims[0] >= 50) {
      printf(
          "ssioutput cycle=%d, comp=%d, write size=%.2fMB, write time=%f, %s",
          cycle, comp, write_size / 1048576.0, write_time, ctime(&now));
      fflush(stdout);
    }
  }

  H5Pclose(dxpl);
#ifndef USE_HDF5_ASYNC
  H5Fflush(m_file_id, H5F_SCOPE_GLOBAL);
#endif

#endif
  return;
}

void ESSI3DHDF5::close_file() {
#ifdef USE_HDF5
  // Updated header with final start,end cycle
  /* hsize_t dim = 2; */
  /* hid_t dataspace_id = H5Screate_simple(1, &dim, NULL); */
  /* hid_t dataset_id = */
  /*     H5Dcreate2(m_file_id, "cycle start, end", H5T_STD_I32LE, dataspace_id,
   */
  /*                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); */
  /* int cycles[2] = {0, m_end_cycle}; */
  /* herr_t ierr = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, */
  /*                        H5P_DEFAULT, cycles); */
  /* ierr = H5Dclose(dataset_id); */
  /* ierr = H5Sclose(dataspace_id); */

#ifdef USE_HDF5_ASYNC
  H5Fclose_async(m_file_id, m_es_id);
#else
  H5Fclose(m_file_id);
#endif
  m_file_id = 0;

#endif
  return;
}

void ESSI3DHDF5::finalize_hdf5() {
#ifdef USE_HDF5_ASYNC
  size_t num_in_progress;
  hbool_t op_failed;
  if (m_es_id > 0)
    H5ESwait(m_es_id, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
#endif
  return;
}
