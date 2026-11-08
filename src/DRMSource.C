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

#include "EW.h"
#include "DRMSource.h"
#include "mpi.h"

#include "Qspline.h"

using namespace std;

//-----------------------------------------------------------------------
// [SW-2] CHANGED constructor.
//  x0, y0, z0 : float_sw4 : The cartesian coordinate location of this source.
//  comp : int : The component this source represents. 1:ux, 2:uy, 3:uz.
//  a_nSamples : number of samples per series in node_motion_data (nt).
//  dt : sample spacing of the series.
//
//  OLD: copied the whole series, prepended a zero sample (the "chopped-off 0th
//       time step"), fitted one quintic spline over all nt+1 samples, and kept
//       6*nt coefficients for the whole run.
//  NEW: stores no data. The global sample sequence is still defined the same way,
//       g[0] = 0 and g[m] = series[m-1] for m = 1..nt (mNpts = nt+1 samples,
//       nt intervals), but coefficients are loaded per window by setWindow().
//-----------------------------------------------------------------------
DRMSource::DRMSource(EW *a_ew,
              float_sw4  x0, float_sw4  y0, float_sw4  z0, int comp,
              int a_nSamples, float_sw4 dt):
  mX0(x0), mY0(y0), mZ0(z0), mDt(dt), mT0(0.0), mPar(NULL),
  mNpts(a_nSamples+1), mComp(comp), m_kOffset(0), m_nInt(0), m_myPoint(false)
{
  // Compute the ijk grid indices and the corresponding grid g for this point
  compute_grid_point( a_ew );
}

//-----------------------------------------------------------------------
DRMSource::DRMSource()
{
   mNpts = 0;
   // [DRMS-3] BUG FIX: mPar (and m_myPoint) were left uninitialised, so the destructor's
   //   delete[] mPar was undefined behaviour for default-constructed objects.
   mPar = NULL;
   m_myPoint = false;
   m_kOffset = 0;   // [SW-2]
   m_nInt = 0;      // [SW-2]
   mT0 = 0.0;       // [SW-2]
}

//-----------------------------------------------------------------------
DRMSource::~DRMSource()
{
  delete[] mPar;
}

// 'Get' functions
//-----------------------------------------------------------------------
float_sw4 DRMSource::getX0() const {return mX0;}
float_sw4 DRMSource::getY0() const {return mY0;}
float_sw4 DRMSource::getZ0() const {return mZ0;}
int DRMSource::getComp() const {return mComp;}

// Ostream return information function
//-----------------------------------------------------------------------
ostream& operator<<( ostream& output, const DRMSource& s )
{
  output << "DRM Location (X,Y,Z) = " << s.mX0 << "," << s.mY0 << "," << s.mZ0 <<
    " at grid point " << s.m_i0 << " " << s.m_j0 << " " << s.m_k0 <<
    " in grid no " << s.m_grid  << endl;
  output << "   Component: " << s.mComp << endl;
  // [SW-2] mPar[0] no longer holds t0; print the window instead
  output << "   t0 = " << s.mT0 << " dt = " << s.mDt << " npts = " << s.mNpts
         << " window intervals = [" << s.m_kOffset << "," << s.m_kOffset+s.m_nInt << ")" << endl;
  return output;
}

// Returns the value of the source function at time t.
// [SW-2] CHANGED: same interval index k and argument as before (same clamping),
//   but k is translated into the current window. The caller (EW::DRMSourcesToUField)
//   guarantees the window contains k; a miss is a logic error and aborts.
//-----------------------------------------------------------------------
float_sw4 DRMSource::getU( float_sw4 t, float_sw4 ewDt )
{
  float_sw4 tstart = mT0;
  float_sw4 freq = 1/mDt;
  int k = static_cast<int>(floor((t-tstart)*freq));

  if( k < 0 ){
    k = 0;
    t = tstart;
  }
  if( k > mNpts-2 ){
    k = mNpts-2;
    t = tstart+(mNpts-1)/freq;
  }

  float_sw4 arg = (t-tstart)*freq-k; // (t-(tstart+k*dt))/dt

  int kl = k - m_kOffset;   // [SW-2] index within the window
  if( kl < 0 || kl >= m_nInt ){
    cerr << "Error: DRMSource::getU interval " << k << " is outside the loaded window ["
         << m_kOffset << "," << m_kOffset+m_nInt << ") at t=" << t << endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  const float_sw4* c = &mPar[6*kl];
  return c[0] + c[1]*arg + c[2]*arg*arg + c[3]*arg*arg*arg +
         c[4]*arg*arg*arg*arg + c[5]*arg*arg*arg*arg*arg;
}

// [SW-2] NEW, replaces spline_interpolation().
//   Fits the same quintic spline (Qspline, 'not-a-knot' ends) as before, but on
//   the samples of one window plus overlap, and keeps only the coefficients of
//   intervals [a_intStart, a_intEnd). Because the spline interpolates the samples
//   exactly at every knot (and each interval ends exactly on the next sample),
//   values at on-grid times are identical to the full-series spline up to rounding;
//   the overlap only matters for off-grid times, where the windowed and global
//   splines differ by an amount that decays quickly with distance from the window edge.
//-----------------------------------------------------------------------
void DRMSource::setWindow( const float_sw4* a_samples, int a_nSamp, int a_firstSample,
                           int a_intStart, int a_intEnd )
{
  // the fitted samples must cover every interval we keep
  if( a_intStart < a_firstSample || a_intEnd > a_firstSample + a_nSamp - 1 || a_nSamp < 7 ){
    cerr << "Error: DRMSource::setWindow bad window: samples [" << a_firstSample << ","
         << a_firstSample+a_nSamp << ") intervals [" << a_intStart << "," << a_intEnd << ")" << endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  // Qspline takes a non-const pointer but does not modify fun
  Qspline quinticspline( a_nSamp, const_cast<float_sw4*>(a_samples),
                         mT0 + a_firstSample*mDt, mDt );
  float_sw4* qsppt = quinticspline.get_polycof_ptr();

  int nInt = a_intEnd - a_intStart;
  if( nInt != m_nInt ){
    delete[] mPar;
    mPar = new float_sw4[6*nInt];
  }
  int first = a_intStart - a_firstSample;   // first kept interval within the fit
  for( int i = 0 ; i < 6*nInt ; i++ )
    mPar[i] = qsppt[6*first + i];
  m_kOffset = a_intStart;
  m_nInt = nInt;
}

// From Source.C
//-----------------------------------------------------------------------
void DRMSource::compute_grid_point( EW* a_ew )
{
  // Sets values (m_i0,m_j0,m_k0) and m_grid.
  // Should be called after the topographic correction of mZ0
  int i, j, k, g;
  int success = a_ew->computeNearestGridPoint2( i, j, k, g, mX0, mY0, mZ0 );
  m_myPoint = success && a_ew->interior_point_in_proc(i, j, g);

  // [DRMS-1] REQUIRED CHANGE: the MPI_Allreduce that used to be here is removed.
  //   It was a collective over m_1d_communicator executed once per DRMSource, so every rank
  //   had to construct every source in the same order (3.9M Allreduces, each preceded by a
  //   full spline fit on every rank, only for non-owners to delete the object afterwards).
  //   readDRMHDF5_data now constructs sources only on the rank that owns the node [DRM-1],
  //   so keeping a collective here would hang. The owner already knows its own indices;
  //   on a non-owner (which never keeps the object) they are set to -9999.
  //   CHECK: grep the rest of SW4 for other 'new DRMSource' call sites that might rely on
  //   the indices being broadcast to non-owning ranks.
  if( m_myPoint )
  {
    m_i0 = i;
    m_j0 = j;
    m_k0 = k;
    m_grid = g;
  }
  else
  {
    m_i0 = m_j0 = m_k0 = m_grid = -9999;
  }

  // if ( (mX0 == 2000.0) && (mY0 == 2000.0) && (mZ0 == 1000.00) && mComp == 1 ){
  //   int rank = a_ew->getRank();
  //   printf("(after)myrank: %i\n", rank);
  //   printf("[i0]: %i\n", m_i0);
  //   printf("[j0]: %i\n", m_j0);
  //   printf("[k0]: %i\n", m_k0);
  // }

   //   std::cout << a_ew->getRank() << " mypoint = " << m_myPoint 
   //             << " (i,j,k)= " << m_i0 <<" " << m_j0 << " " << m_k0 << " grid= "
   //             << m_grid << std::endl;
}