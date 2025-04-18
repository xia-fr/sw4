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
// Constructor,
//  x0, y0, z0 : float_sw4 : The cartesian coordinate location of this source.
//  comp : float_sw4 : The component this source represents. 1:ux, 2:uy, 3:uz.
//
//  pars : pars[1],..pars[npts] should contain the discrete function on a uniform
//         grid with spacing dt=1/freq, and pars[0] is the first time, thus the 
//         grid is t_k = pars[0] + dt*k, k=0,1,..,npts-1
//-----------------------------------------------------------------------
DRMSource::DRMSource(EW *a_ew,
              float_sw4  x0, float_sw4  y0, float_sw4  z0, int comp,
              float_sw4 *pars, int npts, float_sw4 dt):
  mX0(x0), mY0(y0), mZ0(z0), mComp(comp), mNpts(npts), mDt(dt), m_myPoint(false)
{
  // We have to do some adjusting because in the process of exporting as a ssi file
  // and then the extracting of node motions into an hdf5, a 0th time step gets
  // chopped off somewhere which misaligns these node motions with the simulation
  // by one. And we need the timing to match for the in-domain forces.
  mPar = new float_sw4[mNpts + 3];
  mPar[0] = 0.0; // mPar[0] holds the t0 value (assumed to be 0 by default)
  mPar[1] = 0.0; // this is the value of the time series at t = t0
  for( int i = 0 ; i < mNpts+1; i++ )
    mPar[i+2] = pars[i]; // Save the given time series
  mNpts = mNpts+1;

  // if ( (a_ew->getRank() == 5) && (mX0 == 2000.0) && (mY0 == 2000.0) && (mZ0 == 1000.00) && (mComp == 1) ){
  //   printf("rank %i:\n", a_ew->getRank());
  //   for(int ii=95;ii<105;ii++){
  //     printf("[%i]: %.8e\n", ii, mPar[ii]);
  //   }
  // }

  // Apply spline interpolation to the given time series
  spline_interpolation();

  // Correct source location for discrepancy between raw and smoothed topography
  // also sets the ignore flag for sources that are above the topography
  // correct_Z_level( a_ew ); 
  
  // Compute the ijk grid indices and the corresponding grid g for this point
  compute_grid_point( a_ew );
}

//-----------------------------------------------------------------------
DRMSource::DRMSource()
{
   mNpts = 0;
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
  output << "   t0 = " << s.mPar[0] << " dt = " << s.mDt << " npts = " << s.mNpts << endl;
  return output;
}

// Returns the value of the source function at time t.
//-----------------------------------------------------------------------
float_sw4 DRMSource::getU( float_sw4 t, float_sw4 ewDt )
{
  float_sw4 tstart = mPar[0];
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
  //std::cout <<  "t= " << t << " mNpts " << mNpts << " k= " << k << "arg = " << arg <<  std::endl;
  return mPar[6*k+1] + mPar[2+6*k]*arg + mPar[3+6*k]*arg*arg + mPar[4+6*k]*arg*arg*arg +
          mPar[5+6*k]*arg*arg*arg*arg + mPar[6+6*k]*arg*arg*arg*arg*arg; 
}

// Trimmed down version of spline_interpolation() from Source.C
//-----------------------------------------------------------------------
int DRMSource::spline_interpolation( )
{
  // Assume mPar[1], to mPar[npts] contain the function
  // Assume mPar[0] is tstart.
  // Compute the six spline coefficients for each interval and return in mPar[1],to mPar[6*(npts-1)]
  Qspline quinticspline( mNpts, &mPar[1], mPar[0], mDt );
  float_sw4 tstart = mPar[0];
  delete[] mPar;

  int mNpar = 6*(mNpts-1)+1;
  mPar = new float_sw4[mNpar];
  
  mPar[0] = tstart;
  float_sw4* qsppt = quinticspline.get_polycof_ptr();
  for( int i=0 ; i < 6*(mNpts-1) ; i++ )
      mPar[i+1] = qsppt[i];

  return 1;
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

  int inds[4]={-9999,-9999,-9999,-9999}; // if this is not my point, i0,j0,k0=-9999
  if( m_myPoint )
  {
    inds[0] = i;
    inds[1] = j;
    inds[2] = k;
    inds[3] = g;
  }
  int indsg[4]={0,0,0,0};
  MPI_Allreduce( inds, indsg, 4, MPI_INT, MPI_MAX, a_ew->m_1d_communicator );
  m_i0 = indsg[0];
  m_j0 = indsg[1];
  m_k0 = indsg[2];
  m_grid = indsg[3];

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