//-*-c++-*-
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
#ifndef DRMSRC_H
#define DRMSRC_H

#include "sw4.h"

#include <iostream>
#include <vector>
#include <string>

class EW;

class DRMSource{
  friend std::ostream& operator<<(std::ostream& output, const DRMSource& s);

  public:

  // [SW-1] CHANGED: the constructor no longer takes the time series. It only
  //   records where the source is and finds its grid point. The time series is
  //   supplied, one window at a time, through setWindow().
  //   a_nSamples = number of samples in node_motion_data (nt); one leading zero
  //   sample at t = 0 is added internally, as the old constructor did.
  DRMSource(EW * a_ew,
    float_sw4  x0, float_sw4  y0, float_sw4  z0, int comp,
    int a_nSamples, float_sw4 dt);
  ~DRMSource();

  int m_i0, m_j0, m_k0;
  int m_grid;

  float_sw4 getX0() const;
  float_sw4 getY0() const;
  float_sw4 getZ0() const;
  float_sw4 getU( float_sw4 t, float_sw4 ewDt );
  int getComp() const;

  bool myPoint(){ return m_myPoint; }

  // [SW-1] NEW: load spline coefficients for global intervals
  //   [a_intStart, a_intEnd), fitted on the global samples
  //   a_samples[0..a_nSamp-1] = g[a_firstSample .. a_firstSample+a_nSamp-1],
  //   where g[0] = 0 and g[m] = file column m-1 for m >= 1.
  void setWindow( const float_sw4* a_samples, int a_nSamp, int a_firstSample,
                  int a_intStart, int a_intEnd );

  // [SW-1] NEW: total number of spline intervals (global), = number of samples incl. g[0] minus 1
  int numIntervals() const { return mNpts-1; }

  private:
  DRMSource();

  // void correct_Z_level( EW *a_ew );
  void compute_grid_point( EW *a_ew );

  float_sw4 mX0,mY0,mZ0;
  float_sw4 mDt;
  float_sw4 mT0;     // [SW-1] start time of the series (0, as before)
  float_sw4* mPar;   // [SW-1] now holds 6 coefficients per interval of the CURRENT window only
  int mNpts, mComp;  // mNpts = total global samples incl. g[0] (= nt+1, as before)
  int m_kOffset;     // [SW-1] global index of the first interval held in mPar
  int m_nInt;        // [SW-1] number of intervals held in mPar

  bool m_myPoint;
};
#endif