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

  DRMSource(EW * a_ew,
    float_sw4  x0, float_sw4  y0, float_sw4  z0, int comp,
    float_sw4* pars, int npts, float_sw4 dt);
  ~DRMSource();

  int m_i0, m_j0, m_k0;
  int m_grid;

  float_sw4 getX0() const;
  float_sw4 getY0() const;
  float_sw4 getZ0() const;
  float_sw4 getU( float_sw4 t, float_sw4 ewDt );
  int getComp() const;

  bool myPoint(){ return m_myPoint; }

  private:
  DRMSource();
    
  // void correct_Z_level( EW *a_ew );
  void compute_grid_point( EW *a_ew );
  int spline_interpolation( );

  float_sw4 mX0,mY0,mZ0;
  float_sw4 mDt;
  float_sw4* mPar;
  int mNpts, mComp;
  
  bool m_myPoint;
};
#endif