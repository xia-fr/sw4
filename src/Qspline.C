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
#include "Qspline.h"
#include "F77_FUNC.h"

extern "C" { void F77_FUNC(dgbsv,DGBSV)( int*, int*, int*, int*, double*, int*, int*, double*, int*, int* ); }

using namespace std;

//-----------------------------------------------------------------------
Qspline::Qspline( int npts, float_sw4* fun, float_sw4 tmin, float_sw4 dt, int bclow, int bchigh,
		  float_sw4 s1, float_sw4 t1, float_sw4 sn, float_sw4 tn )
{
   // Quintic spline interpolation of a function defined on a uniform grid.
   // 
   // npts - Number of spline points
   // fun  - The function to interpolate, defined at points k=0,..,npts-1
   // tmin, dt - Gives the assumed uniform grid, t_k = tmin+k*dt, k=0,..,npts-1
   // bclow - Lower end boundary condition, 1-'not a knot', condition (full order one sided)
   //                                       2- Give g', and g'' on lower boundary
   // bchigh - Same as bclow, but for the upper boundary.
   // s1 - g'  at lower boundary, used when bclow=2
   // t1 - g'' at lower boundary, used when bclow=2
   // sn - g'  at upper boundary, used when bchigh=2
   // tn - g'' at upper boundary, used when bchigh=2

   // Need npts >=  7
   if( bclow != 1 && bclow != 2 )
   {
      cout << "ERROR in Qspline constructor, bclow = " << bclow << endl;
      cout << "bclow is reset to 1 "  << endl;
      bclow = 1;
   }
   if( bchigh != 1 && bchigh != 2 )
   {
      cout << "ERROR in Qspline constructor, bchigh = " << bchigh << endl;
      cout << "bchigh is reset to 1 "  << endl;
      bchigh = 1;
   }

   int minpts = 0;
   if( bclow == 1 && bchigh == 1 )
      minpts = 7;
   else if( (bclow == 1 && bchigh == 2) || (bclow==2 && bchigh == 1) )
      minpts = 4;
   else if( bclow == 2 && bchigh == 2 )
      minpts = 2;
   if( npts < minpts )
   {
      cout << "ERROR in Qspline, number of interpolation points, " << npts << ", is too small" << endl;
      cout << " Spline not constructed " << endl;
   }
   m_npts = npts;
   m_polcof = NULL;
// AP: Sept-19-2016: NOTE: tmin is only used within this class
   m_tmin = tmin;
   m_dt = dt;
   m_dti=1/dt;
   // setup and solve linear system for uniform grid Hermite 5th order polynomial splines
   // System is a band matrix with kl lower diagonals, and ku upper diagonals
   int kl = 6, ku = 6;
   if( bclow == 2 )
      ku = 3;
   if( bchigh == 2 )
      kl = 3;
// tmp
   // int myRank;
   // MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   // if (myRank == 0)
   // {
   //   cout << "Qspline constructor:" << endl;
   //   cout << "m_tmin = " << m_tmin << " m_dt= " << m_dt << " m_npts= " << m_npts << endl;
   // }
// end tmp
      
   int lda = (2*kl+ku+1);
   double *mat = new double[lda*2*npts];
   double *rhs = new double[2*npts];
#define A(i,j) mat[kl+ku+1+(i)-(j)-1 + lda*((j)-1)]
   for( int i=0 ; i < lda*2*npts ; i++)
      mat[i] = 0;
   for( int i = 1 ; i <= npts-2 ; i++ )
   {
      A(2*i+1,2*i-1) = -8;
      A(2*i+1,2*i)   = -1;
      A(2*i+1,2*i+2) =  6;
      A(2*i+1,2*i+3) =  8;
      A(2*i+1,2*i+4) = -1;
      rhs[2*i] = 20*(fun[i+1]-2*fun[i]+fun[i-1]);
      A(2*i+2,2*i-1) = -7;
      A(2*i+2,2*i)   = -1;
      A(2*i+2,2*i+1) = -16;
      A(2*i+2,2*i+3) = -7;
      A(2*i+2,2*i+4) =  1;
      rhs[2*i+1] = -15*(fun[i+1]-fun[i-1]);
   }
   // Boundary condition
   if( bclow == 1 )
   {
   // 'not a knot' conditions,
      A(1,1) =  -3;
      A(1,2) =  -0.5;
      A(1,4) =   1;
      A(1,5) =   3;
      A(1,6) =  -0.5;
      rhs[0] =  6*(fun[2]-2*fun[1]+fun[0]);

      A(2,3) =  -3;
      A(2,4) =  -0.5;
      A(2,6) =   1;
      A(2,7) =   3;
      A(2,8) =  -0.5;
      rhs[1] =  6*(fun[3]-2*fun[2]+fun[1]);
   }
   else if( bclow == 2 )
   {
      // Give g' and g'' on boundary
      A(1,1) = 1;
      rhs[0] = s1;
      A(2,2) = 1;
      rhs[1] = t1;
   }
   if( bchigh == 1 )
   {
   // 'not a knot' conditions
      A(2*npts,2*npts)   = -0.5;
      A(2*npts,2*npts-1) =  3;
      A(2*npts,2*npts-2) =  1;
      A(2*npts,2*npts-4) = -0.5;
      A(2*npts,2*npts-5) = -3;
      rhs[2*npts-1] = 6*(fun[npts-1]-2*fun[npts-2]+fun[npts-3]);

      A(2*npts-1,2*npts-2) =  -0.5;
      A(2*npts-1,2*npts-3) =  3;
      A(2*npts-1,2*npts-4) =   1;
      A(2*npts-1,2*npts-6) =  -0.5;
      A(2*npts-1,2*npts-7) = -3;
      rhs[2*npts-2] = 6*(fun[npts-2]-2*fun[npts-3]+fun[npts-4]);

   }
   else if( bchigh == 2 )
   {
      A(2*npts-1,2*npts-1) = 1;
      rhs[2*npts-2] = sn;
      A(2*npts,2*npts) = 1;
      rhs[2*npts-1] = tn;
   }

   int n=2*npts, one=1, info;
   int* ipiv = new int[2*npts];
   F77_FUNC(dgbsv,DGBSV)( &n, &kl, &ku, &one, mat, &lda, ipiv, rhs, &n, &info );
   if( info != 0 )
   {
      cout << "ERROR solving dbgsv in Qspline, info = " << info << endl;
      cout << " Spline not constructed " << endl;
      delete[] mat;
      delete[] ipiv;
      delete[] rhs;
      return;
   }
#undef A
   delete[] mat;
   delete[] ipiv;
   
   m_npts = npts;
   if (npts > 1)
      m_polcof = new float_sw4[6*(npts-1)];
   else
      m_polcof = new float_sw4[6];

   for( int i= 0 ; i < npts-1 ; i++ )
   {
      m_polcof[  6*i] = fun[i];
      m_polcof[1+6*i] = rhs[2*i];
      m_polcof[2+6*i] = rhs[2*i+1]*0.5;
      m_polcof[3+6*i] = 10*(fun[i+1]-fun[i]) - 4*rhs[2*i+2]-6*rhs[2*i]-1.5*rhs[2*i+1]+0.5*rhs[2*i+3];
      m_polcof[4+6*i] =-15*(fun[i+1]-fun[i]) + 7*rhs[2*i+2]+8*rhs[2*i]+1.5*rhs[2*i+1]-rhs[2*i+3];
      m_polcof[5+6*i] =  6*(fun[i+1]-fun[i]) - 3*(rhs[2*i]+rhs[2*i+2])+0.5*(rhs[2*i+3]-rhs[2*i+1]);
   }
   delete[] rhs;
}

//-----------------------------------------------------------------------
// void Qspline::evalf( double t, double& f )
// {
//    int k = static_cast<int>(floor((t-m_tmin)/m_dt));
//    if( k < 0 )
//       k = 0;
//    if( k > m_npts-2 )
//       k = m_npts-2;
//    double arg=(t-(m_tmin+k*m_dt))/m_dt;
//    f = m_polcof[6*k] + m_polcof[1+6*k]*arg + m_polcof[2+6*k]*arg*arg + m_polcof[3+6*k]*arg*arg*arg +
//        m_polcof[4+6*k]*arg*arg*arg*arg + m_polcof[5+6*k]*arg*arg*arg*arg*arg; 
// }

// //-----------------------------------------------------------------------
// void Qspline::evaldd( double t, double& f, double& f1, double& f2, double& f3, double& f4 )
// {
//    int k = static_cast<int>(floor((t-m_tmin)/m_dt));
//    if( k < 0 )
//       k = 0;
//    if( k > m_npts-2 )
//       k = m_npts-2;
//    double arg=(t-(m_tmin+k*m_dt))/m_dt;
//    f = m_polcof[6*k] + m_polcof[1+6*k]*arg + m_polcof[2+6*k]*arg*arg + m_polcof[3+6*k]*arg*arg*arg +
//        m_polcof[4+6*k]*arg*arg*arg*arg + m_polcof[5+6*k]*arg*arg*arg*arg*arg; 
//    f1= (m_polcof[1+6*k] + 2*m_polcof[2+6*k]*arg + 3*m_polcof[3+6*k]*arg*arg + 4*m_polcof[4+6*k]*arg*arg*arg+
// 	5*m_polcof[5+6*k]*arg*arg*arg*arg)*m_dti;
//    f2 = (2*m_polcof[2+6*k] + 6*m_polcof[3+6*k]*arg + 12*m_polcof[4+6*k]*arg*arg + 20*m_polcof[5+6*k]*arg*arg*arg)*m_dti*m_dti;
//    f3 = (6*m_polcof[3+6*k] + 24*m_polcof[4+6*k]*arg + 60*m_polcof[5+6*k]*arg*arg)*m_dti*m_dti*m_dti;
//    f4 = (24*m_polcof[4+6*k] + 120*m_polcof[5+6*k]*arg)*m_dti*m_dti*m_dti*m_dti;
// }

//-----------------------------------------------------------------------
Qspline::~Qspline()
{
   if( m_polcof != NULL )
      delete[] m_polcof;
}
