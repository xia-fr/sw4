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
// # This file is part of SW4, Version: 3.0
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

//#include <chrono>

#define SQR(x) ((x)*(x))

using namespace std;

//---------------------------------------------------------------
// 
// 
//---------------------------------------------------------------
bool EW::isEQLConverged()
{
  // Print convergence information
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if( myRank == 0 && m_iter_EQL > 0)
  {
    string indent  = "\n       ";
    string indents =   "       ";

    cout << endl << "Checking for equivalent linear convergence..." << endl;
      
    cout << indent << "------------- EQL Convergence Statistics ----------------"
   << indent << "Convergence criteria: at most " << m_convPercent_EQL << "%" 
   << indent << "                      change in Vs at any node"
   << indent << endl;
   for (int ix = 0; ix < m_vsConv_EQL.size(); ix++)
   {
    if (ix == 0)
      cout << indents <<                0.0 << " <= Vs <  " << m_vsBins_EQL[0] << " : " << m_vsConv_EQL[0] << " %" << endl;
    else if (ix == m_vsConv_EQL.size()-1)
      cout << indents << m_vsBins_EQL[ix-1] << " <= Vs                   : " << m_vsConv_EQL[ix] << " %" << endl;
    else
      cout << indents << m_vsBins_EQL[ix-1] << " <= Vs <  " << m_vsBins_EQL[ix] << " : " << m_vsConv_EQL[ix] << " %" << endl;
   }
    cout  << indents << "---------------------------------------------------------" << endl << endl;
  }

  // Check for convergence
  m_conv_EQL = true;
  for (int ix = 0; ix < m_vsConv_EQL.size(); ix++)
  {
    if (m_vsConv_EQL[ix] > m_convPercent_EQL)
      m_conv_EQL = false;

    // Reset
    m_vsConv_EQL[ix] = 0.0;
  }

  if (m_iter_EQL >= m_iterLim_EQL)
    m_conv_EQL = true;

  // Reset equivalent linear arrays
  for( int g = 0 ; g < mNumberOfGrids; g++)
  {
    mEmax[g].set_to_minusOne();
    for (int a = 0; a < m_srctype_EQL; a++)
    {
      U_EQL[g][a].set_to_zero();
    }
  }

  // turn off eql strain tracking once we are running
  // the final iteration
  if (m_conv_EQL)
    m_use_EQL = false; 
    
  return m_conv_EQL;
} // end isEQLConverged

//---------------------------------------------------------------
// 
// 
//---------------------------------------------------------------
void EW::calculateEQLUpdate(vector<Source*> & a_Sources)
{
  if (m_myRank == 0 && m_iterLim_EQL != 0)
    cout << endl << "Calculating updated equivalent linear material properties..." << endl;

  // If there is more than one grid, we apply smoothing to the material properties at the
  // boundaries between grids to ensure good material property compatibility at bounds
  if (mNumberOfGrids > 1)
    preprocessInterfaceEQL(mEmax);

  // Iterate over every item in grid
  for( int g = 0 ; g < mNumberOfGrids; g++)
  {
    //vector<float_sw4> maxvsConv(m_vsConv_EQL.size(), 0.0);
    float_sw4 maxvsConv0 = 0.0;
    float_sw4 maxvsConv1 = 0.0;
    float_sw4 maxvsConv2 = 0.0;
    float_sw4 maxvsConv3 = 0.0;
    #pragma omp parallel for reduction(max:maxvsConv0,maxvsConv1,maxvsConv2,maxvsConv3)
    for( int k = m_kStart[g] ; k <= m_kEnd[g]; k++ )
    {
      for( int j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
      {
        for( int i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
        {
          // Save previous EQL mu value
          float_sw4 vs_old = sqrt(mMu[g](i,j,k)/mRho[g](i, j, k));

          // Reset the material arrays
          //mQs[g](i,j,k) = mQsOrig_EQL[g](i,j,k);
          //mMu[g](i,j,k) = mMuOrig_EQL[g](i,j,k);
          mQs[g](i,j,k) = -1.0;
          mQp[g](i,j,k) = -1.0;
          mMu[g](i,j,k) = -1.0;

          // If is a point in the EQL domain
          if ( (i >= m_iStartIntEQL[g] && i <= m_iEndIntEQL[g]) &&
               (j >= m_jStartIntEQL[g] && j <= m_jEndIntEQL[g]) &&
               (k >= m_kStartIntEQL[g] && k <= m_kEndIntEQL[g]))
          {
            float_sw4 ggmax = calculateDarendeli(i, j, k, g, a_Sources);

            // Calculate equivalent linear mu as a product of the calculated
            // ggmax value, and the viscoelastic mu (in mMuOrig_EQL)
            mMu[g](i, j, k) = mMuOrig_EQL[g](i, j, k) * ggmax;

            // Calculate and save convergence statistics based on vs
            float_sw4 vs = sqrt(mMuOrig_EQL[g](i, j, k)/mRho[g](i, j, k));
            float_sw4 vs_new = sqrt(mMu[g](i, j, k)/mRho[g](i, j, k));

            if(m_iterLim_EQL!=0) 
            {
              if (vs < m_vsBins_EQL[0])
                maxvsConv0 = max(maxvsConv0, abs(vs_new-vs_old)/vs_old*100);
              else if (vs < m_vsBins_EQL[1])
                maxvsConv1 = max(maxvsConv1, abs(vs_new-vs_old)/vs_old*100);
              else if (vs < m_vsBins_EQL[2])
                maxvsConv2 = max(maxvsConv2, abs(vs_new-vs_old)/vs_old*100);
              else
                maxvsConv3 = max(maxvsConv3, abs(vs_new-vs_old)/vs_old*100);
            }

            // Debug
            //float_sw4 x, y, z;
            //getCoordinates(i, j, k, g, x, y, z);
            //if (x == 5000.0 && y == 5000.0 && z == 2500.0)
            //{
            //  cout << "i" << i << " j" << j << " k" << k << endl;
            //  cout << "x" << x << " y" << y << " z" << z << endl;
            //  cout << "grid " <<  g << endl;
            //  cout << "eref " <<  mEmax[g](i, j, k) << endl;
            //  cout << "vs " << vs_new << endl;
            //  cout << "ggmax " << ggmax << endl;
            //  cout << " " << endl;
            //}
            
          } // End if point in EQL domain
          
        } // End i
      } // End j
    } // End k

    // Communicate max values for convergence calculation
    if(m_iterLim_EQL!=0)
    {
      MPI_Allreduce(&maxvsConv0,&m_vsConv_EQL[0],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
      MPI_Allreduce(&maxvsConv1,&m_vsConv_EQL[1],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
      MPI_Allreduce(&maxvsConv2,&m_vsConv_EQL[2],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
      MPI_Allreduce(&maxvsConv3,&m_vsConv_EQL[3],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
    }

    // First: communicate ghost points between processors in x and y dir.
    communicate_array( mQs[g], g );
    communicate_array( mQp[g], g );
    communicate_array( mMu[g], g );

    // Extrapolate to define material properties above the free surface (topography)
    if (g == mNumberOfGrids-1)
    {
      extrapolateInZ(g, mMu[g], true, false );
      extrapolateInZ(g, mQs[g], true, false );
      extrapolateInZ(g, mQp[g], true, false );
    }
    // OR
    // Extrapolate to ghost points at bottom of domain
    // (every grid does this one for continuity)
    extrapolateInZ(g, mMu[g], false, true );
    extrapolateInZ(g, mQs[g], false, true );
    extrapolateInZ(g, mQp[g], false, true );

    // // Debug: print out the ghost points
    // if (proc_zero())
    // {
    //   cout << "1(inner): " << mMu[g](1, m_iEnd[g], m_jStartIntEQL[g], m_kStartIntEQL[g]) << endl;
    //   cout << "2(inner): " << mMu[g](1, m_iEnd[g]-1, m_jStartIntEQL[g], m_kStartIntEQL[g]) << endl;
    //   cout << "1(inner2): " << mMu[g](1, m_iEnd[g], m_jStartIntEQL[g], m_kStartIntEQL[g]+1) << endl;
    //   cout << "2(inner2): " << mMu[g](1, m_iEnd[g]-1, m_jStartIntEQL[g], m_kStartIntEQL[g]+1) << endl;
    //   cout << "1(above?): " << mMu[g](1, m_iStartIntEQL[g], m_jStartIntEQL[g], m_kStart[g]) << endl;
    //   cout << "2(above?): " << mMu[g](1, m_iEndIntEQL[g], m_jStartIntEQL[g], m_kStart[g]) << endl;
    //   cout << "1(corner?): " << mMu[g](1, m_iEnd[g], m_jStartIntEQL[g], m_kStart[g]) << endl;
    //   cout << "2(corner?): " << mMu[g](1, m_iEnd[g]-1, m_jStartIntEQL[g], m_kStart[g]) << endl;
    // }

  } // End g

  // If there is more than one grid, we apply smoothing to the material properties at the
  // boundaries between grids to ensure good material property compatibility at bounds
  if (mNumberOfGrids > 1)
  {
    processInterfaceEQL(mMu);
    processInterfaceEQL(mQs);
    processInterfaceEQL(mQp);
  }

  // Extrapolate to ghost points in x and y, if they were not set by the previous routines.
  // e.g. ghost points on the outside boundary of the model and corner points
  // after this point there should be no "-1" values left in any grid
  extrapolateInXY_EQL(mMu);
  extrapolateInXY_EQL(mQs);
  extrapolateInXY_EQL(mQp);

  // Check materials
  if (m_iterLim_EQL != 0)
    check_materials_EQL();

} // end calculateEQLUpdate

//---------------------------------------------------------------
// 
// 
//---------------------------------------------------------------
float_sw4 EW::calculateDarendeli(int i, int j, int k, int g, vector<Source*> & a_Sources)
{
  // Pressures using standard gravity and depth of node
  // - makes an assumption that the vertical pressure on the node
  //   is from a column made up of the same material as the node itself
  //   - FUTURE WORK: more precise calculation of pressure using exact densities
  //                  of nodes in the column
  // - assumed to be 1 atm at depth=0
  // - Units of density assumed to be kg/m^3

  float_sw4 x, y, z;
  getCoordinates(i, j, k, g, x, y, z);
  float_sw4 ptDepth;
  getDepth(x, y, z, ptDepth); // Gets total depth below free surface including topography

  float_sw4 sigma_v0 = mRho[g](i, j, k)*9.80665*ptDepth + 101325; // vertical overburden pressure (units: Pa)

  // Plasticity index as a fn. of shear wave velocity
  float_sw4 vs = sqrt(mMuOrig_EQL[g](i, j, k)/mRho[g](i, j, k));
  float_sw4 PL;
  if (vs <= 200)
    PL = 10.0;
  else if (vs <= 360)
    PL = 5.0;
  else
    PL = 0.0;  

  // Overconsolidation ratio
  float_sw4 sigma_p0 = 0.106*pow(vs, 1.47) * 1000; // Mayne, Robertson, Lunne (1998) "Clay stress history evaluated from seismic piezocone tests"  # noqa: E501,E226
  float_sw4 OCR = sigma_p0 / sigma_v0; // Units: Pa

  // Lateral earth pressure coefficient at rest
  // from OCR using the empirical formula by Mayne & Kulhawy (1982)
  // - Makes assumption that internal effective friction angle is
  //   30 deg.
  float_sw4 k0 = (1 - sin(M_PI/6.0)) * pow(OCR, sin(M_PI/6.0));

  // mean confining pressure (units: atm)
  float_sw4 sigma_0 = (2.0 * k0 + 1.0) / 3.0 * sigma_v0 / 101325; 

  // values from darendeli
  float_sw4 nr_cycle = 10.0;
  float_sw4 frq = 1.0;
  float_sw4 phi1 = 0.0352;
  float_sw4 phi2 = 0.0010;
  float_sw4 phi3 = 0.3246;
  float_sw4 phi4 = 0.3483;
  float_sw4 phi5 = 0.9190;
  float_sw4 phi6 = 0.8005;
  float_sw4 phi7 = 0.0129;
  float_sw4 phi8 = -0.1069;
  float_sw4 phi9 = -0.2889;
  float_sw4 phi10 = 0.2919;
  float_sw4 phi11 = 0.6329;
  float_sw4 phi12 = -0.0057;
  float_sw4 a = phi5;

  float_sw4 c1 = -1.1143*(a*a) + 1.8618*a + 0.2523; // Darendeli (2001), page 226
  float_sw4 c2 = 0.0805*(a*a) - 0.0710*a - 0.0095;
  float_sw4 c3 = -0.0005*(a*a) + 0.0002*a + 0.0003;
  float_sw4 b = phi11 + phi12*log(nr_cycle);

  float_sw4 gamma_r = (phi1 + phi2*PL*(pow(OCR, phi3))) * pow(sigma_0, phi4);
  
  // Grab relevant strain value
  // Calculate cap to Eref based on user defined global min. vs
  float_sw4 muMin = mRho[g](i,j,k)*m_vsMin*m_vsMin;
  float_sw4 ggmax_min = muMin/mMuOrig_EQL[g](i,j,k);
  //float_sw4 maxEref = gamma_r * pow((1/ggmax_min)-1, 1/a);
  float_sw4 Eref = mEmax[g](i, j, k) * 0.65 * 100; // strain in percent

  float_sw4 ggmax = max(1.0 / (1 + pow(Eref/gamma_r, a)), ggmax_min);
  ggmax = min(1.0, ggmax);

  // Attenuation 
  float_sw4 xi = 1.0 / (2.0 * (mQsOrig_EQL[g](i, j, k))) * 100; // units: %
  if (Eref > (gamma_r/100000.0)) // threshold for D_masing_1 calculation
  {
    float_sw4 D_masing_1 = ( (100.0 / M_PI) 
                          * ( 
                                4.0
                                * (Eref-gamma_r * log((Eref+gamma_r)/gamma_r)) 
                                / (Eref*Eref / (Eref+gamma_r)) - 2.0 
                             ) 
                         );

    // D_masing_1 is the masing calculation for damping at a=1.0, makes no physical sense
    // if it is negative (numerical artifact from equation limitations for small eref)
    if ( (D_masing_1 > 0.0) )
    {
      float_sw4 D_masing = c1*D_masing_1 
                         + c2*D_masing_1*D_masing_1 
                         + c3*D_masing_1*D_masing_1*D_masing_1; // units: %

      float_sw4 D_min = (phi6 + phi7 * M_PI * pow(OCR, phi8)) * pow(sigma_0, phi9) * (1 + phi10 * log(frq));
      float_sw4 D_adj = b * pow(ggmax, 0.1) * D_masing; // units: %
      xi = xi+D_adj;
    } 
  }

  // Apply attenuation value using (Q = 1/(2*xi)) and (Qp = 2*Qs)
  mQs[g](i, j, k) = 1.0 / (2.0 * (xi/100.0));
  mQp[g](i, j, k) = 2.0 * mQs[g](i, j, k);

  return ggmax;

} // end isEQLConverged

//---------------------------------------------------------------
// Given a grid and the location in ijk, return the xyz coordinates.
//---------------------------------------------------------------
void EW::getCoordinates(int i, int j, int k, int g,
                             float_sw4 &x, float_sw4 &y, float_sw4 &z)
{
  if (g < mNumberOfCartesianGrids) // Cartesian grids
  {
    x = (i-1)*mGridSize[g];
    y = (j-1)*mGridSize[g];
    z = m_zmin[g]+(k-1)*mGridSize[g];
  }
  else // is a Curvilinear grid
  {
    x = mX[g](i,j,k);
    y = mY[g](i,j,k);
    z = mZ[g](i,j,k);
  }
}

//---------------------------------------------------------------
// Accounts for the fact that sometimes the input in SW4 is
// velocity or acceleration. This is calculated for some points
// outside the EQL domain, since the strain calculation for EQL 
// may rely on displacements that are one node outside the 
// boundary of the EQL domain.
//
// Displacement is calculated from velocity using forward euler
// - u(t + dt) = u(t) + v(t)*dt
//
// Displacement is calculated from acceleration using a second
// order Taylor Series expansion
// - u(t + dt) = u(t) + v(t)*dt + 0.5*a(t)*dt^2
//
// If U_EQL is of size [gx1 x (EQL domain size)], then
// U_EQL[g][0] is displacement.
//
// If U_EQL is of size [gx2 x (EQL domain size)], then
// U_EQL[g][0] is displacement and U_EQL[g][1] is velocity.
//---------------------------------------------------------------
void EW::calcEQLDispl(vector<Sarray> &U)
{
  // Iterate over every item in grid
  for( int g = 0 ; g < mNumberOfGrids; g++)
  {
    #pragma omp parallel for     
    for( int k = m_kStart[g] ; k <= m_kEnd[g]; k++ )
    {
      for( int j = m_jStart[g] ; j <= m_jEnd[g]; j++ )
      {
        for( int i = m_iStart[g] ; i <= m_iEnd[g] ; i++ )
        {
          // Debugging
          /*
          if ((i == 75) && (j == 75) && (k == 75))
          {
            cout << "i" << i << " j" << j << " k" << k << endl;
            cout << "Uex  " << U_EQL[g][0](1, i,j,k) << endl;
            cout << "Uey  " << U_EQL[g][0](2, i,j,k) << endl;
            cout << "Uez  " << U_EQL[g][0](3, i,j,k) << endl;
            cout << "Ux  " << U[g](1, i,j,k) << endl;
            cout << "Uy  " << U[g](2, i,j,k) << endl;
            cout << "Uz  " << U[g](3, i,j,k) << endl;
          }
          */
          
          if (m_srctype_EQL == 1) // Source (Up) is velocity
          {
            // U            := Velocity at current time step
            // U_EQL[g][0]  := Displacement at current time step

            // Calculate displacement at new time step
            // u(t + dt) = u(t) + v(t)*dt (forward euler)
            U_EQL[g][0](1, i,j,k) += U[g](1, i,j,k)*mDt;
            U_EQL[g][0](2, i,j,k) += U[g](2, i,j,k)*mDt;
            U_EQL[g][0](3, i,j,k) += U[g](3, i,j,k)*mDt;
            
            // Updated:
            // U_EQL[g][0]  := Displacement at new time step
            // U_EQL[g][1]  := Velocity at new time step
          }
          else // m_srctype_EQL == 2,  Source (Up) is acceleration
          {
            // U            := Acceleration at current time step
            // U_EQL[g][0]  := Displacement at current time step
            // U_EQL[g][1]  := Velocity at current time step

            // Calculate displacement at new time step using current values
            // u(t + dt) = u(t) + v(t)*dt + 0.5*a(t)*dt^2 (Taylor expansion)
            U_EQL[g][0](1, i,j,k) += U_EQL[g][1](1, i,j,k)*mDt
                                          + 0.5*U[g](1, i,j,k)*mDt*mDt;
            U_EQL[g][0](2, i,j,k) += U_EQL[g][1](2, i,j,k)*mDt
                                          + 0.5*U[g](2, i,j,k)*mDt*mDt;
            U_EQL[g][0](3, i,j,k) += U_EQL[g][1](3, i,j,k)*mDt
                                          + 0.5*U[g](3, i,j,k)*mDt*mDt;

            // Calculate velocity at new time step
            U_EQL[g][1](1, i,j,k) += U[g](1, i,j,k)*mDt;
            U_EQL[g][1](2, i,j,k) += U[g](2, i,j,k)*mDt;
            U_EQL[g][1](3, i,j,k) += U[g](3, i,j,k)*mDt;

            // Updated:
            // U_EQL[g][0]  := Displacement at new time step
            // U_EQL[g][1]  := Velocity at new time step
          }
        }
      }
    }
    communicate_array( U_EQL[g][0], g );
    if (m_srctype_EQL == 2) communicate_array( U_EQL[g][1], g );
  }

} // end calcEQLDispl

//---------------------------------------------------------------
// 
// 
//---------------------------------------------------------------
bool EW::updateEmax(vector<Sarray*> &U)
{
  int supergrid_pts;
  bool updated = false;
  for( int g = 0 ; g < mNumberOfGrids; g++)
  {
    supergrid_pts = static_cast<int>(m_supergrid_width/mGridSize[g]);
    #pragma omp parallel for     
    for( int k = m_kStartIntEQL[g] ; k <= m_kEndIntEQL[g]; k++ )
    {
      for( int j = m_jStartIntEQL[g] ; j <= m_jEndIntEQL[g]; j++ )
      {
        for( int i = m_iStartIntEQL[g] ; i <= m_iEndIntEQL[g] ; i++ )
        {
          float_sw4 x, y, z;
          float_sw4 ptDepth;
          float_sw4 vs = sqrt(mMuOrig_EQL[g](i, j, k)/mRho[g](i, j, k));
          getCoordinates(i, j, k, g, x, y, z);
          getDepth(x, y, z, ptDepth); // Gets total depth below free surface including topography

          // if the shear wave velocity is greater than vslim,
          // or if too close to the source,
          // or if depth below the set limit,
          // or if its a supergrid point
          // then don't apply eql (aka, set max strain = 0)
          if ( (vs > m_vslim_eql) ||
               (ptDepth > m_max_depth_eql) ||
               (z > m_max_z_eql) ||
               (m_min_dist_to_srcs[g](i,j,k) < m_src_Dmin) ||
               (i <= supergrid_pts) ||
               (i >= m_global_nx[g]-supergrid_pts) ||
               (j <= supergrid_pts) ||
               (j >= m_global_ny[g]-supergrid_pts) ||
               (g == 0 && k >= m_global_nz[g]-supergrid_pts)
             )
          { 
            mEmax[g](i, j, k) = 0.0;
          }
          else
          { 
            if (g < mNumberOfCartesianGrids) // must be a Cartesian grid
            {
              //       int i=m_i0, j=m_j0, k=m_k0, g=m_grid0;
              float_sw4 factor = 1.0 / (2 * mGridSize[g]);
              float_sw4 duydx = (U[g][0](2, i + 1, j, k) - U[g][0](2, i - 1, j, k)) * factor;
              float_sw4 duzdx = (U[g][0](3, i + 1, j, k) - U[g][0](3, i - 1, j, k)) * factor;
              float_sw4 duxdy = (U[g][0](1, i, j + 1, k) - U[g][0](1, i, j - 1, k)) * factor;
              float_sw4 duzdy = (U[g][0](3, i, j + 1, k) - U[g][0](3, i, j - 1, k)) * factor;
              float_sw4 duxdz = (U[g][0](1, i, j, k + 1) - U[g][0](1, i, j, k - 1)) * factor;
              float_sw4 duydz = (U[g][0](2, i, j, k + 1) - U[g][0](2, i, j, k - 1)) * factor;
              float_sw4 duxdx = (U[g][0](1, i + 1, j, k) - U[g][0](1, i - 1, j, k)) * factor;
              float_sw4 duydy = (U[g][0](2, i, j + 1, k) - U[g][0](2, i, j - 1, k)) * factor;
              float_sw4 duzdz = (U[g][0](3, i, j, k + 1) - U[g][0](3, i, j, k - 1)) * factor;

              // strain tensor values
              float_sw4 E11 = (duxdx);
              float_sw4 E22 = (duydy);
              float_sw4 E33 = (duzdz);
              float_sw4 E12 = (0.5 * (duydx + duxdy));
              float_sw4 E13 = (0.5 * (duzdx + duxdz));
              float_sw4 E23 = (0.5 * (duydz + duzdy));

              // strain tensor invariants
              float_sw4 I1 = (E11 + E22 + E33);
              float_sw4 I2 = (E11 * E22 + E22 * E33 + E33 * E11 - E12 * E12 - E13 * E13 - E23 * E23);
              float_sw4 I3 = (E11 * E22 * E33 - E11 * E23 * E23 - E22 * E13 * E13 - E33 * E12 * E12 + 2.0 * E12 * E13 * E23);

              // second invariant of strain deviator tensor
              float_sw4 J2 = (1.0/3.0) * I1*I1 - I2;

              // approximation of max shear strain using J2
              float_sw4 J2_est = sqrt(4*J2);

              // if (m_iterLim_EQL != 0)
              //     mEmax[g](i, j, k) = max(J2_est, mEmax[g](i, j, k));
              // else
              //     mEmax[g](i, j, k) = J2_est;
              if (m_iterLim_EQL != 0)
              {
                if (J2_est > mEmax[g](i, j, k))
                {
                  mEmax[g](i, j, k) = J2_est;
                  updated = true;
                }
              }
              // else
              //   mEmax[g](i, j, k) = J2_est;
            }
            else // must be curvilinear
            {
              //       int i=m_i0, j=m_j0, k=m_k00, g=m_grid0;
              float_sw4 factor = 0.5 / sqrt(mJ[g](i, j, k));
              //      float_sw4 duxdq = (U[g](1,i+1,j,k) - U[g](1,i-1,j,k));
              float_sw4 duydq = (U[g][0](2, i + 1, j, k) - U[g][0](2, i - 1, j, k));
              float_sw4 duzdq = (U[g][0](3, i + 1, j, k) - U[g][0](3, i - 1, j, k));
              float_sw4 duxdr = (U[g][0](1, i, j + 1, k) - U[g][0](1, i, j - 1, k));
              //      float_sw4 duydr = (U[g](2,i,j+1,k) - U[g](2,i,j-1,k));
              float_sw4 duzdr = (U[g][0](3, i, j + 1, k) - U[g][0](3, i, j - 1, k));
              float_sw4 duxds = (U[g][0](1, i, j, k + 1) - U[g][0](1, i, j, k - 1));
              float_sw4 duyds = (U[g][0](2, i, j, k + 1) - U[g][0](2, i, j, k - 1));
              float_sw4 duzds = (U[g][0](3, i, j, k + 1) - U[g][0](3, i, j, k - 1));
              float_sw4 duzdy = (mMetric[g](1, i, j, k) * duzdr + mMetric[g](3, i, j, k) * duzds) * factor;
              float_sw4 duydz = (mMetric[g](4, i, j, k) * duyds) * factor;
              float_sw4 duxdz = (mMetric[g](4, i, j, k) * duxds) * factor;
              float_sw4 duzdx = (mMetric[g](1, i, j, k) * duzdq + mMetric[g](2, i, j, k) * duzds) * factor;
              float_sw4 duydx = (mMetric[g](1, i, j, k) * duydq + mMetric[g](2, i, j, k) * duyds) * factor;
              float_sw4 duxdy = (mMetric[g](1, i, j, k) * duxdr + mMetric[g](3, i, j, k) * duxds) * factor;
              float_sw4 duxdx = (mMetric[g](1, i, j, k) * (U[g][0](1, i + 1, j, k) - U[g][0](1, i - 1, j, k)) +
                                  mMetric[g](2, i, j, k) * (U[g][0](1, i, j, k + 1) - U[g][0](1, i, j, k - 1))) *
                              factor;
              float_sw4 duydy = (mMetric[g](1, i, j, k) * (U[g][0](2, i, j + 1, k) - U[g][0](2, i, j - 1, k)) +
                                  mMetric[g](3, i, j, k) * (U[g][0](2, i, j, k + 1) - U[g][0](2, i, j, k - 1))) *
                              factor;
              float_sw4 duzdz = (mMetric[g](4, i, j, k) * (U[g][0](3, i, j, k + 1) - U[g][0](3, i, j, k - 1))) * factor;

              // strain tensor values
              float_sw4 E11 = (duxdx);
              float_sw4 E22 = (duydy);
              float_sw4 E33 = (duzdz);
              float_sw4 E12 = (0.5 * (duydx + duxdy));
              float_sw4 E13 = (0.5 * (duzdx + duxdz));
              float_sw4 E23 = (0.5 * (duydz + duzdy));

              // invariants
              float_sw4 I1 = (E11 + E22 + E33);
              float_sw4 I2 = (E11 * E22 + E22 * E33 + E33 * E11 - E12 * E12 - E13 * E13 - E23 * E23);
              float_sw4 I3 = (E11 * E22 * E33 - E11 * E23 * E23 - E22 * E13 * E13 - E33 * E12 * E12 + 2.0 * E12 * E13 * E23);

              // second invariant of strain deviator tensor
              float_sw4 J2 = (1.0/3.0) * I1*I1 - I2;

              // approximation of max shear strain using J2
              float_sw4 J2_est = sqrt(4*J2);

              if (m_iterLim_EQL != 0)
              {
                if (J2_est > mEmax[g](i, j, k))
                {
                  mEmax[g](i, j, k) = J2_est;
                  updated = true;
                }
              }
              // else
              //   mEmax[g](i, j, k) = J2_est;
            }
          }
        }
      }
    }
    //communicate_array( mEmax[g], g );
  } // end for all grids
  return updated;
} // end updateEmax

//-----------------------------------------------------------------------
// Preprocesses the strain field to smooth out values near
// interfaces between two grids
//-----------------------------------------------------------------------
void EW::preprocessInterfaceEQL(vector<Sarray>& field)
{
  // Each grid interacts with the grid at the interface above it, 
  // the top-most grid does not need to do anything
  for( int g = 0 ; g < mNumberOfGrids-1; g++)
  {
    // Get the grid spacings
    float_sw4 thisGrid_h = mGridSize[g];
    float_sw4 aboveGrid_h = mGridSize[g+1];

    // First, interpolate vertically across the shared boundary using two points
    // north and south of the boundary point from each grid.
    // We set the values to be in the finer grid, since they'll overwrite the
    // coarser grid later on at these overlap areas
    int ni=m_iEndIntEQL[g]-m_iStartIntEQL[g]+1, nj=m_jEndIntEQL[g]-m_jStartIntEQL[g]+1;
    int g0_ib=m_iStartInt[g], g0_jb=m_jStartInt[g], g1_ib=m_iStartInt[g+1], g1_jb=m_jStartInt[g+1];
    int g0_k=m_kStartInt[g], g1_k=m_kEndInt[g+1];
    for( int jx = 0; jx <= nj; jx++ )
    for( int ix = 0; ix <= ni; ix++ )
    {
      // Aka if this is interface between a cartesian and curvilinear grid
      if (thisGrid_h == aboveGrid_h)
      {
        field[g+1](g1_ib+ix, g1_jb+jx, g1_k) = (field[g](g0_ib+ix, g0_jb+jx, g0_k+1)
                                               +field[g](g0_ib+ix, g0_jb+jx, g0_k+2)
                                               +field[g+1](g1_ib+ix, g1_jb+jx, g1_k-1)
                                               +field[g+1](g1_ib+ix, g1_jb+jx, g1_k-2))/4.0;
      }
      else // aboveGrid_h = 1/2 * thisGrid_h
      {
        field[g+1](g1_ib+(2*ix), g1_jb+(2*jx), g1_k) = 
                            (field[g](g0_ib+ix, g0_jb+jx, g0_k+1)
                            +field[g](g0_ib+ix, g0_jb+jx, g0_k+2)
                            +field[g+1](g1_ib+(2*ix), g1_jb+(2*jx), g1_k-1)
                            +field[g+1](g1_ib+(2*ix), g1_jb+(2*jx), g1_k-2))/4.0;
      }
    }
    communicate_array(field[g+1], g+1);

    // Next, if the grid spacings aren't the same, interpolate the floating values
    // on the finer grid using the new values
    // Aka if this is interface between two grids with different discretization width
    if (thisGrid_h != aboveGrid_h)
    {
      interpolateEQL(field, g, g1_k);
    }
  } // end for g
  communicate_arrays(field);
}

//-----------------------------------------------------------------------
// NOTE: this gets passed the index for grid g, but it acts on grid *g+1*
// - Assumes that the points coinciding with the coarse grid have already
//   been set by previous routines
//-----------------------------------------------------------------------
void EW::interpolateEQL(vector<Sarray>& field, int g, int g_k)
{
  int ni=m_iEndInt[g]-m_iStartInt[g]+1, nj=m_jEndInt[g]-m_jStartInt[g]+1;
  int g1_ib=m_iStartInt[g+1], g1_jb=m_jStartInt[g+1];

  // 1. The center points
  for( int jx = 1; jx <= nj; jx++ )
  for( int ix = 1; ix <= ni; ix++ )
  {
    field[g+1](g1_ib+(2*ix-1), g1_jb+(2*jx-1), g_k) = 
      (field[g+1](g1_ib+(2*ix-2), g1_jb+(2*jx-2), g_k)
      +field[g+1](g1_ib+(2*ix-2), g1_jb+(2*jx-0), g_k)
      +field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-2), g_k)
      +field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-0), g_k)
      )/4.0;
  }

  // 2. odd points in ix
  for( int jx = 0; jx <= nj; jx++ )
  for( int ix = 1; ix <= ni; ix++ )
  {
    if (jx == 0)
    {
      field[g+1](g1_ib+(2*ix-1), g1_jb+(2*jx-0), g_k) = 
        (field[g+1](g1_ib+(2*ix-2), g1_jb+(2*jx-0), g_k)
        +field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-0), g_k)
        +field[g+1](g1_ib+(2*ix-1), g1_jb+(2*jx+1), g_k)
        )/3.0;
    }
    else if (jx == nj)
    {
      field[g+1](g1_ib+(2*ix-1), g1_jb+(2*jx-0), g_k) = 
        (field[g+1](g1_ib+(2*ix-2), g1_jb+(2*jx-0), g_k)
        +field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-0), g_k)
        +field[g+1](g1_ib+(2*ix-1), g1_jb+(2*jx-1), g_k)
        )/3.0;
    }
    else
    {
      field[g+1](g1_ib+(2*ix-1), g1_jb+(2*jx-0), g_k) = 
        (field[g+1](g1_ib+(2*ix-2), g1_jb+(2*jx-0), g_k)
        +field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-0), g_k)
        +field[g+1](g1_ib+(2*ix-1), g1_jb+(2*jx+1), g_k)
        +field[g+1](g1_ib+(2*ix-1), g1_jb+(2*jx-1), g_k)
        )/4.0;
    }
  }

  // 3. odd points in jx
  for( int jx = 1; jx <= nj; jx++ )
  for( int ix = 0; ix <= ni; ix++ )
  {
    if (ix == 0)
    {
      field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-1), g_k) = 
        (field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-2), g_k)
        +field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-0), g_k)
        +field[g+1](g1_ib+(2*ix+1), g1_jb+(2*jx-1), g_k)
        )/3.0;
    }
    else if (ix == ni)
    {
      field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-1), g_k) = 
        (field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-2), g_k)
        +field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-0), g_k)
        +field[g+1](g1_ib+(2*ix-1), g1_jb+(2*jx-1), g_k)
        )/3.0;
    }
    else
    {
      field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-1), g_k) = 
        (field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-2), g_k)
        +field[g+1](g1_ib+(2*ix-0), g1_jb+(2*jx-0), g_k)
        +field[g+1](g1_ib+(2*ix+1), g1_jb+(2*jx-1), g_k)
        +field[g+1](g1_ib+(2*ix-1), g1_jb+(2*jx-1), g_k)
        )/4.0;
    }
  }
}

//-----------------------------------------------------------------------
void EW::processInterfaceEQL(vector<Sarray>& field)
{
  // Each grid interacts with the grid at the interface above it, 
  // the top-most grid does not need to do anything
  for( int g = 0 ; g < mNumberOfGrids-1; g++)
  {
    // Get the grid spacings
    float_sw4 thisGrid_h = mGridSize[g];
    float_sw4 aboveGrid_h = mGridSize[g+1];

    // First, set the interface points and ghost points on grid g
    // to be the values on grid g+1 (which is either equal in
    // discretization or finer, so it overrides)
    int ni=m_iEndInt[g]-m_iStartInt[g]+1, nj=m_jEndInt[g]-m_jStartInt[g]+1;
    int g0_ib=m_iStartInt[g], g0_jb=m_jStartInt[g], g1_ib=m_iStartInt[g+1], g1_jb=m_jStartInt[g+1];
    int g0_k=m_kStartInt[g], g1_k=m_kEndInt[g+1];
    for( int kx = 0; kx <= m_ghost_points; kx++)
    for( int jx = 0; jx <= nj; jx++ )
    for( int ix = 0; ix <= ni; ix++ )
    {
      // Aka if this is interface between a cartesian and curvilinear grid
      if (thisGrid_h == aboveGrid_h)
      {
        field[g](g0_ib+ix, g0_jb+jx, g0_k-kx) = field[g+1](g1_ib+ix, g1_jb+jx, g1_k-kx);
      }
      else // aboveGrid_h = 1/2 * thisGrid_h
      {
        field[g](g0_ib+ix, g0_jb+jx, g0_k-kx) = field[g+1](g1_ib+(2*ix), g1_jb+(2*jx), g1_k-(2*kx));
      }
    }

    // Second, vertically interpolate the points on grid g near the interface
    // to smooth out the transition just a little bit
    vector<float_sw4> tmps;
    tmps.resize(6);
    for( int jx = 0; jx <= nj; jx++ )
    for( int ix = 0; ix <= ni; ix++ )
    {
      // Calculate all the averaged values first
      for( int kx = 1; kx <= 6; kx++)
      {
        // There will always be min. two ghost points
        tmps[kx-1] = (field[g](g0_ib+ix, g0_jb+jx, g0_k+kx-3)
                      +field[g](g0_ib+ix, g0_jb+jx, g0_k+kx-2)
                      +field[g](g0_ib+ix, g0_jb+jx, g0_k+kx-1)
                      +field[g](g0_ib+ix, g0_jb+jx, g0_k+kx)
                      +field[g](g0_ib+ix, g0_jb+jx, g0_k+kx+1)
                      +field[g](g0_ib+ix, g0_jb+jx, g0_k+kx+2)
                      +field[g](g0_ib+ix, g0_jb+jx, g0_k+kx+3)
                     )/7.0;
      }
      // Then set them
      for( int kx = 1; kx <= 6; kx++)
        field[g](g0_ib+ix, g0_jb+jx, g0_k+kx) = tmps[kx-1];
    }

    // // Third, set the finer mesh ghost points equal to the coarser mesh
    // // points where they overlap (in layer 2)
    // for( int jx = 0; jx < nj; jx++ )
    // for( int ix = 0; ix < ni; ix++ )
    // {
    //   // If this is interface between a cartesian and curvilinear grid
    //   // and the grid spacing is the 'same'
    //   if (thisGrid_h == aboveGrid_h)
    //   {
    //     field[g+1](g1_ib+ix, g1_jb+jx, g1_k+1) = field[g](g0_ib+ix, g0_jb+jx, g0_k+1);
    //     field[g+1](g1_ib+ix, g1_jb+jx, g1_k+2) = field[g](g0_ib+ix, g0_jb+jx, g0_k+2);
    //   }
    //   else // aboveGrid_h = 1/2 * thisGrid_h
    //   {
    //     field[g+1](g1_ib+(2*ix), g1_jb+(2*jx), g1_k+2) = field[g](g0_ib+ix, g0_jb+jx, g0_k+1);
    //   }
    // }

    // // Only need to do this if the two grid discretizations are different
    // // aka, aboveGrid_h = 1/2 * thisGrid_h
    // if (thisGrid_h != aboveGrid_h)
    // {
    //   // Fourth, interpolate between fine mesh ghost points in layer 2 of the ghost points.
    //   interpolateEQL(field, g, g1_k+2);

    //   // Fifth, vertically interpolate the first layer of ghost points
    //   ni=m_iEndIntEQL[g+1]-m_iStartIntEQL[g+1]+1; // nj and ni are now for the finer mesh
    //   nj=m_jEndIntEQL[g+1]-m_jStartIntEQL[g+1]+1;
    //   for( int jx = 0; jx < nj; jx++ )
    //   for( int ix = 0; ix < ni; ix++ )
    //   {
    //       field[g+1](g1_ib+ix, g1_jb+jx, g1_k+1) 
    //         = (field[g+1](g1_ib+ix, g1_jb+jx, g1_k+0) 
    //           +field[g+1](g1_ib+ix, g1_jb+jx, g1_k+2) ) / 2.0;
    //   }
    // }

    // // Lastly, vertically interpolate the points on grid g+1 near the interface
    // // to smooth out the transition just a little bit
    // for( int jx = 0; jx < nj; jx++ ) // nj and ni are now for the finer mesh
    // for( int ix = 0; ix < ni; ix++ )
    // {
    //   // Calculate all the averaged values first
    //   for( int kx = 1; kx <= 6; kx++)
    //   {
    //     // There will always be min. two ghost points
    //     tmps[kx-1] = (field[g+1](g1_ib+ix, g1_jb+jx, g1_k-kx+3)
    //                   +field[g+1](g1_ib+ix, g1_jb+jx, g1_k-kx+2)
    //                   +field[g+1](g1_ib+ix, g1_jb+jx, g1_k-kx+1)
    //                   +field[g+1](g1_ib+ix, g1_jb+jx, g1_k-kx)
    //                   +field[g+1](g1_ib+ix, g1_jb+jx, g1_k-kx-1)
    //                   +field[g+1](g1_ib+ix, g1_jb+jx, g1_k-kx-2)
    //                   +field[g+1](g1_ib+ix, g1_jb+jx, g1_k-kx-3)
    //                  )/7.0;
    //   }
    //   // Then set them
    //   for( int kx = 1; kx <= 6; kx++)
    //     field[g+1](g1_ib+ix, g1_jb+jx, g1_k-kx) = tmps[kx-1];
    // }
  }
}

//-----------------------------------------------------------------------
void EW::extrapolateInXY_EQL( vector<Sarray>& field )
{
   for( int g= 0; g < mNumberOfGrids ; g++ )
   {
    // Extrapolate in i
    #pragma omp parallel for
    for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
    for( int i=m_iStart[g] ; i < m_iStartIntEQL[g] ; i++ )
    {
      if( field[g](i,j,k) == -1 )
        field[g](i,j,k) = field[g](m_iStartIntEQL[g],j,k);
    }
    #pragma omp parallel for
    for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
    for( int j=m_jStart[g] ; j <= m_jEnd[g] ; j++ )
    for( int i=m_iEndInt[g]+1 ; i <= m_iEnd[g] ; i++ )
    {
      if( field[g](i,j,k) == -1 )
        field[g](i,j,k) = field[g](m_iEndInt[g],j,k);
    }

    //Extrapolate in j
    #pragma omp parallel for
    for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
    for( int j=m_jStart[g] ; j < m_jStartIntEQL[g] ; j++ )
    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
    {
      if( field[g](i,j,k) == -1 )
        field[g](i,j,k) = field[g](i,m_jStartIntEQL[g],k);
    }
    #pragma omp parallel for
    for( int k=m_kStart[g] ; k <= m_kEnd[g] ; k++ )
    for( int j=m_jEndInt[g]+1 ; j <= m_jEnd[g] ; j++ )
    for( int i=m_iStart[g] ; i <= m_iEnd[g] ; i++ )
    {
      if( field[g](i,j,k) == -1 )
        field[g](i,j,k) = field[g](i,m_jEndInt[g],k);
    }
   } // end for g
}

//-----------------------------------------------------------------------
void EW::check_materials_EQL()
{
  //---------------------------------------------------------------
  // Verify that the density is nonzero and positive in the 
  // internal grid points
  //---------------------------------------------------------------

   // Minimum allowed  cp/cs, positive definite operator requires cp/cs > sqrt(4/3) = 1.155...
   // lambda >0 requires cp/cs > sqrt(2)
   // const float_sw4 mincpcsratio = 1.2;
  const float_sw4 mincpcsratio = sqrt(4.0/3);
  const float_sw4 la_min_fact = mincpcsratio*mincpcsratio-2;
  
  float_sw4 mins[8],maxs[8];

  float_sw4 lmin = localMin(mRho);
  MPI_Allreduce(&lmin,&mins[0],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
  lmin = localMinVp_EQL();  
  MPI_Allreduce(&lmin,&mins[1],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
  lmin = localMinVs_EQL();  
  MPI_Allreduce(&lmin,&mins[2],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
  lmin = localMin_EQL(mMu);
  MPI_Allreduce(&lmin,&mins[3],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
  lmin = localMin_EQL(mLambda);
  MPI_Allreduce(&lmin,&mins[4],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
  
  CHECK_INPUT(mins[2] >= 0.0,
          "Error: the material data has s velocities that are negative.");

  lmin = localMinVpOverVs_EQL();  
  MPI_Allreduce(&lmin,&mins[5],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);

  float_sw4 lmax = localMax(mRho);
  MPI_Allreduce(&lmax,&maxs[0],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  lmax = localMaxVp_EQL();  
  MPI_Allreduce(&lmax,&maxs[1],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  lmax = localMaxVs_EQL();  
  MPI_Allreduce(&lmax,&maxs[2],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  lmax = localMax_EQL(mMu);
  MPI_Allreduce(&lmax,&maxs[3],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  lmax = localMax_EQL(mLambda);
  MPI_Allreduce(&lmax,&maxs[4],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  lmax = localMaxVpOverVs_EQL();  
  MPI_Allreduce(&lmax,&maxs[5],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);

  // EQL: added EQL specific trackers to list
  if( usingAttenuation() && !m_twilight_forcing)
  {
      lmin = localMin_EQL(mQs);
      MPI_Allreduce(&lmin,&mins[6],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
      lmin = localMin_EQL(mQp);
      MPI_Allreduce(&lmin,&mins[7],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);

      lmax = localMax_EQL(mQs);
      MPI_Allreduce(&lmax,&maxs[6],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
      lmax = localMax_EQL(mQp);
      MPI_Allreduce(&lmax,&maxs[7],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
  }
  
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if( myRank == 0 )
  {
    string indent  = "\n       ";
    string indents =   "       ";
      
    cout << indent << "--------- EQL material properties ranges -------------"
//            << indent << "  For grid [" << m_level << "][" << m_gridnr << "]" << endl
	 << indent << mins[0] << " kg/m^3 <=  Density <= " << maxs[0] << " kg/m^3"
	 << indent << mins[1] << " m/s    <=  Vp      <= " << maxs[1] << " m/s"
	 << indent << mins[5] << "        <=  Vp/Vs   <= " << maxs[5]
	 << indent << mins[4] << " Pa     <=  lambda  <= " << maxs[4] << " Pa" 
   << indent 
   << indent << mins[3] << " Pa     <=  mu      <= " << maxs[3] << " Pa"
   << indent << mins[2] << " m/s    <=  Vs      <= " << maxs[2] << " m/s" 
   << indent << endl;

    if( usingAttenuation() && !m_twilight_forcing)
    {
      cout << indents << "Using attenuation "
           << indent << mins[6] << "        <=  Qs      <= " << maxs[6] << "  "
           << indent << mins[7] << "        <=  Qp      <= " << maxs[7] << "  " << endl;
    }
    cout  << indents << "------------------------------------------------------" << endl;
  }
   
  if( mins[0] <= 0.0 )
  {
    for (int g = 0; g < mNumberOfGrids; g++)
#pragma omp parallel for
      for( int k=m_kStartIntEQL[g] ; k <= m_kEnd[g] ; k++ )
	for( int j=m_jStartIntEQL[g] ; j <= m_jEnd[g] ; j++ )
	  for( int i=m_iStartIntEQL[g] ; i <= m_iEnd[g] ; i++ )
	  {
	    CHECK_INPUT( mRho[g](i,j,k) > 0., "Density= " << mRho[g](i,j,k)<< " in grid g= " << g << " at point " 
			  << " (" << i <<","<<j<<","<<k<<") ");
	  }
  }
   
  if( mins[3] < 0.0 )
  {
    for (int g = 0; g < mNumberOfGrids; g++)
#pragma omp parallel for
      for( int k=m_kStartIntEQL[g] ; k <= m_kEnd[g] ; k++ )
	for( int j=m_jStartIntEQL[g] ; j <= m_jEnd[g] ; j++ )
	  for( int i=m_iStartIntEQL[g] ; i <= m_iEnd[g] ; i++ )
	  {
	    CHECK_INPUT( mMu[g](i,j,k) >= 0., "mu= " << mMu[g](i,j,k)<< " in grid g= " << g << " at point " 
			  << " (" << i <<","<<j<<","<<k<<") ");
	  }
  }
  if( mins[4] <= 0.0 )
  {
    for (int g = 0; g < mNumberOfGrids; g++)
#pragma omp parallel for
      for( int k=m_kStartIntEQL[g] ; k <= m_kEnd[g] ; k++ )
	for( int j=m_jStartIntEQL[g] ; j <= m_jEnd[g] ; j++ )
	  for( int i=m_iStartIntEQL[g] ; i <= m_iEnd[g] ; i++ )
	  {
	     CHECK_INPUT( mLambda[g](i,j,k) >= la_min_fact*mMu[g](i,j,k), "lambda= " << mLambda[g](i,j,k)<< " in grid g= " << g << " at point " 
			 << " (" << i <<","<<j<<","<<k<<") "<<"mMu="<<mMu[g](i,j,k) );
	  }
  }
  if( m_use_attenuation && !m_twilight_forcing)
  {
     if( mins[6] <= 0.0 )
     {
	for (int g = 0; g < mNumberOfGrids; g++)
#pragma omp parallel for
	   for( int k=m_kStartIntEQL[g] ; k <= m_kEnd[g] ; k++ )
	      for( int j=m_jStartIntEQL[g] ; j <= m_jEnd[g] ; j++ )
		 for( int i=m_iStartIntEQL[g] ; i <= m_iEnd[g] ; i++ )
		 {
		    CHECK_INPUT( mQs[g](i,j,k) >= 0., "Qs= " << mQs[g](i,j,k)<< " in grid g= " << g << " at point " 
				 << " (" << i <<","<<j<<","<<k<<") ");
		 }
     }
     if( mins[7] <= 0.0 )
     {
	for (int g = 0; g < mNumberOfGrids; g++)
#pragma omp parallel for
	   for( int k=m_kStartIntEQL[g] ; k <= m_kEnd[g] ; k++ )
	      for( int j=m_jStartIntEQL[g] ; j <= m_jEnd[g] ; j++ )
		 for( int i=m_iStartIntEQL[g] ; i <= m_iEnd[g] ; i++ )
		 {
		    CHECK_INPUT( mQp[g](i,j,k) >= 0., "Qp= " << mQp[g](i,j,k)<< " in grid g= " << g << " at point " 
				 << " (" << i <<","<<j<<","<<k<<") ");
		 }
     }
  }

   CHECK_INPUT(mins[0] > 0.0,
		        "Error: the material data has density values less than or equal to zero.");
   CHECK_INPUT(mins[1] > 0.0,
           "Error: the material data has p velocities that are less than or equal to zero.");
   CHECK_INPUT(mins[3] >= 0.0,
           "Error: mu has values that are negative.");
   //   CHECK_INPUT(mins[4] > 0.0,
   //           "Error: lambda has values that are negative or zero.");
  
  VERIFY2(mins[5] >= mincpcsratio,
	  "Error: vp/vs is smaller than set limit, " << mincpcsratio );

// check material ranges on each grid
   if (mVerbose >= 3)
   {
   float_sw4 minRho, maxRho, minMu, maxMu, minLambda, maxLambda;
   
   for( int g = 0 ; g < mNumberOfGrids; g++)
   {
     minRho=1e38, minMu=1e38, minLambda=1e38;
     maxRho=0, maxMu=0, maxLambda=0;
     
    #pragma omp parallel for reduction(min:minMu,minLambda,minRho) reduction(max:maxMu,maxLambda,maxRho)
    for( int k = m_kStartIntEQL[g] ; k <= m_kEnd[g]; k++ )
      for( int j = m_jStartIntEQL[g] ; j <= m_jEnd[g]; j++ )
      for( int i = m_iStartIntEQL[g] ; i <= m_iEnd[g] ; i++ )
      {
        if (mMu[g](i,j,k) < minMu) minMu=mMu[g](i,j,k);
        if (mMu[g](i,j,k) > maxMu) maxMu=mMu[g](i,j,k);
        if (mLambda[g](i,j,k) < minLambda) minLambda=mLambda[g](i,j,k);
        if (mLambda[g](i,j,k) > maxLambda) maxLambda=mLambda[g](i,j,k);
        if (mRho[g](i,j,k) < minRho) minRho=mRho[g](i,j,k);
        if (mRho[g](i,j,k) > maxRho) maxRho=mRho[g](i,j,k);
      }
// communicate min & max
     MPI_Allreduce(&minRho,&mins[0],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
     MPI_Allreduce(&maxRho,&maxs[0],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
     MPI_Allreduce(&minMu,&mins[1],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
     MPI_Allreduce(&maxMu,&maxs[1],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
     MPI_Allreduce(&minLambda,&mins[2],1,m_mpifloat,MPI_MIN,m_cartesian_communicator);
     MPI_Allreduce(&maxLambda,&maxs[2],1,m_mpifloat,MPI_MAX,m_cartesian_communicator);
// printout results
     if (proc_zero())
     {
       printf("Grid #%i:, %e <= Rho <= %e, %e <= Mu <= %e, %e <= Lambda <= %e\n", g, mins[0], maxs[0], mins[1], maxs[1], mins[2], maxs[2]);
     }
     
   } // end for all grids
   }
// end mVerbose >= 3 

// evaluate min(Cs) and max(sqrt(Cp^2+2*Cs^2))for each grid  
   if (mVerbose >= 1)
   {
     double Cs, C_hat, minCs, maxCs, minC_hat, maxC_hat;
   
     for( int g = 0 ; g < mNumberOfGrids; g++)
     {
       minCs=1e100, minC_hat=1e100;
       maxCs=0, maxC_hat=0;
       
       for( int k = m_kStartIntEQL[g] ; k <= m_kEnd[g]; k++ )
	 for( int j = m_jStartIntEQL[g] ; j <= m_jEnd[g]; j++ )
	   for( int i = m_iStartIntEQL[g] ; i <= m_iEnd[g] ; i++ )
	   {
// take square root after computing the min and max
	     Cs = mMu[g](i,j,k)/mRho[g](i,j,k);
	     C_hat = (mLambda[g](i,j,k) + 4*mMu[g](i,j,k))/mRho[g](i,j,k);

	     
	     if (Cs < minCs) minCs=Cs;
	     if (Cs > maxCs) maxCs=Cs;

	     if (C_hat < minC_hat) minC_hat=C_hat;
	     if (C_hat > maxC_hat) maxC_hat=C_hat;
	   }
       minCs = sqrt(minCs);
       maxCs = sqrt(maxCs);
       minC_hat = sqrt(minC_hat);
       maxC_hat = sqrt(maxC_hat);
       
// communicate min & max
       MPI_Allreduce(&minCs,&mins[0],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
       MPI_Allreduce(&maxCs,&maxs[0],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
       MPI_Allreduce(&minC_hat,&mins[1],1,MPI_DOUBLE,MPI_MIN,m_cartesian_communicator);
       MPI_Allreduce(&maxC_hat,&maxs[1],1,MPI_DOUBLE,MPI_MAX,m_cartesian_communicator);
// printout results
       if (mVerbose >=2 && proc_zero())
       {
	 printf("Material model info, Grid g=%i: %e <= Cs <= %e, %e <= C-hat <= %e, h[g]/max(C-hat) = %e\n",
		g, mins[0], maxs[0], mins[1], maxs[1], mGridSize[g]/maxs[1]);
       }
     } // end for all grids
   }
// end mVerbose >= 1
}

//-----------------------------------------------------------------------
float_sw4 EW::localMin_EQL(std::vector<Sarray> & a_field) 
{
  float_sw4 lmin_all = a_field[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]);
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmin=1e38;
#pragma omp parallel for reduction(min:lmin)
      for (int k = m_kStartIntEQL[g]; k <= m_kEndIntEQL[g]; k++ )
        {
          for (int j = m_jStartIntEQL[g]; j <= m_jEndIntEQL[g]; j++ )
            {      
              for (int i = m_iStartIntEQL[g]; i <= m_iEndIntEQL[g]; i++ )
                {
                  if (a_field[g](i,j,k) < lmin)
                    {
                      lmin = a_field[g](i,j,k);
		      //		      cout << "lmin = " << lmin << " at " << i << " " << j << " " << k << endl;
                    }
                }
            }
        }
      lmin_all = lmin < lmin_all?lmin:lmin_all;
    }

  return lmin_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMax_EQL(std::vector<Sarray> & a_field) 
{
  float_sw4 lmax_all = a_field[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]);
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmax=-1e38;
#pragma omp parallel for reduction(max:lmax)
      for (int k = m_kStartIntEQL[g]; k <= m_kEndIntEQL[g]; k++ )
        {
          for (int j = m_jStartIntEQL[g]; j <= m_jEndIntEQL[g]; j++ )
            {      
              for (int i = m_iStartIntEQL[g]; i <= m_iEndIntEQL[g]; i++ )
                {
                  if (a_field[g](i,j,k) > lmax)
                    {
                      lmax = a_field[g](i,j,k);
                    }
                }
            }
        }
      lmax_all = lmax > lmax_all?lmax:lmax_all;
    }

  return lmax_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMinVp_EQL() 
{
  float_sw4 lmin_all = sqrt((2.*mMu[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0])+mLambda[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]))
			 /mRho[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmin=1e38;
#pragma omp parallel for reduction(min:lmin)
      for (int k = m_kStartIntEQL[g]; k <= m_kEndIntEQL[g]; k++ )
        {
          for (int j = m_jStartIntEQL[g]; j <= m_jEndIntEQL[g]; j++ )
            {      
              for (int i = m_iStartIntEQL[g]; i <= m_iEndIntEQL[g]; i++ )
                {
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k)) < lmin)
                    {
                      lmin = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k));
                    }
                }
            }
        }
      lmin_all = lmin < lmin_all?lmin:lmin_all;
    }
  return lmin_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMaxVp_EQL() 
{
  float_sw4 lmax_all = sqrt((2.*mMu[0](m_iStart[0],m_jStart[0],m_kStart[0])+mLambda[0](m_iStart[0],m_jStart[0],m_kStart[0]))
			 /mRho[0](m_iStart[0],m_jStart[0],m_kStart[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmax=-1e38;
#pragma omp parallel for reduction(max:lmax)
      for (int k = m_kStartIntEQL[g]; k <= m_kEndIntEQL[g]; k++ )
        {
          for (int j = m_jStartIntEQL[g]; j <= m_jEndIntEQL[g]; j++ )
            {      
              for (int i = m_iStartIntEQL[g]; i <= m_iEndIntEQL[g]; i++ )
                {
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k)) > lmax)
                    {
                      lmax = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k));
                    }
                }
            }
        }
      lmax_all = lmax > lmax_all?lmax:lmax_all;
    }

  return lmax_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMinVs_EQL() 
{
  float_sw4 lmin_all = sqrt(mMu[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0])/mRho[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmin=1e38;
#pragma omp parallel for reduction(min:lmin)
      for (int k = m_kStartIntEQL[g]; k <= m_kEndIntEQL[g]; k++ )
        {
          for (int j = m_jStartIntEQL[g]; j <= m_jEndIntEQL[g]; j++ )
            {      
              for (int i = m_iStartIntEQL[g]; i <= m_iEndIntEQL[g]; i++ )
                {
                  if (sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) < lmin)
                    {
                      lmin = sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
      lmin_all = lmin < lmin_all?lmin:lmin_all;
    }

  return lmin_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMaxVs_EQL() 
{
  float_sw4 lmax_all = sqrt(mMu[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0])/mRho[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmax=-1e38;
#pragma omp parallel for reduction(max:lmax)
      for (int k = m_kStartIntEQL[g]; k <= m_kEndIntEQL[g]; k++ )
        {
          for (int j = m_jStartIntEQL[g]; j <= m_jEndIntEQL[g]; j++ )
            {      
              for (int i = m_iStartIntEQL[g]; i <= m_iEndIntEQL[g]; i++ )
                {
                  if ( sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) > lmax)
                    {
                      lmax = sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
      lmax_all = lmax > lmax_all?lmax:lmax_all;
    }
  return lmax_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMinVpOverVs_EQL() 
{
  float_sw4 lmin_all = sqrt((2.*mMu[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0])+mLambda[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]))
			 /mRho[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]))/sqrt(mMu[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0])
									     /mRho[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmin=1e38;
#pragma omp parallel for reduction(min:lmin)
      for (int k = m_kStartIntEQL[g]; k <= m_kEndIntEQL[g]; k++ )
        {
          for (int j = m_jStartIntEQL[g]; j <= m_jEndIntEQL[g]; j++ )
            {      
              for (int i = m_iStartIntEQL[g]; i <= m_iEndIntEQL[g]; i++ )
                {
// Unneccessary to divided by rho in Vp and Vs because it cancels
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) < lmin)
                    {
                      lmin = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
      lmin_all = lmin < lmin_all?lmin:lmin_all;
    }
  return lmin_all; 
}

//-----------------------------------------------------------------------
float_sw4 EW::localMaxVpOverVs_EQL() 
{
  float_sw4 lmax_all = sqrt((2.*mMu[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0])+mLambda[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]))
			 /mRho[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]))/sqrt(mMu[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0])
									     /mRho[0](m_iStartIntEQL[0],m_jStartIntEQL[0],m_kStartIntEQL[0]));
 
  for (int g = 0; g < mNumberOfGrids; g++)
    {
       float_sw4 lmax=-1e38;
#pragma omp parallel for reduction(max:lmax)
      for (int k = m_kStartIntEQL[g]; k <= m_kEndIntEQL[g]; k++ )
        {
          for (int j = m_jStartIntEQL[g]; j <= m_jEndIntEQL[g]; j++ )
            {      
              for (int i = m_iStartIntEQL[g]; i <= m_iEndIntEQL[g]; i++ )
                {
// Unneccessary to divided by rho in Vp and Vs because it cancels
                  if (sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k)) > lmax)
                    {
                      lmax = sqrt((2.*mMu[g](i,j,k)+mLambda[g](i,j,k))/mRho[g](i,j,k))/sqrt(mMu[g](i,j,k)/mRho[g](i,j,k));
                    }
                }
            }
        }
      lmax_all = lmax > lmax_all?lmax:lmax_all;
    }
  return lmax_all; 
}

//---------------------------------------------------------------
// 
// 
//---------------------------------------------------------------
bool EW::updateEmax(vector<Sarray> &U)
{
  int supergrid_pts;
  bool updated = false;
  for( int g = 0 ; g < mNumberOfGrids; g++)
  {
    supergrid_pts = static_cast<int>(m_supergrid_width/mGridSize[g]);
    #pragma omp parallel for     
    for( int k = m_kStartIntEQL[g] ; k <= m_kEndIntEQL[g]; k++ )
    {
      for( int j = m_jStartIntEQL[g] ; j <= m_jEndIntEQL[g]; j++ )
      {
        for( int i = m_iStartIntEQL[g] ; i <= m_iEndIntEQL[g] ; i++ )
        {
          float_sw4 x, y, z;
          float_sw4 ptDepth;
          float_sw4 vs = sqrt(mMuOrig_EQL[g](i, j, k)/mRho[g](i, j, k));
          getCoordinates(i, j, k, g, x, y, z);
          getDepth(x, y, z, ptDepth); // Gets total depth below free surface including topography

          // if the shear wave velocity is greater than vslim,
          // or if too close to the source,
          // or if depth below the set limit,
          // or if its a supergrid point
          // then don't apply eql (aka, set max strain = 0)
          if ( (vs > m_vslim_eql) ||
               (ptDepth > m_max_depth_eql) ||
               (z > m_max_z_eql) ||
               (m_min_dist_to_srcs[g](i,j,k) < m_src_Dmin) ||
               (i <= supergrid_pts) ||
               (i >= m_global_nx[g]-supergrid_pts) ||
               (j <= supergrid_pts) ||
               (j >= m_global_ny[g]-supergrid_pts) ||
               (g == 0 && z >= m_global_nz[g]-supergrid_pts)
             )
          { 
            mEmax[g](i, j, k) = 0.0;
          }
          else
          {
            if (g < mNumberOfCartesianGrids) // must be a Cartesian grid
            {
              float_sw4 factor = 1.0 / (2 * mGridSize[g]);
              float_sw4 duydx = (U[g](2, i + 1, j, k) - U[g](2, i - 1, j, k)) * factor;
              float_sw4 duzdx = (U[g](3, i + 1, j, k) - U[g](3, i - 1, j, k)) * factor;
              float_sw4 duxdy = (U[g](1, i, j + 1, k) - U[g](1, i, j - 1, k)) * factor;
              float_sw4 duzdy = (U[g](3, i, j + 1, k) - U[g](3, i, j - 1, k)) * factor;
              float_sw4 duxdz = (U[g](1, i, j, k + 1) - U[g](1, i, j, k - 1)) * factor;
              float_sw4 duydz = (U[g](2, i, j, k + 1) - U[g](2, i, j, k - 1)) * factor;
              float_sw4 duxdx = (U[g](1, i + 1, j, k) - U[g](1, i - 1, j, k)) * factor;
              float_sw4 duydy = (U[g](2, i, j + 1, k) - U[g](2, i, j - 1, k)) * factor;
              float_sw4 duzdz = (U[g](3, i, j, k + 1) - U[g](3, i, j, k - 1)) * factor;

              // strain tensor values
              float_sw4 E11 = (duxdx);
              float_sw4 E22 = (duydy);
              float_sw4 E33 = (duzdz);
              float_sw4 E12 = (0.5 * (duydx + duxdy));
              float_sw4 E13 = (0.5 * (duzdx + duxdz));
              float_sw4 E23 = (0.5 * (duydz + duzdy));

              // strain tensor invariants
              float_sw4 I1 = (E11 + E22 + E33);
              float_sw4 I2 = (E11 * E22 + E22 * E33 + E33 * E11 - E12 * E12 - E13 * E13 - E23 * E23);
              float_sw4 I3 = (E11 * E22 * E33 - E11 * E23 * E23 - E22 * E13 * E13 - E33 * E12 * E12 + 2.0 * E12 * E13 * E23);
              
              // second invariant of strain deviator tensor
              float_sw4 J2 = (1.0/3.0) * I1*I1 - I2;

              // approximation of max shear strain using J2
              float_sw4 J2_est = sqrt(4*J2);

              if (m_iterLim_EQL != 0)
              {
                if (J2_est > mEmax[g](i, j, k))
                {
                  mEmax[g](i, j, k) = J2_est;
                  updated = true;
                }
              }
              else
                mEmax[g](i, j, k) = J2_est;
            }
            else // must be curvilinear
            {
              //       int i=m_i0, j=m_j0, k=m_k00, g=m_grid0;
              float_sw4 factor = 0.5 / sqrt(mJ[g](i, j, k));
              //      float_sw4 duxdq = (U[g](1,i+1,j,k) - U[g](1,i-1,j,k));
              float_sw4 duydq = (U[g](2, i + 1, j, k) - U[g](2, i - 1, j, k));
              float_sw4 duzdq = (U[g](3, i + 1, j, k) - U[g](3, i - 1, j, k));
              float_sw4 duxdr = (U[g](1, i, j + 1, k) - U[g](1, i, j - 1, k));
              //      float_sw4 duydr = (U[g](2,i,j+1,k) - U[g](2,i,j-1,k));
              float_sw4 duzdr = (U[g](3, i, j + 1, k) - U[g](3, i, j - 1, k));
              float_sw4 duxds = (U[g](1, i, j, k + 1) - U[g](1, i, j, k - 1));
              float_sw4 duyds = (U[g](2, i, j, k + 1) - U[g](2, i, j, k - 1));
              float_sw4 duzds = (U[g](3, i, j, k + 1) - U[g](3, i, j, k - 1));
              float_sw4 duzdy = (mMetric[g](1, i, j, k) * duzdr + mMetric[g](3, i, j, k) * duzds) * factor;
              float_sw4 duydz = (mMetric[g](4, i, j, k) * duyds) * factor;
              float_sw4 duxdz = (mMetric[g](4, i, j, k) * duxds) * factor;
              float_sw4 duzdx = (mMetric[g](1, i, j, k) * duzdq + mMetric[g](2, i, j, k) * duzds) * factor;
              float_sw4 duydx = (mMetric[g](1, i, j, k) * duydq + mMetric[g](2, i, j, k) * duyds) * factor;
              float_sw4 duxdy = (mMetric[g](1, i, j, k) * duxdr + mMetric[g](3, i, j, k) * duxds) * factor;
              float_sw4 duxdx = (mMetric[g](1, i, j, k) * (U[g](1, i + 1, j, k) - U[g](1, i - 1, j, k)) +
                                  mMetric[g](2, i, j, k) * (U[g](1, i, j, k + 1) - U[g](1, i, j, k - 1))) *
                              factor;
              float_sw4 duydy = (mMetric[g](1, i, j, k) * (U[g](2, i, j + 1, k) - U[g](2, i, j - 1, k)) +
                                  mMetric[g](3, i, j, k) * (U[g](2, i, j, k + 1) - U[g](2, i, j, k - 1))) *
                              factor;
              float_sw4 duzdz = (mMetric[g](4, i, j, k) * (U[g](3, i, j, k + 1) - U[g](3, i, j, k - 1))) * factor;

              // strain tensor values
              float_sw4 E11 = (duxdx);
              float_sw4 E22 = (duydy);
              float_sw4 E33 = (duzdz);
              float_sw4 E12 = (0.5 * (duydx + duxdy));
              float_sw4 E13 = (0.5 * (duzdx + duxdz));
              float_sw4 E23 = (0.5 * (duydz + duzdy));

              // invariants
              float_sw4 I1 = (E11 + E22 + E33);
              float_sw4 I2 = (E11 * E22 + E22 * E33 + E33 * E11 - E12 * E12 - E13 * E13 - E23 * E23);
              float_sw4 I3 = (E11 * E22 * E33 - E11 * E23 * E23 - E22 * E13 * E13 - E33 * E12 * E12 + 2.0 * E12 * E13 * E23);

              // second invariant of strain deviator tensor
              float_sw4 J2 = (1.0/3.0) * I1*I1 - I2;

              // approximation of max shear strain using J2
              float_sw4 J2_est = sqrt(4*J2);

              if (m_iterLim_EQL != 0)
              {
                if (J2_est > mEmax[g](i, j, k))
                {
                  mEmax[g](i, j, k) = J2_est;
                  updated = true;
                }
              }
              else
                mEmax[g](i, j, k) = J2_est;
            }
          }
        }
      }
    }
    communicate_array( mEmax[g], g );
  } // end for all grids
  return updated;
} // end updateEmax
