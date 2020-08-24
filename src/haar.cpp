/*
    haar.cpp
    Copyright (C) 2014 John Burkardt
    Copyright (C) 2019 Richard Knight


    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software Foundation,
    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <plugin.h>

using namespace std;

#include "haar.hpp"

//****************************************************************************80

void haar_1d (csnd::Csound* csound, csnd::Vector<MYFLT> x)

//****************************************************************************80
//
//  Purpose:
//
//    HAAR_1D computes the Haar transform of a vector.
//
//  Discussion:
//
//    For the classical Haar transform, N should be a power of 2.
//    However, this is not required here.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vector.
//
//    Input/output, double X[N], on input, the vector to be transformed.
//    On output, the transformed vector.
//
{ 
  int n = x.len();
  int i;
  int k;
  double s;
  double *y;

  s = sqrt ( 2.0 );

  //y = new double[n];
  y = (double*) csound->malloc(sizeof(double) * n);
//
//  Initialize.
//
  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
  }
//
//  Determine K, the largest power of 2 such that K <= N.
//
  k = 1;
  while ( k * 2 <= n )
  {
    k = k * 2;
  }
  
  while ( 1 < k )
  {
    k = k / 2;
    for ( i = 0; i < k; i++ )
    {
      y[i]   = ( x[2*i] + x[2*i+1] ) / s;
      y[i+k] = ( x[2*i] - x[2*i+1] ) / s;
    }
    for ( i = 0; i < k * 2; i++ )
    {
      x[i] = y[i];
    }
  }

  //delete [] y;
  csound->free(y);
  
  return;
}
//****************************************************************************80

void haar_1d_inverse (csnd::Csound* csound, csnd::Vector<MYFLT> x)

//****************************************************************************80
//
//  Purpose:
//
//    HAAR_1D_INVERSE computes the inverse Haar transform of a vector.
//
//  Discussion:
//
//    For the classical Haar transform, N should be a power of 2.
//    However, this is not required here.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vector.  
//
//    Input/output, double X[N], on input, the vector to be transformed.
//    On output, the transformed vector.
//
{
    
  int n = x.len(); 
  int i;
  int k;
  double s;
  double *y;

  s = sqrt ( 2.0 );

  //y = new double[n];
  y = (double*) csound->malloc(sizeof(double) * n);
//
//  Initialize.
//
  for ( i = 0; i < n; i++ )
  {
    y[i] = 0.0;
  }

  k = 1;
  while ( k * 2 <= n )
  {
    for ( i = 0; i < k; i++ )
    {
      y[2*i]   = ( x[i] + x[i+k] ) / s;
      y[2*i+1] = ( x[i] - x[i+k] ) / s;
    }
    for ( i = 0; i < k * 2; i++ )
    {
      x[i] = y[i];
    }
    k = k * 2;
  }

 //delete [] y;
  csound->free(y);
  return;
}



/*


//****************************************************************************80

void haar_2d ( int m, int n, double u[] )

//****************************************************************************80
//
//  Purpose:
//
//    HAAR_2D computes the Haar transform of an array.
//
//  Discussion:
//
//    For the classical Haar transform, M and N should be a power of 2.
//    However, this is not required here.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the dimensions of the array.
//
//    Input/output, double U[M*N], the array to be transformed.
//
{
  int i;
  int j;
  int k;
  double s;
  double *v;

  s = sqrt ( 2.0 );

  v = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = u[i+j*m];
    }
  }
//
//  Determine K, the largest power of 2 such that K <= M.
//
  k = 1;
  while ( k * 2 <= m )
  {
    k = k * 2;
  }
//
//  Transform all columns.
//
  while ( 1 < k )
  {
    k = k / 2;

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < k; i++ )
      {
        v[i  +j*m] = ( u[2*i+j*m] + u[2*i+1+j*m] ) / s;
        v[k+i+j*m] = ( u[2*i+j*m] - u[2*i+1+j*m] ) / s;
      }
    }
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < 2 * k; i++ )
      {
        u[i+j*m] = v[i+j*m];
      }
    }
  }
//
//  Determine K, the largest power of 2 such that K <= N.
//
  k = 1;
  while ( k * 2 <= n )
  {
    k = k * 2;
  }
//
//  Transform all rows.
//
  while ( 1 < k )
  { 
    k = k / 2;

    for ( j = 0; j < k; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        v[i+(  j)*m] = ( u[i+2*j*m] + u[i+(2*j+1)*m] ) / s;
        v[i+(k+j)*m] = ( u[i+2*j*m] - u[i+(2*j+1)*m] ) / s;
      }
    }

    for ( j = 0; j < 2 * k; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        u[i+j*m] = v[i+j*m];
      }
    }
  }
  delete [] v;

  return;
}
//****************************************************************************80

void haar_2d_inverse ( int m, int n, double u[] )

//****************************************************************************80
//
//  Purpose:
//
//    HAAR_2D_INVERSE inverts the Haar transform of an array.
//
//  Discussion:
//
//    For the classical Haar transform, M and N should be a power of 2.
//    However, this is not required here.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the dimensions of the array.
//
//    Input/output, double U[M*N], the array to be transformed.
//
{
  int i;
  int j;
  int k;
  double s;
  double *v;

  s = sqrt ( 2.0 );

  v = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = u[i+j*m];
    }
  }
//
//  Inverse transform of all rows.
//
  k = 1;

  while ( k * 2 <= n )
  {
    for ( j = 0; j < k; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        v[i+(2*j  )*m] = ( u[i+j*m] + u[i+(k+j)*m] ) / s;
        v[i+(2*j+1)*m] = ( u[i+j*m] - u[i+(k+j)*m] ) / s;
      }
    }

    for ( j = 0; j < 2 * k; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        u[i+j*m] = v[i+j*m];
      }
    }
    k = k * 2;
  }
//
//  Inverse transform of all columns.
//
  k = 1;

  while ( k * 2 <= m )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < k; i++ )
      {
        v[2*i  +j*m] = ( u[i+j*m] + u[k+i+j*m] ) / s;
        v[2*i+1+j*m] = ( u[i+j*m] - u[k+i+j*m] ) / s;
      }
    }

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < 2 * k; i++ )
      {
        u[i+j*m] = v[i+j*m];
      }
    }
    k = k * 2;
  }
  delete [] v;

  return;
}

 * 
 */