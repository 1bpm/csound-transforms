/*
    cosine.cpp
    Copyright (C) 2015 John Burkardt
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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <plugin.h>

using namespace std;

#include "cosine.hpp"

//****************************************************************************80

void cosine_transform_data (csnd::Vector<MYFLT> d, csnd::Vector<MYFLT> c)

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_TRANSFORM_DATA does a cosine transform on a vector of data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of data points.
//
//    Input, double D[N], the vector of data.
//
//    Output, double COSINE_TRANSFORM_DATA[N], the transform coefficients.
//
{
  int n = d.len();
  double angle;
  int i;
  int j;
  const double r8_pi = 3.141592653589793;


  for ( i = 0; i < n; i++ )
  {
    c[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      angle = r8_pi * ( double ) ( i * ( 2 * j + 1 ) ) / ( double ) ( 2 * n );
      c[i] = c[i] + cos ( angle ) * d[j];
    }
    c[i] = c[i] * sqrt ( 2.0 / ( double ) ( n ) );
  }
}
//****************************************************************************80

void cosine_transform_inverse (csnd::Vector<MYFLT> c, csnd::Vector<MYFLT> d)

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_TRANSFORM_INVERSE does an inverse cosine transform.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of data points.
//
//    Input, double C[N], the vector of transform coefficients
//
//    Output, double COSINE_TRANSFORM_INVERSE[N], the original data.
//
{
  int n = c.len();
  double angle;
  int i;
  int j;
  double r8_pi = 3.141592653589793;


  for ( i = 0; i < n; i++ )
  {
    d[i] = c[0] / 2.0;
    for ( j = 1; j < n; j++ )
    {
      angle = r8_pi * ( double ) ( ( 2 * i + 1 ) * j ) / ( double ) ( 2 * n );
      d[i] = d[i] + cos ( angle ) * c[j];
    }
    d[i] = d[i] * sqrt ( 2.0 / ( double ) ( n ) );
  }
}
//****************************************************************************80
