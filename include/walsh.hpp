/*
    walsh.hpp
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

#ifndef WALSH_HPP
#define WALSH_HPP
#include <plugin.h>
void ffwt (csnd::Csound* csound, csnd::Vector<MYFLT> x);
void fwt (csnd::Csound* csound, csnd::Vector<MYFLT> x);
void haar (csnd::Csound* csound, csnd::Vector<MYFLT> x);
void haarin (csnd::Csound* csound, csnd::Vector<MYFLT> x);
void hnorm (csnd::Vector<MYFLT> x);
int i4_log_2 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
void r8vec_copy ( int n, double a1[], csnd::Vector<MYFLT> a2);
void r8vec_copy ( int n, csnd::Vector<MYFLT> a1, double a2[]);
double *r8vec_copy_new (csnd::Csound* csound, int n, double a1[] );
void walsh (csnd::Csound* csound, csnd::Vector<MYFLT> x);


#endif /* WALSH_HPP */

