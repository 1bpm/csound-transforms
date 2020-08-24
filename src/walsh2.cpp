/*
    walsh2.cpp
    Copyright (C) 2018 Mottakin Chowdhury
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
#include <plugin.h>

const int N = 1<<16;


void fwtanal(csnd::Vector<MYFLT> io) {
    int n = io.len();
    for (int d = 1; d < n; d <<= 1) {
        for (int i = 0, m = d<<1; i < n; i += m) {
            for (int j = 0; j < d; j++) { /// Don't forget modulo if required
                MYFLT x = io[i+j], y = io[i+j+d];
                io[i+j] = (x+y), io[i+j+d] = (x-y);	// xor
                // io[i+j] = x+y; // and
                // io[i+j+d] = x+y; // or
            }
        }
    }
}


void fwtsynth(csnd::Vector<MYFLT> io) {
    int n = io.len();
    for (int d = 1; d < n; d <<= 1) {
        for (int i = 0, m = d<<1; i < n; i += m) {
            for (int j = 0; j < d; j++) { /// Don't forget modulo if required
                MYFLT x = io[i+j]; 
                MYFLT y = io[i+j+d];
                 /// Modular inverse if required here
                ///////io[i+j] = (x+y)>>1, io[i+j+d] = (x-y)>>1; // xor
                io[i+j] = (x+y)/2;
                io[i+j+d] = (x-y)/2;
                // io[i+j] = x-y; // and
                // io[i+j+d] = y-x; // or
            }
        }
    }
}