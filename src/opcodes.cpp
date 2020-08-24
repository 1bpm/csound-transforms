/*
    opcodes.cpp
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
#include "cosine.hpp"
#include "haar.hpp"
#include "walsh.hpp"
#include "walsh2.hpp"


struct tfcosine : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "k[]";
    
    int init() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        out.init(csound, in.len());
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        cosine_transform_data(in, out);
        return OK;
    }
};

struct tfcosineinv : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "k[]";
    
    int init() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        out.init(csound, in.len());
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        cosine_transform_inverse(in, out);
        return OK;
    }
};

struct tfhaar1 : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "k[]";
    
    int init() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        out.init(csound, in.len());
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        std::copy(in.begin(), in.end(), out.begin());
        haar_1d(csound, out);
        return OK;
    }
};

struct tfhaar1inv : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "k[]";
    
    int init() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        out.init(csound, in.len());
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        std::copy(in.begin(), in.end(), out.begin());
        haar_1d_inverse(csound, out);
        return OK;
    }
};

struct tfhaar2 : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "k[]";
    
    int init() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        out.init(csound, in.len());
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        std::copy(in.begin(), in.end(), out.begin());
        haar(csound, out);
        return OK;
    }
};

struct tfhaar2inv : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "k[]";
    
    int init() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        out.init(csound, in.len());
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        std::copy(in.begin(), in.end(), out.begin());
        haarin(csound, out);
        return OK;
    }
};

struct tfwalsh1 : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "k[]";
    
    int init() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        out.init(csound, in.len());
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        std::copy(in.begin(), in.end(), out.begin());
        walsh(csound, out);
        return OK;
    }
};

struct tfwalsh2 : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "k[]";
    
    int init() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        out.init(csound, in.len());
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        std::copy(in.begin(), in.end(), out.begin());
        fwtanal(out);
        return OK;
    }
};

struct tfwalsh2inv : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "k[]";
    
    int init() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        out.init(csound, in.len());
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        std::copy(in.begin(), in.end(), out.begin());
        fwtsynth(out);
        return OK;
    }
};


#include <modload.h>

void csnd::on_load(csnd::Csound *csound) {
    csnd::plugin<tfcosine>(csound, "tfcosine", csnd::thread::ik);
    csnd::plugin<tfcosineinv>(csound, "tfcosineinv", csnd::thread::ik);
    csnd::plugin<tfhaar1>(csound, "tfhaar1", csnd::thread::ik);
    csnd::plugin<tfhaar1inv>(csound, "tfhaar1inv", csnd::thread::ik);
    csnd::plugin<tfhaar2>(csound, "tfhaar2", csnd::thread::ik);
    csnd::plugin<tfhaar2inv>(csound, "tfhaar2inv", csnd::thread::ik);
    csnd::plugin<tfwalsh1>(csound, "tfwalsh1", csnd::thread::ik);
    csnd::plugin<tfwalsh2>(csound, "tfwalsh2", csnd::thread::ik);
    csnd::plugin<tfwalsh2inv>(csound, "tfwalsh2inv", csnd::thread::ik);
}
