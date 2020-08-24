# Introduction

These opcodes present variations on the following transforms:
    Haar, 
    Walsh/Hadamard 
    Cosine

At current the opcodes will function simply and similar to the fft and fftinv 
opcodes in Csound - ie, the user must implement windowing and ksmps in respect 
to the desired analysis.

These are experimental and exact usage/manipulation of the functions is subject
to some testing so there may be further development to make the opcodes more
usable.

Tested on Linux with Csound 6.14


# Installation

    mkdir build && cd build
    cmake ..
    make
    sudo make install
     

# Opcodes provided
Each of the following opcodes receives and emits a k-rate array. The size of the
input array should be a power of two. The output size will match the input.


### tfcosine
Cosine transform. Notably slow (probably not for realtime usage)

### tfcosineinv
Inverse cosine transform. As Above

### tfhaar1
Haar transform, algorithm 1

### tfhaar1inv
Inverse Haar transform, algorithm 1

### tfhaar2
Haar transform, algorithm 2

### tfhaar2inv
Inverse Haar transform, algorithm 1

### tfwalsh1
Walsh/Hadamard transform, algorithm 1

### tfwalsh2
Walsh/Hadamard transform, algorithm 2

### tfwalsh2inv
Inverse Walsh/Hadamard transform, algorithm 2




# Credits

Code, inspiration and relevant licensing has been derived from the following:
    https://people.sc.fsu.edu/~jburkardt/f_src/walsh/walsh.html
    https://people.sc.fsu.edu/~jburkardt/f_src/haar/haar.html
    https://people.sc.fsu.edu/~jburkardt/f_src/cosine_transform/cosine_transform.html
    https://github.com/mochow13/competitive-programming-library/blob/master/Math/Fast%20Walsh-Hadamard%20Transform.cpp


Todo / possibly forthcoming
    Implement Daubechies wavelets, sine transform.

