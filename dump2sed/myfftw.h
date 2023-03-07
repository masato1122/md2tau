#include <stdio.h>
#include <stdlib.h>
#include "fftw3.h"
#include <math.h>

namespace MYFFTW{
    
    void fftw1d(int n1, double *in1, double *out1);
    void fftw3d(int n1, int n2, int n3, double *in1, double *out1);

}

