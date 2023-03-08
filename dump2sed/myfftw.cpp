#include "myfftw.h"
#include <iostream>

namespace MYFFTW{
    
    void fftw1d(int n1, double *in1, double *out1){
        
        fftw_complex *in, *out;
        int i, j, k, num;
        
        //---- prepare matrices
        in = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) *(n1));
        out = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) *(n1));
        fftw_plan plan = fftw_plan_dft_1d(n1, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        //--- input data
        for(i=0; i<n1; i++){
            in[i][0] = in1[i];
            in[i][1] = 0.;
        }
        
        //----- 1D Fourier transform
        fftw_execute(plan);
        
        //---- output
        for(i=0;i<n1; i++)
            out1[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
        
        //---- delete matrices
        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);
    
    }
    
    void fftw3d(int n1, int n2, int n3, double *in1, double *out1)
    {
        fftw_complex *in, *out;
        fftw_plan plan;
        int i, j, k, num;
        
        //---- prepare matrices
        in = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) *(n1*n2*n3));
        out = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) *(n1*n2*n3));
        plan = fftw_plan_dft_3d(n1, n2, n3, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        //--- input data
        for(i=0; i<n1*n2*n3; i++){
            in[i][0] = in1[i];
            in[i][1] = 0.;
        }
        
        //----- 3D Fourier transform
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        
        //---- output
        for(i=0;i<n1*n2*n3; i++) out1[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
        //for(i=0;i<n1*n2*n3; i++) out1[i] = out[i][0];

        //---- delete matrices
        fftw_free(in);
        fftw_free(out);
    }

    int main(){
        int i, j, k;
        int n1, n2, n3, num1, num2;
        n1 = n2 = n3 = 16;

        double *in  = new double[n1 * n2 * n3];    // should be deleted
        double *out = new double[n1 * n2 * n3];    // should be deleted

        //--- make input data
        for(i=0; i<n1; i++){
            for(j=0; j<n2; j++){
              for(k=0; k<n3; k++){
                num1 = i * n2 * n3 + j * n3 + k;
                in[num1] = cos(2.*M_PI * (6.*i/n1 + 8.*j/n2 + 4.*k/n3));
                //in[num1] += cos(2.*M_PI * (5.*i/n1 + 2.*j/n2 + 10.*k/n3));
              }
            }
        }
        
        //--- 3D Fourier transform
        fftw3d(n1, n2, n3, in, out);
        
        //----- Output results
        char ofile[] = "result3.txt";
        FILE *fp = fopen(ofile, "w");
        fprintf(fp, "#(i, j, k) IN OUT\n");
        for(i=0; i < n1; i++){
            for(j=0; j < n2; j++){
              for(k=0; k < n3; k++){
                num1 = i*n2*n3 + j*n3 + k;
                if(fabs(out[num1])>0.00001){
                  printf("%3d\t%3d\t%3d \t%lf \t%lf \n", i, j, k, in[num1], out[num1]);
                  fprintf(fp, "%3d\t%3d\t%3d \t%lf \t%lf \n", i, j, k, in[num1], out[num1]);
                }
              }
              //fprintf(fp, "\n");
            }
        }
        fclose(fp);
        printf("Output: %s\n", ofile);

        return 0;
    }
}

