const char* dgemm_desc = "Dlaiton, best way dgemm.";
#include <math.h>
#include <stdlib.h>

void square_dgemm(int n, double* A, double* B, double* C){
  register int i,j,k,jn,kjn,kn;
  for(j=0;j<n;j++){
      jn=j*n;
      for(k=0;k<n;k++){
        kjn=k+jn;
        kn=k*n;
        for(i=0;i<n;i++){
          C[i+jn]+=A[i+kn]*B[kjn];        
        }
      }
    }
}