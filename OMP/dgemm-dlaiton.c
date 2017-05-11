const char* dgemm_desc = "Dlaiton, best way dgemm.";
#include <math.h>
#include <stdlib.h>
#include <omp.h>
/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. 
  //A+=B (MAIS RAPIDO)
  //A=A*B (MAIS RAPIDO)
  //float (MAIS RAPIDO)
  //A*=1 (MAIS RAPIDO)
  //var externa (MAIS RAPIDO)
  //normal (MAIS RAPIDO)
  //salvar (MUITO MAIS RAPIDO)
 */
/*    
int cacheDouble;
void fillWithZeros(double* O,double* N,int o, int n){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      if(i<o&&j<o)
        N[i+j*n]=O[i+j*o];
      else
        N[i+j*n]=0;
    }
}

void implicito(int n, double* A, double* B, double* C){
  int i,j,k;
  for(j=0;j<n;j++){
    for(k=0;k<n;k++){
       for(i=0;i<n;i++){
          C[i+j*n]+=A[i+k*n]*B[k+j*n];        
        }
      }
    }
}

void sum(double* A, double* B, double* C, int n, int signal){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      C[i+j*n]=A[i+j*n]+B[i+j*n]*signal;

}
void mult(double* A, double* B, double* C,int n){
  int i,j,k;
  double t;
  for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        t=C[i+j*n];
        for(k=0;k<n;k++){
          t+=A[i+k*n]*B[k+j*n];        
        }
        C[i+j*n]=t;
      }
    }
}

void strassen(int n, double* A, double* B, double* C){
  //base
  if(n==1){
    C[0]=A[0]*B[0];
    return;
  }
  //spliting
  int ns2=n/2,i,j;
  int ns22=ns2*ns2;
  int sizeofns22=ns22*sizeof(double);
  double* A11=(double*)malloc(sizeofns22),*A12=(double*)malloc(sizeofns22),*A21=(double*)malloc(sizeofns22),*A22=(double*)malloc(sizeofns22);
  double* B11=(double*)malloc(sizeofns22),*B12=(double*)malloc(sizeofns22),*B21=(double*)malloc(sizeofns22),*B22=(double*)malloc(sizeofns22);
  double* C11=(double*)malloc(sizeofns22),*C12=(double*)malloc(sizeofns22),*C21=(double*)malloc(sizeofns22),*C22=(double*)malloc(sizeofns22);
  for(i=0;i<ns2;i++)
    for(j=0;j<ns2;j++){
      A11[i+j*ns2]=A[(i)+(j)*n];
      A12[i+j*ns2]=A[(i)+(j+ns2)*n];
      A21[i+j*ns2]=A[(i+ns2)+(j)*n];
      A22[i+j*ns2]=A[(i+ns2)+(j+ns2)*n];
      B11[i+j*ns2]=B[(i)+(j)*n];
      B12[i+j*ns2]=B[(i)+(j+ns2)*n];
      B21[i+j*ns2]=B[(i+ns2)+(j)*n];
      B22[i+j*ns2]=B[(i+ns2)+(j+ns2)*n];
    }
  double* M1=(double*)malloc(sizeofns22),*M2=(double*)malloc(sizeofns22),*M3=(double*)malloc(sizeofns22),*M4=(double*)malloc(sizeofns22),*M5=(double*)malloc(sizeofns22),*M6=(double*)malloc(sizeofns22),*M7=(double*)malloc(sizeofns22);
  double* T1=(double*)malloc(sizeofns22),*T2=(double*)malloc(sizeofns22);
  //coefficients
  if(ns2<=cacheDouble){
    sum(A11,A22,T1,ns2,1);
    sum(B11,B22,T2,ns2,1);
    implicito(ns2,T1,T2,M1);
    sum(A21,A22,T1,ns2,1);
    implicito(ns2,T1,B11,M2);
    sum(B12,B22,T1,ns2,-1);
    implicito(ns2,A11,T1,M3);
    sum(B21,B11,T1,ns2,-1);
    implicito(ns2,A22,T1,M4);
    sum(A11,A12,T1,ns2,1);
    implicito(ns2,T1,B22,M5);
    sum(A21,A11,T1,ns2,-1);
    sum(B11,B12,T2,ns2,1);
    implicito(ns2,T1,T2,M6);
    sum(A12,A22,T1,ns2,-1);
    sum(B21,B22,T2,ns2,1);
    implicito(ns2,T1,T2,M7);

  }else{
    sum(A11,A22,T1,ns2,1);
    sum(B11,B22,T2,ns2,1);
    strassen(ns2,T1,T2,M1);
    sum(A21,A22,T1,ns2,1);
    strassen(ns2,T1,B11,M2);
    sum(B12,B22,T1,ns2,-1);
    strassen(ns2,A11,T1,M3);
    sum(B21,B11,T1,ns2,-1);
    strassen(ns2,A22,T1,M4);
    sum(A11,A12,T1,ns2,1);
    strassen(ns2,T1,B22,M5);
    sum(A21,A11,T1,ns2,-1);
    sum(B11,B12,T2,ns2,1);
    strassen(ns2,T1,T2,M6);
    sum(A12,A22,T1,ns2,-1);
    sum(B21,B22,T2,ns2,1);
    strassen(ns2,T1,T2,M7);
  }
  //results
  sum(M1,M4,T1,ns2,1);
  sum(M7,M5,T2,ns2,-1);
  sum(T1,T2,C11,ns2,1);
  sum(M3,M5,C12,ns2,1);
  sum(M2,M4,C21,ns2,1);
  sum(M1,M2,T1,ns2,-1);
  sum(M3,M6,T2,ns2,1);
  sum(T1,T2,C22,ns2,1);
  //compiling
  for(i=0;i<ns2;i++)
    for(j=0;j<ns2;j++){
      C[i+j*n]=C11[i+j*ns2];
      C[(i)+(j+ns2)*n]=C12[i+j*ns2];
      C[(i+ns2)+(j)*n]=C21[i+j*ns2];
      C[(i+ns2)+(j+ns2)*n]=C22[i+j*ns2];
    }

}



void square_dgemm(int n, double* A, double* B, double* C){
  int cacheLine = 64; //Linha de cache em Kb
  cacheDouble = (cacheLine*1000)/sizeof(double); //Quantos doubles cabem na linha
  int newSize=pow(2,ceil(log(n)/log(2))),i,j;
  double AS[newSize*newSize], BS[newSize*newSize], CS[newSize*newSize];
  fillWithZeros(A,AS,n,newSize);
  fillWithZeros(B,BS,n,newSize);
  fillWithZeros(C,CS,n,newSize);
  strassen(newSize,AS,BS,CS);
  for(i=0;i<n;i++)
  for(j=0;j<n;j++)
    C[i+j*n]=CS[i+j*newSize];
}

*/
//------------------------------------------------------------//
//--Implicito-~=6,77%
 // -1-ijk-~=3.17238%
   //   -ikj-~=5.62568%
   //   -jik-~=3.11340%
   //   -jki-~=11.6985%
   //   -kji-~=11.6412%
   //   -kij-~=5.77109%
void square_dgemm(int n, double* A, double* B, double* C){
int test=1, j, k;
omp_set_num_threads(test);  
for(k=0;k<n;k++){
      for(j=0;j<n;j++){
        #pragma omp parallel 
	{
	  for(int i=omp_get_thread_num();i<n;i+=test)
       		C[i+j*n]+=A[i+k*n]*B[k+j*n];  
	}
      }
  }
}
/*
void square_dgemm(int n, double* A, double* B, double* C){
  int i,j,k;
  double t;
  for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        t=C[i+j*n];
        for(k=0;k<n;k++){
          t+=A[i+k*n]*B[k+j*n];        
        }
        C[i+j*n]=t;
      }
    }
}
//------------------------------------------------------------//
//--Explicito-
      -ijk-~=2.89752%
      -ikj-~=7.17026%
      -jik-~=2.92621%
      -jki-~=10.9781%
      -kji-~=10.9628%
      -kij-~=6.46127%
void square_dgemm(int n, double* A, double* B, double* C){
  int i,j,k;
  int n2=n*n;
  double *BT=malloc(n2*sizeof(double));
  for(i=0;i<n;i++)
      for(j=0;j<n;j++)
        BT[i+j*n]=B[i*n+j];
  double t;
  for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        t=C[i+j*n];
        for(k=0;k<n;k++){
          t+=A[i+k*n]*BT[j+k*n];        
        }
        C[i+j*n]=t;
      }
    }
    free(BT);
  }
//------------------------------------------------------------//
//--Explicito Dinamico-
      -ijk-~=2.82894%
      -ikj-~=6.0828%
      -jik-~=2.8175%
      -jki-~=8.19757%
      -kji-~=7.92369%
      -kij-~=5.45646%
void square_dgemm(int n, double* A, double* B, double* C){
  int i,j,k;
  int n2=n*n;
  double *BT=malloc(n2*sizeof(double));
  double t;
  for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        t=C[i+j*n];
        for(k=0;k<n;k++){
          if(i==0){
            BT[j+k*n]=B[k+j*n];
          }
          t+=A[i+k*n]*BT[j+k*n];        
        }
        C[i+j*n]=t;
      }
    }
    free(BT);
  }
//------------------------------------------------------------//
//--Camadas-
      -ijk-~=0.668939%
      -ikj-~=0.538015%
      -jik-~=1.27227%
      -jki-~=1.17059%
      -kji-~=0.490923%
      -kij-~=0.522082%
void square_dgemm(int n, double* A, double* B, double* C){
  int i,j,k;
  double a,b;
  int n3=n*n*n;
  int n2=n*n;
  double *S=calloc(n3,sizeof(double));
  for(j=0;j<n;j++)
    for(i=0;i<n;i++){
      a=A[i+j*n],b=B[j+i*n];
      for(k=0;k<n;k++){
        if(S[i+k*n+j*n2]!=0){
          S[i+k*n+j*n2]*=a;//second
          C[i+k*n]+=S[i+k*n+j*n2];//soma no valor final
        }else
          S[i+k*n+j*n2]=a;//first

        if(S[k+i*n+j*n2]!=0){
          S[k+i*n+j*n2]*=b;//second
          C[k+i*n]+=S[k+i*n+j*n2];
        }else
          S[k+i*n+j*n2]=b;//first
      }
    }
  free(S);
}
//------------------------------------------------------------//
//--Achatando-
      -ijk-~=0.534283%
      -ikj-~=0.478002%
      -jik-~=0.815438%
      -jki-~=0.816721%
      -kji-~=0.493425%
      -kij-~=0.467199%
void square_dgemm (int n, double* A, double* B, double* C){
  int i,j,k;
  int n3=n*n*n;
  int n2=n*n;
  double *S=malloc(n3*sizeof(double));
  double a,b;
  for(i=0;i<n3;i++)
    S[i]=1;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      a=A[i+j*n];
      b=B[j+i*n];
      for(k=0;k<n;k++){
          S[i+k*n+j*n2]*=a;
          S[k+i*n+j*n2]*=b;
      }
    }
  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      C[i+j*n]=0;
      for(k=0;k<n;k++){
        C[i+j*n]+=S[i+j*n+k*n2];
    }
  }
  free(S);
}
//------------------------------------------------------------//
//--Strassen- ~=0.055%
void fillWithZeros(double* O,double* N,int o, int n){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      if(i<o&&j<o)
        N[i+j*n]=O[i+j*o];
      else
        N[i+j*n]=0;
    }
}

void sum(double* A, double* B, double* C, int n, int signal){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      C[i+j*n]=A[i+j*n]+B[i+j*n]*signal;

}
void mult(double* A, double* B, double* C,int n){
  int i,j,k;
  double t;
  for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        t=C[i+j*n];
        for(k=0;k<n;k++){
          t+=A[i+k*n]*B[k+j*n];        
        }
        C[i+j*n]=t;
      }
    }
}

void strassen(int n, double* A, double* B, double* C){
  //base
  if(n==1){
    C[0]=A[0]*B[0];
    return;
  }
  //spliting
  int ns2=n/2,i,j;
  int ns22=ns2*ns2;
  int sizeofns22=ns22*sizeof(double);
  double* A11=(double*)malloc(sizeofns22),*A12=(double*)malloc(sizeofns22),*A21=(double*)malloc(sizeofns22),*A22=(double*)malloc(sizeofns22);
  double* B11=(double*)malloc(sizeofns22),*B12=(double*)malloc(sizeofns22),*B21=(double*)malloc(sizeofns22),*B22=(double*)malloc(sizeofns22);
  double* C11=(double*)malloc(sizeofns22),*C12=(double*)malloc(sizeofns22),*C21=(double*)malloc(sizeofns22),*C22=(double*)malloc(sizeofns22);
  for(i=0;i<ns2;i++)
    for(j=0;j<ns2;j++){
      A11[i+j*ns2]=A[(i)+(j)*n];
      A12[i+j*ns2]=A[(i)+(j+ns2)*n];
      A21[i+j*ns2]=A[(i+ns2)+(j)*n];
      A22[i+j*ns2]=A[(i+ns2)+(j+ns2)*n];
      B11[i+j*ns2]=B[(i)+(j)*n];
      B12[i+j*ns2]=B[(i)+(j+ns2)*n];
      B21[i+j*ns2]=B[(i+ns2)+(j)*n];
      B22[i+j*ns2]=B[(i+ns2)+(j+ns2)*n];
    }
  double* M1=(double*)malloc(sizeofns22),*M2=(double*)malloc(sizeofns22),*M3=(double*)malloc(sizeofns22),*M4=(double*)malloc(sizeofns22),*M5=(double*)malloc(sizeofns22),*M6=(double*)malloc(sizeofns22),*M7=(double*)malloc(sizeofns22);
  double* T1=(double*)malloc(sizeofns22),*T2=(double*)malloc(sizeofns22);
  //coefficients
  sum(A11,A22,T1,ns2,1);
  sum(B11,B22,T2,ns2,1);
  strassen(ns2,T1,T2,M1);
  sum(A21,A22,T1,ns2,1);
  strassen(ns2,T1,B11,M2);
  sum(B12,B22,T1,ns2,-1);
  strassen(ns2,A11,T1,M3);
  sum(B21,B11,T1,ns2,-1);
  strassen(ns2,A22,T1,M4);
  sum(A11,A12,T1,ns2,1);
  strassen(ns2,T1,B22,M5);
  sum(A21,A11,T1,ns2,-1);
  sum(B11,B12,T2,ns2,1);
  strassen(ns2,T1,T2,M6);
  sum(A12,A22,T1,ns2,-1);
  sum(B21,B22,T2,ns2,1);
  strassen(ns2,T1,T2,M7);
  //results
  sum(M1,M4,T1,ns2,1);
  sum(M7,M5,T2,ns2,-1);
  sum(T1,T2,C11,ns2,1);
  sum(M3,M5,C12,ns2,1);
  sum(M2,M4,C21,ns2,1);
  sum(M1,M2,T1,ns2,-1);
  sum(M3,M6,T2,ns2,1);
  sum(T1,T2,C22,ns2,1);
  //compiling
  for(i=0;i<ns2;i++)
    for(j=0;j<ns2;j++){
      C[i+j*n]=C11[i+j*ns2];
      C[(i)+(j+ns2)*n]=C12[i+j*ns2];
      C[(i+ns2)+(j)*n]=C21[i+j*ns2];
      C[(i+ns2)+(j+ns2)*n]=C22[i+j*ns2];
    }

}

void square_dgemm(int n, double* A, double* B, double* C){
  int newSize=pow(2,ceil(log(n)/log(2))),i,j;
  double AS[newSize*newSize], BS[newSize*newSize], CS[newSize*newSize];
  fillWithZeros(A,AS,n,newSize);
  fillWithZeros(B,BS,n,newSize);
  fillWithZeros(C,CS,n,newSize);
  strassen(newSize,AS,BS,CS);
  for(i=0;i<n;i++)
  for(j=0;j<n;j++)
    C[i+j*n]=CS[i+j*newSize];
}
//------------------------------------------------------------//
//--New Strassen- ~=4.860%
int cacheDouble;
void fillWithZeros(double* O,double* N,int o, int n){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      if(i<o&&j<o)
        N[i+j*n]=O[i+j*o];
      else
        N[i+j*n]=0;
    }
}

void implicito(int n, double* A, double* B, double* C){
  int i,j,k;
  for(j=0;j<n;j++){
    for(k=0;k<n;k++){
       for(i=0;i<n;i++){
          C[i+j*n]+=A[i+k*n]*B[k+j*n];        
        }
      }
    }
}

void sum(double* A, double* B, double* C, int n, int signal){
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      C[i+j*n]=A[i+j*n]+B[i+j*n]*signal;

}
void mult(double* A, double* B, double* C,int n){
  int i,j,k;
  double t;
  for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        t=C[i+j*n];
        for(k=0;k<n;k++){
          t+=A[i+k*n]*B[k+j*n];        
        }
        C[i+j*n]=t;
      }
    }
}

void strassen(int n, double* A, double* B, double* C){
  //base
  if(n==1){
    C[0]=A[0]*B[0];
    return;
  }
  //spliting
  int ns2=n/2,i,j;
  int ns22=ns2*ns2;
  int sizeofns22=ns22*sizeof(double);
  double* A11=(double*)malloc(sizeofns22),*A12=(double*)malloc(sizeofns22),*A21=(double*)malloc(sizeofns22),*A22=(double*)malloc(sizeofns22);
  double* B11=(double*)malloc(sizeofns22),*B12=(double*)malloc(sizeofns22),*B21=(double*)malloc(sizeofns22),*B22=(double*)malloc(sizeofns22);
  double* C11=(double*)malloc(sizeofns22),*C12=(double*)malloc(sizeofns22),*C21=(double*)malloc(sizeofns22),*C22=(double*)malloc(sizeofns22);
  for(i=0;i<ns2;i++)
    for(j=0;j<ns2;j++){
      A11[i+j*ns2]=A[(i)+(j)*n];
      A12[i+j*ns2]=A[(i)+(j+ns2)*n];
      A21[i+j*ns2]=A[(i+ns2)+(j)*n];
      A22[i+j*ns2]=A[(i+ns2)+(j+ns2)*n];
      B11[i+j*ns2]=B[(i)+(j)*n];
      B12[i+j*ns2]=B[(i)+(j+ns2)*n];
      B21[i+j*ns2]=B[(i+ns2)+(j)*n];
      B22[i+j*ns2]=B[(i+ns2)+(j+ns2)*n];
    }
  double* M1=(double*)malloc(sizeofns22),*M2=(double*)malloc(sizeofns22),*M3=(double*)malloc(sizeofns22),*M4=(double*)malloc(sizeofns22),*M5=(double*)malloc(sizeofns22),*M6=(double*)malloc(sizeofns22),*M7=(double*)malloc(sizeofns22);
  double* T1=(double*)malloc(sizeofns22),*T2=(double*)malloc(sizeofns22);
  //coefficients
  if(ns2<=cacheDouble){
    sum(A11,A22,T1,ns2,1);
    sum(B11,B22,T2,ns2,1);
    implicito(ns2,T1,T2,M1);
    sum(A21,A22,T1,ns2,1);
    implicito(ns2,T1,B11,M2);
    sum(B12,B22,T1,ns2,-1);
    implicito(ns2,A11,T1,M3);
    sum(B21,B11,T1,ns2,-1);
    implicito(ns2,A22,T1,M4);
    sum(A11,A12,T1,ns2,1);
    implicito(ns2,T1,B22,M5);
    sum(A21,A11,T1,ns2,-1);
    sum(B11,B12,T2,ns2,1);
    implicito(ns2,T1,T2,M6);
    sum(A12,A22,T1,ns2,-1);
    sum(B21,B22,T2,ns2,1);
    implicito(ns2,T1,T2,M7);

  }else{
    sum(A11,A22,T1,ns2,1);
    sum(B11,B22,T2,ns2,1);
    strassen(ns2,T1,T2,M1);
    sum(A21,A22,T1,ns2,1);
    strassen(ns2,T1,B11,M2);
    sum(B12,B22,T1,ns2,-1);
    strassen(ns2,A11,T1,M3);
    sum(B21,B11,T1,ns2,-1);
    strassen(ns2,A22,T1,M4);
    sum(A11,A12,T1,ns2,1);
    strassen(ns2,T1,B22,M5);
    sum(A21,A11,T1,ns2,-1);
    sum(B11,B12,T2,ns2,1);
    strassen(ns2,T1,T2,M6);
    sum(A12,A22,T1,ns2,-1);
    sum(B21,B22,T2,ns2,1);
    strassen(ns2,T1,T2,M7);
  }
  //results
  sum(M1,M4,T1,ns2,1);
  sum(M7,M5,T2,ns2,-1);
  sum(T1,T2,C11,ns2,1);
  sum(M3,M5,C12,ns2,1);
  sum(M2,M4,C21,ns2,1);
  sum(M1,M2,T1,ns2,-1);
  sum(M3,M6,T2,ns2,1);
  sum(T1,T2,C22,ns2,1);
  //compiling
  for(i=0;i<ns2;i++)
    for(j=0;j<ns2;j++){
      C[i+j*n]=C11[i+j*ns2];
      C[(i)+(j+ns2)*n]=C12[i+j*ns2];
      C[(i+ns2)+(j)*n]=C21[i+j*ns2];
      C[(i+ns2)+(j+ns2)*n]=C22[i+j*ns2];
    }

}



void square_dgemm(int n, double* A, double* B, double* C){
  int cacheLine = 64; //Linha de cache em Kb
  cacheDouble = (cacheLine*1000)/sizeof(double); //Quantos doubles cabem na linha
  int newSize=pow(2,ceil(log(n)/log(2))),i,j;
  double AS[newSize*newSize], BS[newSize*newSize], CS[newSize*newSize];
  fillWithZeros(A,AS,n,newSize);
  fillWithZeros(B,BS,n,newSize);
  fillWithZeros(C,CS,n,newSize);
  strassen(newSize,AS,BS,CS);
  for(i=0;i<n;i++)
  for(j=0;j<n;j++)
    C[i+j*n]=CS[i+j*newSize];
}
//------------------------------------------------------------//
*/
