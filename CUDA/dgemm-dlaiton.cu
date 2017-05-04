const char* dgemm_desc = "Dlaiton, best way dgemm.";

#define THREADSPERBLOCK 256


__global__ 
void dlaitonComeback(int n,float *A,float *B,float *C){
	int tid=threadIdx.x+(blockDim.x*blockIdx.x);
	if(tid>=n*n) return;
	int i=tid%n;
	int j=tid/n;
	for(int k=0;k<n;k++){
        C[i+j*n]+=A[i+k*n]*B[k+j*n];        
    }
}


__host__
void square_dgemm(int n, double* A, double* B, double* C){
	float *d_A,*d_B,*d_C; 
	int n2=n*n;
	int size=n*sizeof(float);
	cudaMalloc((void**)&d_A,size);
    cudaMalloc((void**)&d_B,size);
    cudaMalloc((void**)&d_C,size);
    cudaMemcpy(d_A,A,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_B,B,size,cudaMemcpyHostToDevice);
	cudaMemcpy(d_C,C,size,cudaMemcpyHostToDevice);
	dlaitonComeback<<<(n2+THREADSPERBLOCK-1)/THREADSPERBLOCK,THREADSPERBLOCK>>>(n,d_A,d_B,d_C);
	cudaMemcpy(C,d_C,size,cudaMemcpyDeviceToHost);
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);
}