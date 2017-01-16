//******************************************************
// Project
// Names: Anthony Enem
// Parallel Programming Date: 12/05/16
//******************************************************
// This program implements the cooley tukey fft algorithm
// and computes the values fro X_k from 0 to N pin parallel. 
//******************************************************

#include<cuda.h>
#include<sys/time.h>
#include<stdio.h>
#include<math.h>
#include<time.h>

//define struct for complex numbers
typedef struct complex{
    double real;
    double imaginary;
} complex;

#define N 16384
#define PI 3.14159265359

//******************************************************
// Parameters: two complex structs a and b
// This function multiplies 2 complex structs and returns
// the result as a complex struct
//******************************************************

__device__
complex multiplyComplex(complex a, complex b)
{
    complex result;
    result.real = a.real*b.real - a.imaginary*b.imaginary;
    result.imaginary = a.real*b.imaginary + b.real*a.imaginary;

    return result;
}

//******************************************************
// Parameters: two complex structs a and b
// This function adds 2 complex structs and returns
// the result as a complex struct
//******************************************************
__device__
complex addComplex(complex a, complex b)
{
    complex result = {a.real+b.real, a.imaginary+b.imaginary};
    return result;
}

//******************************************************
// Parameters: two complex structs a and b
// This function subtracts 2 complex structs and returns
// the result as a complex struct
//******************************************************
__device__
complex subtractComplex(complex a, complex b)
{
    complex result = {a.real-b.real, a.imaginary-b.imaginary};
    return result;
}

//******************************************************
// Each thread for this kernel computes the even and odd portions of the
// fft formula
//******************************************************
__global__
void compute_Even_Odd(complex* xInput_d, complex *even_d, complex* odd_d)
{	
	int idx = blockIdx.x*1024 + threadIdx.x;
	int k = idx / (N/2);
	int m = idx % (N/2);
	
	//compute exponent
    double exp = -4.0 * PI * m * k / N;
	
	complex e_part;    

    //compute part with e^(i * exponent)
    e_part.real = cos(exp);
    e_part.imaginary = sin(exp);
        
    //save even part
	even_d[idx] = multiplyComplex(xInput_d[2*m], e_part);
	
    //save odd part
    odd_d[idx] = multiplyComplex(xInput_d[2*m+1], e_part);
}

//******************************************************
// Parameters: pointer to complex twiddle and integer k
// This function computes the twiddle factor for given
// N and k values
//******************************************************
__device__
void compute_twiddle(complex* twiddle, int* k)
{
    double exp = -2.0 * PI * (*k) / N;
	twiddle->real = cos(exp);
	twiddle->imaginary = sin(exp);
}

//******************************************************
// This computes the summation for the fft algorithm 
// using regression
//******************************************************
__global__
void sum_Regression(complex* even_d, complex* odd_d, complex* result_d)
{	
	int n = N/2;
	int idx = blockIdx.x*1024 + threadIdx.x;
	int k = idx / (n/2); //we only need half the number of threads for regression
	int m = idx % (n/2);
	int z = k*n + m;
	
	while(n > 1)
	{
		__syncthreads();		
		if(m < n/2)
		{
			even_d[z] = addComplex(even_d[z], even_d[z + n/2]);
			odd_d[z] = addComplex(odd_d[z], odd_d[z + n/2]);
		}		
		n = n/2;
	}
	if(m == 0)
	{
		compute_twiddle(result_d + k, &k);
		//calculate twiddle * O_k
		odd_d[z] = multiplyComplex(odd_d[z], result_d[k]);
		
		//compute X_k and X_k + N/2
		result_d[k] = addComplex(even_d[z], odd_d[z]);
		result_d[k + N/2] = subtractComplex(even_d[z], odd_d[z]);
	}
}

int main()
{	
	complex* xInput, *result;
	
	complex* xInput_d, *result_d;
	complex *even_d, *odd_d;

	//allocate memory
	xInput = (complex *)malloc(N * sizeof(complex));
	result = (complex *)malloc(N * sizeof(complex));
	
	cudaMalloc((void**)&xInput_d, N*sizeof(complex));
	cudaMalloc((void**)&result_d, N*sizeof(complex));
	
	//set defined values from 0 to 7 for xInput
	xInput[0].real = 3.6;
	xInput[0].imaginary = 2.6;
	xInput[1].real = 2.9;
	xInput[1].imaginary = 6.3;
	xInput[2].real = 5.6;
	xInput[2].imaginary = 4;
	xInput[3].real = 4.8;
	xInput[3].imaginary = 9.1;
	xInput[4].real = 3.3;
	xInput[4].imaginary = 0.4;
	xInput[5].real = 5.9;
	xInput[5].imaginary = 4.8;
	xInput[6].real = 5;
	xInput[6].imaginary = 2.6;
	xInput[7].real = 4.3;
	xInput[7].imaginary = 4.1;
	
	//set defined values for x[8] to x[N-1] as 0
	for (int i = 8; i < N; i++)
	{
		xInput[i].real = 0;
		xInput[i].imaginary = 0;
	}
	
	clock_t st = clock();
	cudaMemcpy(xInput_d, xInput, sizeof(complex)*N, cudaMemcpyHostToDevice);
	
	int m_blocks, m_threads;
	m_threads = (N/2)*(N/2); //k goes 0 to N - 1 (but same values used 2ice),
							//m goes 0 to N/2-1 for each k value
	
	cudaMalloc((void**)&even_d, sizeof(complex)*m_threads);
	cudaMalloc((void**)&odd_d, sizeof(complex)*m_threads);
	
	m_blocks = (m_threads)/1024; //number of blocks we need
	
	if(m_blocks == 0){
		m_blocks = 1;
	}
	else{
		m_threads = 1024;
	}
	
	//Use half the number of blocks to compute results using regression
	dim3 dimSumGrid(m_blocks/2, 1);
	dim3 dimSumBlock(m_threads, 1);
	if(m_blocks == 1)
	{
		dimSumGrid = dim3(m_blocks, 1);
		dimSumBlock = dim3(m_threads/2, 1);
	}
		
	
	cudaEvent_t start, stop;
	float elapsedTime;

	cudaEventCreate(&start);
	cudaEventRecord(start,0);
	
	compute_Even_Odd<<<m_blocks, m_threads>>>(xInput_d, even_d, odd_d);
	
	sum_Regression<<<dimSumGrid, dimSumBlock>>>(even_d, odd_d, result_d);
	
	clock_t end = clock();
	cudaEventCreate(&stop);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start,stop);
	
	printf("Cuda time : %f milli secs\n", elapsedTime);	
	
	cudaMemcpy(result, result_d, sizeof(complex)*N, cudaMemcpyDeviceToHost);
	
	cudaFree(xInput_d);
	cudaFree(result_d);
	cudaFree(even_d);
	cudaFree(odd_d);
	
	printf("PARALLEL VERSION\n\n");
	printf("TOTAL PROCESSED SAMPLES: %d\n", N);
	
	for(int i = 0; i < 8; i++)
	{
		printf("==================\n");
		printf("XR[%d]: %f\n", i, result[i].real);
		printf("XI[%d]: %f\n", i, result[i].imaginary);
	}
	
	printf("==================\n");
	
	//printf("\nC clock Time: %f\n", ((double)(end - st)) / CLOCKS_PER_SEC);
	
	
	return 0;
	
}
