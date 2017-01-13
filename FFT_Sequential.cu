//******************************************************
// Assignment #1
// Names: Anthony Enem and Cavaughn Browne
// Parallel Programming Date: 10/10/16
//******************************************************
// This program implements the cooley tukey fft algorithm
// and computes the values fro X_k from 0 to N. The program
// should be compiled using "gcc -o cooley_fft Browne_Enem_A1.c"
// Then it should be run with "./cooley_fft N" where N is the
// maximum value N in the fft formula
//******************************************************

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

double PI = acos(-1);
#define N 16384


//define struct for complex numbers
typedef struct complex{
    double real;
    double imaginary;
} complex;

//create complex pointer for input
complex* x_input = NULL;

//******************************************************
// Parameters: two complex structs a and b
// This function multiplies 2 complex structs and returns
// the result as a complex struct
//******************************************************
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
complex subtractComplex(complex a, complex b)
{
    complex result = {a.real-b.real, a.imaginary-b.imaginary};
    return result;
}

//******************************************************
// Parameters: integer N, integer k, pointers to complex 
//		structs E_k and O_k
// This function computes the even and odd portions of the
// fft formula and returns them by reference in the complex
// pointers
//******************************************************
void compute_Even_Odd(int k, complex* E_k, complex* O_k)
{
    if(k >= N/2){
        k -= N/2;
    }

    E_k->real = 0;
    E_k->imaginary = 0;
    O_k->real = 0;
    O_k->imaginary = 0;

    int N_half = N/2;
    complex e_part;
    double exp;

    for(int m = 0; m <= N/2 - 1; m++)
    {
        //compute exponent
        exp = -2.0 * PI * m * k / N_half;

        //compute part with e^(i * exponent)
        e_part.real = cos(exp);
        e_part.imaginary = sin(exp);
        
        //add even part
        *E_k = addComplex(
            *E_k,
            multiplyComplex(x_input[2*m], e_part)
        );

        //add odd part
        *O_k = addComplex(
            *O_k,
            multiplyComplex(x_input[2*m+1], e_part)
        );

    }

}

//******************************************************
// Parameters: integer N and integer k
// This function computes the twiddle factor for given
// N and k values
//******************************************************
complex compute_twiddle(int k)
{
    double exp = -2.0 * PI * k / N;
    complex twiddle = { cos(exp), sin(exp) };

    return twiddle;        
}


//Main function
int main(int argc, char** argv)
{
	//make sure user enters a value for N
	/*if (argc != 2)
	{
		printf("Wrong usage!\n");
		printf("Must specify value for N!\n");

		return 1;
	}*/
	//get N from command line argument

	//allocate memory for N complex inputs
	x_input = (complex *)malloc(N * sizeof(complex));
	complex* result = (complex*)malloc(N *sizeof(complex));
	
	//set defined values from 0 to 7 for x_input
	x_input[0].real = 3.6;
	x_input[0].imaginary = 2.6;
	x_input[1].real = 2.9;
	x_input[1].imaginary = 6.3;
	x_input[2].real = 5.6;
	x_input[2].imaginary = 4;
	x_input[3].real = 4.8;
	x_input[3].imaginary = 9.1;
	x_input[4].real = 3.3;
	x_input[4].imaginary = 0.4;
	x_input[5].real = 5.9;
	x_input[5].imaginary = 4.8;
	x_input[6].real = 5;
	x_input[6].imaginary = 2.6;
	x_input[7].real = 4.3;
	x_input[7].imaginary = 4.1;

	//set defined values for x[8] to x[N-1] as 0
	for (int i = 8; i < N; i++)
	{
		x_input[i].real = 0;
		x_input[i].imaginary = 0;
	}

    complex twiddle_factor;

    //allocate memory for E_k and O_k
    complex* E_k = (complex *)malloc(sizeof(complex));
    complex* O_k = (complex *)malloc(sizeof(complex));  	
	
	clock_t st = clock();
	
	cudaEvent_t start, stop;
	float elapsedTime;

	cudaEventCreate(&start);
	cudaEventRecord(start,0);
	

    for(int k = 0; k < N/2; k++)
    {
        //get twiddle factor
        twiddle_factor = compute_twiddle(k);

        //compute even and odd portions
        compute_Even_Odd(k, E_k, O_k);

        //add up results
        result[k] = addComplex(*E_k, multiplyComplex(twiddle_factor, *O_k));
		result[k + N/2] = subtractComplex(*E_k, multiplyComplex(twiddle_factor, *O_k));
    }
	
	clock_t end = clock();
	cudaEventCreate(&stop);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&elapsedTime, start,stop);
	
	//PRINT FIRST 8 RESULTS
	
	printf("SEQUENTIAL VERSION\n\n");
	
	printf("TOTAL PROCESSED SAMPLES: %d\n", N);	
	for(int i = 0; i < 8; i++)
	{	
		printf("==================\n");
		printf("XR[%d]: %f\n", i, result[i].real);
		printf("XI[%d]: %f\n", i, result[i].imaginary);
	}

	
	printf("==================\n");
	
	printf("Cuda time : %f secs\n", elapsedTime/1000);
	
	printf("\nC clock Time: %f\n", (double)(end - st) / CLOCKS_PER_SEC);

    return 0;

}