#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define PI acos(-1)

typedef struct complex{
    double real;
    double imaginary;
} complex;

const complex x_input[8] = 
    {
        {3.6, 2.6}, 
        {2.9, 6.3}, 
        {5.6, 4}, 
        {4.8, 9.1},
        {3.3, 0.4}, 
        {5.9, 4.8},
        {5, 2.6}, 
        {4.3, 4.1}
    };

//Function to multiply two complex numbers
complex multiplyComplex(complex a, complex b)
{

}

//Function to add complex numbers
complex addComplex(complex a, complex b)
{

}

//Function to solve E_k (even portion) of X_k
complex solve_Even(int N, int k)
{
    if(k >= N/2){
        k -= N/2;
    }

    int N_half = N/2;
    complex E_k = {0, 0};
    complex e_part;
    double exp;

    for(int m = 0; m <= N/2 - 1; m++)
    {
        exp = -2.0 * PI * m * k / N_half;
        //printf("%f\n", exp);
        e_part.real = cos(exp);
        e_part.imaginary = sin(exp);

        E_k = addComplex(
            E_k,
            multiplyComplex(x_input[2*m], e_part)
        );
        //printf("%f, %f\n", E_k.real, E_k.imaginary);
    }

    return E_k;
}

void solve_Xk(int N, int k)
{
    complex twiddle;
    complex X_k;


    
}


int main(int argc, char** argv)
{
    if(argc != 2)
    {
        printf("Wrong usage!\n");
        printf("Must specify value for N!\n");

        return 1;
    }
    int N = atoi(argv[1]);

    printf("N: %d", N);    

    return 0;

}