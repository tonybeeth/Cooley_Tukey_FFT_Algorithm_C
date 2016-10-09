#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define PI acos(-1)

//define struct for complex numbers
typedef struct complex{
    double real;
    double imaginary;
} complex;

//create constant complex array for input
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
    complex result;
    result.real = a.real*b.real - a.imaginary*b.imaginary;
    result.imaginary = a.real*b.imaginary + b.real*a.imaginary;

    return result;
}

//Function to add  two complex numbers
complex addComplex(complex a, complex b)
{
    complex result = {a.real+b.real, a.imaginary+b.imaginary};
    return result;
}

//Function to solve even and odd portions of formula
void compute_Even_Odd(int N, int k, complex* E_k, complex* O_k)
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
        //printf("m: %d", m);

        //compute exponent
        exp = -2.0 * PI * m * k / N_half;
        //printf("%f\n", exp);

        //compute part with e^(i * exponent)
        e_part.real = cos(exp);
        e_part.imaginary = sin(exp);
        //printf("e_part: %f, %f\n", e_part.real, e_part.imaginary);
        
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

        //printf("Even: %f, %f\n",E_k->real, E_k->imaginary);
        //printf("Odd: %f, %f\n\n", O_k->real, O_k->imaginary);
    }

}

//function to compute twiddle factor
complex compute_twiddle(int N, int k)
{
    double exp = -2.0 * PI * k / N;
    complex twiddle = { cos(exp), sin(exp) };

    return twiddle;        
}


int main(int argc, char** argv)
{
    if(argc != 2)
    {
        printf("Wrong usage!\n");
        printf("Must specify value for N!\n");

        return 1;
    }
    //get N from command line argument
    int N = atoi(argv[1]);

    complex result;
    complex twiddle_factor;

    //allocate memory for E_k and O_k
    complex* E_k = malloc(sizeof(complex));
    complex* O_k = malloc(sizeof(complex));
    //printf("Even: %f", solve_Even_Odd(8, 0).imaginary);    

    for(int k = 0; k < 8; k++)
    {
        //get twiddle factor
        twiddle_factor = compute_twiddle(N, k);

        //compute even and odd portions
        compute_Even_Odd(N, k, E_k, O_k);

        //add up results
        result = addComplex(*E_k, multiplyComplex(twiddle_factor, *O_k));
        printf("X[%d]: %f, %f\n\n", k, result.real, result.imaginary);
    }

    return 0;

}