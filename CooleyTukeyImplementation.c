//Assigntment #1
//Name: Cavaughn Browne and Anthony Enem
//Parallel Programming Date: 10/9/2016

//Reads N sets of data from a file called data.dat and processes them using 
//the FFT-Cooley Tukey Algorithm

#include <stdio.h>
#include <math.h>

#define N 8
#define MAX_K 8
#define PI 3.141592653589793

struct ComplexNum 
{
	double real;
	double imag;
};

int main()
{
	
	FILE *fp;

	fopen_s(&fp, "data.dat", "r");
	
	struct ComplexNum Xk;
	struct ComplexNum x[N];
	struct ComplexNum evenP;
	struct ComplexNum oddP;
	double c, s, realPart, imgPart;

	for (int i = 0; i < N; i++)
	{
		fscanf_s(fp, "%lf", &x[i].real);
		fscanf_s(fp, "%lf", &x[i].imag);

	}

	printf("TOTAL PROCESSED SAMPLES: %d\n", N);

	for (int k = 0; k < MAX_K; k++)
	{
		double theta = (-2 * PI * k) / (N / 2);

		evenP.real = 0;
		evenP.imag = 0;
		oddP.real = 0;
		oddP.imag = 0;

		for (int m = 0; m < N / 2; m++)
		{
			//Even
			c = cos(theta * m);
			s = sin(theta * m);
			realPart = (x[2 * m].real *c) - ((x[2 * m].imag * s));
			evenP.real += realPart;
			imgPart = (x[2 * m].real *s) + ((x[2 * m].imag * c));
			evenP.imag += imgPart;

			//Odd
			realPart = (x[(2 * m) + 1].real *c) - ((x[(2 * m) + 1].imag * s));
			oddP.real += realPart;
			imgPart = (x[(2 * m) + 1].real *s) + ((x[(2 * m) + 1].imag * c));
			oddP.imag += imgPart;
		}

		Xk.real = evenP.real + (cos(theta / 2) * oddP.real) - (sin(theta / 2) * oddP.imag);
		Xk.imag = evenP.imag + (cos(theta / 2) * oddP.imag) + (sin(theta / 2) * oddP.real);

		printf("XR[%d] : %lf\nXI[%d] : %lf\n", k, Xk.real, k, Xk.imag);
	}
}