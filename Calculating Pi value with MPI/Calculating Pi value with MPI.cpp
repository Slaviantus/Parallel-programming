// Calculating Pi value with MPI.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <stdio.h>
#include <mpi.h>
#include <iostream>
#include <cmath>

using namespace std;

double f(double x)
{
	return 1 / (1 + x * x);
}


int main(int args, char* argv[])
{

	/*____________ Average Rectangle Method ____________*/

	long double pi, sum = 0, term, h;
	const double realPi = 3.141592653589793238462643;
	int myrank, nprocs, n, i;
	double t1, t2, dt;

	MPI_Init(&args, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if (myrank == 0)
	{
		cout << "The amout of Process nprocs = " << nprocs << endl;
		cout << "Enter the number of iterations = " << endl;
		scanf_s("%d", &n);
	}

	t1 = MPI_Wtime();

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	h = 1.0 / n;

	for (i = myrank + 1; i <= n; i += nprocs)
	{
		sum += f(h * (i - 0.5));
	}

	term = 4 * h * sum;

	MPI_Reduce(&term, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	t2 = MPI_Wtime();
	dt = t2 - t1;

	if (myrank == 0)
	{
		cout << "============  Average Rectangle Method  ============" << endl;
		cout << "Computed value of pi = " << pi << "   Difference = " << abs(pi - realPi) << "   Time = " << dt << " Seconds" << endl;
	}

	/*____________ Trapezium Method ____________*/

	sum = 0;
	double a = 0;
	double b = 1;
	t1 = MPI_Wtime();
	double internalTerm = (f(a) + f(b)) / 2;
	for (i = myrank + 1; i <= n - 1; i += nprocs)
	{
		sum += f(a + i * h);
	}
	MPI_Reduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	pi = 4 * h * (internalTerm + pi);

	t2 = MPI_Wtime();
	dt = t2 - t1;
	if (myrank == 0)
	{
		cout << "============  Trapezium Method  ============" << endl;
		cout << "Computed value of pi = " << pi << "   Difference = " << abs(pi - realPi) << "   Time = " << dt << " Seconds" << endl;
	}
	/*____________ Simpson Method ____________*/
	sum = 0;
	double secondSum = 0;
	double result1 = 0;
	double result2 = 0;
	t1 = MPI_Wtime();
	internalTerm = (f(a) + f(b)) / 2;

	for (i = myrank + 1; i <= n - 1; i += nprocs)
	{
		sum += f(a + i * h);
	}
	MPI_Reduce(&sum, &result1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	for (i = myrank + 1; i <= n; i += nprocs)
	{
		secondSum += f(a + (i - 0.5) * h);
	}
	MPI_Reduce(&secondSum, &result2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	pi = 4 * h * (internalTerm + 2 * result2 + result1) / 3;

	t2 = MPI_Wtime();
	dt = t2 - t1;

	if (myrank == 0)
	{
		cout << "============  Simpson Method  ============" << endl;
		cout << "Computed value of pi = " << pi << "   Difference = " << abs(pi - realPi) << "   Time = " << dt << " Seconds" << endl;
	}
	MPI_Finalize();
	return 0;
}
