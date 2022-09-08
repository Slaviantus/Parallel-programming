// Matrix multiplication with openMP.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include<iostream>
#include<omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<ctime>

#defineN2000  // set the size of matrix 
usingnamespacestd;
/*____________ Showing Array ____________*/
voidshowMatrix(double* matrix, constintn)
{
	for (inti = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << matrix[i * n + j] << " ";
		}

		cout << endl;
	}
}
/*____________ Filling Array ____________*/
voidfillMatrix(double* matrix, constintn)
{
	for (inti = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "Fill the element of " << i + 1 << " row and " << j + 1 << " column" << endl;
			cin >> matrix[i * n + j];
		}
	}
}
/*____________ Creating Random Array ____________*/
voidrandomMatrix(double* matrix, constintn)
{
	int minimum = -10;
	int range = 20;


	for (inti = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			matrix[i * n + j] = minimum + rand() % range;
		}
	}
}
int main()
{
	#pragmacomment(linker, "/STACK:100000000")
		inti;
	doubleStartTimer;
	doubleFinishTimer;
	doublewtime;
	double A[N * N];
	double B[N * N];
	double C[N * N];
	/*____________ Initialisation C ____________*/
	for (inti = 0; i < N * N; i++)
	{
		C[i] = 0;
	}
	/*____________ filling random Arrays ____________*/
	randomMatrix(A, N);
	srand(time(NULL));
	randomMatrix(B, N);
	/*____________ Showing Arrays ____________*/

	cout << "===============  Array A  ===============" << endl;
	//showMatrix(A, N);
	cout << "===============  Array B  ===============" << endl;
	//showMatrix(B, N);

	/*____________ OpenMp Multiplication ____________*/

	cout << "===============  Threads  ===============" << endl;
	omp_set_num_threads(4); // set the number of Processes
	#pragmaomp parallel
	{
		StartTimer = omp_get_wtime();

	#pragmaompforschedule(guided, 10)
		for (i = 0; i < N; i++)
		{

			for (int k = 0; k < N; k++)
			{
				for (int j = 0; j < N; j++)
				{
					C[N * i + k] += A[N * i + j] * B[N * j + k];
				}
			}

		}

		FinishTimer = omp_get_wtime();
	}
	wtime = abs(FinishTimer - StartTimer);


	/*____________ Output result ____________*/

	cout << "===============  Result of multiplication  ===============" << endl;
	showMatrix(C, N);
	cout << "===============  Max time of threads  ===============" << endl;
	cout << wtime * 1000 << " milliseconds" << endl;

	getchar();
	return 0;
}
