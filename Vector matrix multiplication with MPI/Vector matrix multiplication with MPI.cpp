// Vector matrix multiplication with MPI.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include<stdio.h>
#include<mpi.h>
#include<fstream>
#include<iostream>
#include<cmath>
#include<conio.h>

usingnamespacestd;


/*____________ Reading file ____________*/

intbinRead(doubleoutputArray[], conststring& fileName, constunsignedintlength)
{
	ifstreaminputData;
	inputData.open(fileName);
	if (inputData)
	{
		inputData.read(reinterpret_cast<char*>(outputArray), sizeof(double) * length);
		inputData.close();
	}
	else
	{
		return -1;
	}

}

/*____________ Writing file ____________*/

intbinWrite(constdoubleinputArray[], conststring& fileName, constunsignedintlength)
{
	ofstreamoutputData;
	outputData.open(fileName);
	if (outputData)
	{
		outputData.write(reinterpret_cast<constchar*>(inputArray), sizeof(double) * length);
		outputData.close();
		return 0;
	}
	else
	{
		return -1;
	}
}


/*____________ Showing Array ____________*/

voidshowMatrix(double* matrix, constintn)
{
	cout << "===============  Array A  ===============" << endl;

	for (inti = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << matrix[i * n + j] << " ";
		}

		cout << endl;
	}
}


/*____________ Showing Vector ____________*/

voidshowVector(double* matrix, constintn)
{
	for (inti = 0; i < n; i++)
	{
		cout << matrix[i] << endl;
	}
}


/*____________ Filling Array ____________*/

voidfillMatrix(double* matrix, constintn)
{
	cout << "===============  Fill the matrix A  ===============" << endl;

	for (inti = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "Fill the element of " << i + 1 << " row and " << j + 1 << " column" << endl;
			cin >> matrix[i * n + j];
		}
	}

}


/*____________ FillingVector ____________*/

voidfillVector(double* matrix, constintn)
{
	cout << "===============  Fill the matrix B  ===============" << endl;

	for (inti = 0; i < n; i++)
	{
		cout << "Fill the element of " << i + 1 << " column" << endl;
		cin >> matrix[i];
	}

}



int main(intargs, char* argv[])
{
	int n = 3;
	double A[9];
	double Result[9];
	double Result1[9];
	double B[3];
	double C[3];
	double p[3];
	doubleendResult[3];
	doubleStartTimer, FinishTimer;
	intmyrank, nprocs, j, i, k;

	/*____________ Filling Arrays ____________*/
	MPI_Init(&args, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (myrank == 0)
	{
		cout << "Amount of processes nprocs = " << nprocs << endl;
		fillMatrix(A, n);
		fillVector(B, n);
		for (inti = 0; i < n; i++)
		{
			C[i] = 0;
			p[i] = 0;
		}
		for (inti = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				Result[i * n + j] = 0;
				Result1[i * n + j] = 0;
			}
		}
	}

	/*____________ WritingArraysintofile ____________*/
	if (myrank == 0)
	{
		interrorA;
		interrorB;
		errorA = binWrite(A, "MatrixA.bin", 9);
		errorB = binWrite(B, "MatrixB.bin", 3);

		/*____________ Reading Arrays from file ____________*/

		if ((errorA == 0) && (errorB == 0))
		{
			errorA = binRead(A, "MatrixA.bin", 9);
			errorB = binRead(B, "MatrixB.bin", 3);
		}
	}
	/*____________ ShowingArrays ____________*/

	if (myrank == 0)
	{
		showMatrix(A, n);
		cout << "===============  Array B  ===============" << endl;
		showVector(B, n);
	}

	MPI_Bcast(B, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(A, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(C, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(p, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(Result, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(Result1, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&StartTimer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&FinishTimer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/*____________ ParallelMultiplication ____________*/

	StartTimer = MPI_Wtime();

	for (j = myrank; j < n; j += nprocs)
	{
		for (inti = 0; i < n; i++)
		{
			Result[i * n + j] = A[i * n + j] * B[j];
		}
	}

	MPI_Reduce(Result, Result1, 9, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Bcast(Result1, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (k = myrank; k < n; k += nprocs)
	{
		for (int l = 0; l < n; l++)
		{
			p[k] += Result1[k * n + l];
		}
		C[k] = p[k];
	}
	MPI_Reduce(C, endResult, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	FinishTimer = MPI_Wtime();

	/*____________ Writing Result into file ____________*/
	if (myrank == 0)
	{
		interrorC;
		errorC = binWrite(C, "Result.bin", 3);

		/*____________ Reading Result from file ____________*/

		if (errorC == 0)
		{
			errorC = binRead(A, "Result.bin", 3);
		}
	}
	/*____________ Showing Result ____________*/

	if (myrank == 0)
	{
		cout << "===============  Result Vector  ===============" << endl;
		showVector(endResult, n);
		cout << "_______________ Time of processes _______________" << endl;
	}

	cout << abs(FinishTimer - StartTimer) << " Seconds" << endl;

	MPI_Finalize();

	return 0;
}
