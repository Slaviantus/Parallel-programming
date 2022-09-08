// System od dif equations solver with openMP.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include<iostream>
#include"functionsParallel.h"
#include"conio.h"
#include<omp.h>
#include<stdio.h>
#include<stdlib.h>

usingnamespacestd;

doubleF_v(intnum, double* g, doubleV, double* r, doublem, doubleh, doublen, doubleI_ext)
{
	doubleI_k = g_k * n * n * n * n * (V + V_k);
	doubleI_na = g_na * m * m * m * h * (V + V_na);
	doubleI_leak = g_l * (V + V_l);
	doubleI_syn = 0.0;

	for (int j = 0; j < N; j++)
		I_syn = I_syn - r[j] * g[j] * (E_syn - V);

	double res = -(I_k + I_na + I_leak + I_syn + I_ext);
	return res;
}
doubleF_m(doubleV, doublem)
{
	double alpha = (0.1 * (-V + 25.0)) / (exp((-V + 25.0) / 10.0) - 1.0);
	double beta = 4 * exp(-V / 18.0);
	double res = alpha * (1.0 - m) - beta * m;
	return res;
}
doubleF_h(doubleV, doubleh)
{
	double alpha = 0.07 * exp(-V / 20.0);
	double beta = 1 / (exp((-V + 30.0) / 10.0) + 1.0);
	double res = alpha * (1 - h) - beta * h;
	return res;
}
doubleF_n(doubleV, doublen)
{
	double alpha = (0.01 * (-V + 10.0)) / (exp((-V + 10.0) / 10.0) - 1.0);
	double beta = 0.125 * exp(-V / 80.0);
	double res = alpha * (1 - n) - beta * n;
	return res;
}
doubleF_r(doubleV, doubler)
{
	double res;
	double T = 1 / (1 + exp(-(V - 20.0) / 2.0));
	res = a * T * (1 - r) - b * r;
	return res;
}
doubleF_gi(doubler_pre, doubler_post, doubleg)
{
	return 0.0;
}



constint N = 1000;
constdoubledt = 0.01;
double* G = newdouble[N * N];
doubleVold[N], Vnew[N];
double m[N], h[N], n[N];
doubleI_ext[N];
doublerold[N], rnew[N];


constint seed = 1;


voidGenIExt()
{
	inti;
	for (i = 0; i < N; i++)
		I_ext[i] = -10.0 - (double)rand() / ((double)RAND_MAX);
}



voidInit()
{
	inti, j;
	for (i = 0; i < N; i++)
	{
		rold[i] = m[i] = n[i] = h[i] = 0.0;

		for (j = 0; j < N; j++)
			G[i * N + j] = 0.001;

		Vold[i] = -9.0;
	}
	srand(seed);
}




void Solve(intTime, intnumOfThreads)
{
	inti, count;
	#pragmacomment(linker, "/STACK:100000000")
		omp_set_num_threads(numOfThreads); // set the number of using threads

	// Циклповремени
	for (count = 0; count < Time; count++)
	{

		GenIExt();
		#pragmaomp parallel
		{
			// Циклпонейронам
#pragmaompforschedule(guided, 10)
				for (i = 0; i < N; i++)
				{
					Vnew[i] = Vold[i] + F_v(i, &G[i * N], Vold[i], rold, m[i], h[i], n[i], I_ext[i]) * dt;
					m[i] = m[i] + F_m(Vold[i], m[i]) * dt;
					h[i] = h[i] + F_h(Vold[i], h[i]) * dt;
					n[i] = n[i] + F_n(Vold[i], n[i]) * dt;
					rnew[i] = rold[i] + F_r(Vold[i], rold[i]) * dt;
				}
		}

			// Копирование
			#pragmaomp parallel
		{
#pragmaompforschedule(guided, 10)
				for (i = 0; i < N; i++)
				{
					Vold[i] = Vnew[i];
					rold[i] = rnew[i];
				}
		}

	}

}

void Print()
{
	for (inti = 0; i < N; i++)
	{

		cout << Vold[i] << " " << m[i] << " " << h[i] << " " << n[i] << " " << rold[i] << endl;
	}

}




int main(intargc, char* argv[])
{
	intnumOfThreads = 1;
	cout << "Set the number of threads" << endl;
	cout << "________________________________________" << endl;
	cin >> numOfThreads;



	clock_t t1, t2;

	Init();

	t1 = clock();

	Solve(2000, numOfThreads);

	t2 = clock();

	Print();

	cout << "________________________________________" << endl;
	cout << "time = " << double(t2 - t1) / CLOCKS_PER_SEC << " seconds" << endl;

	_getch();
	return 0;
}
