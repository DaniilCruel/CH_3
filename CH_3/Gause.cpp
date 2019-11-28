#include "head.h"


bool fixz(double** &A,int u,int n)
{
	double max = A[u][u];
	long str, stl;
	str = stl = u;

	for(int j = u; j < n; j++) // находим максимальное
		for(int a = 0; a < n; a++)
			if(fabs(A[j][a]) > fabs(max))
			{
				max = A[j][a];
				str = j;
				stl = a;
			}

	if(str != u)
	{
		double* N = A[u];  // меняем строки
		A[u] = A[str];
		A[str] = N;
	}

	if(stl < n-1 && stl != u)
	{
		for(long g = 0; g < n; g++)
		{
			double z = A[g][u];		// меняем столбцы
			A[g][u] = A[g][stl];
			A[g][stl] = z;
		}
	}

	if(A[u][u] != 0)
		return true;

	return false;
}

bool Gaus(double** Ar, int n, double *X)
{
	double **A = copyA(Ar, n);

	for(int u = 0; u < n; u++)
	{
		if(!fixz(A,u,n))
			return false;

		double d = A[u][u];
		for(long q = u; q < n+1; q++)
			A[u][q]/=d;

		for(int o = u+1; o < n; o++)
		{
			double del = A[o][u];
			for(long t = u; t < n+1; t++)
				A[o][t]-=A[u][t]*del;
		}
	}

	if(X == NULL)
		X = new double[n];

	X[n-1] = A[n-1][n];

	for(int v = n-2; v > -1; v--)
	{
		double s = 0;
		X[v] = A[v][n];
		for(int y = n-1; y > v; y--)
			s+=A[v][y]*X[y];
		X[v]-=s;
	}

    freeA(A,n);
	return true;
}

double** mallocA(int n)
{
	double** r = new double*[n];
	for(int i = 0; i < n; i++)
		r[i]=new double[n+1];
	return r;
}
void freeA(double** A, int n)
{
    for(int i = 0; i < n; i++)
		delete[] A[i];
	delete[] A;
}
double** copyA(double** A, int n)
{
	double** r = mallocA(n);
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n+1; j++)
			r[i][j] = A[i][j];

	return r;
}

void print(double **T, int n)
{
	cout << endl << "Matrix: \n";
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n+1; j++)
			cout << T[i][j] << "  ";
		cout << endl;
	}
}
