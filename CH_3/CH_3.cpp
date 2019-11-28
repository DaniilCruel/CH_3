#include "head.h"


int main()
{

	setlocale(0, "rus");
	int n = 2;
	double* u = new double[n];


	u[0] = 0;
	u[1] = -0.412;


	cout << "\nЯвный метод:\n";
	EilerExplicit(u, n);
	double* u1 = new double[n+1];
	u1[0] = 10;  u1[1] = 22;  u1[2] = 9;
	cout << "\nНеявный метод:\n";
	EilerImplicit_(u1, n+1);


	delete[]u;
	return 0;
}
