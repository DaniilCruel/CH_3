#include "head.h"




double func(double* u, int n, double tk)
{
	if (tk == 0)
	{
		tk = 0.000000001;
	}
	const int a = 1;
	switch (n)
	{
		case 0:
		return -u[0]*u[1]+sin(tk)/tk;
		break;
		case 1:
			return -u[1] * u[1] + a * tk / (1 + tk * tk);
		break;
	default:
		return 0;
		break;
	}
		
}

void EilerExplicit(double* u, int n)
{
	double Tau = 0;
	double Eps = 0.001, T = 1, TauMax = 0.5;
	double tk = 0;
	double* yk = new double[n];

	for(int i = 0; i < n; i++)
		yk[i] = u[i];

	ofstream out;
	out.open("Eiler_Explicit.txt", ios::trunc);
	
	out<< "t"  << setw(16);
	out<< "u1" << setw(16);
	out<< "u2" << endl << "=================================================\n";



	int Shag = 0;
	do{
		tk+=Tau;
		for(int i = 0; i < n; i++) { 
			if(i == 0)
				out<<setiosflags(ios::fixed) << tk << setw(15);
			out<< setiosflags(ios::fixed) <<setw(15)<<yk[i];
			if(i == n-1)
				out<< endl;
		}

		for (int t = 0; t < n - 1; t++)
		{
			if (Eps / (fabs(func(yk, t,tk)) + Eps / TauMax) > Eps / (fabs(func(yk, t + 1,tk)) + Eps / TauMax))
				Tau = Eps / (fabs(func(yk, t + 1,tk)) + Eps / TauMax);
			else
				Tau = Eps / (fabs(func(yk, t,tk)) + Eps / TauMax);
		}
		for (int i = 0; i < n; i++) 
		{ 
			yk[i] = yk[i] + Tau * func(yk, i,tk);
		
		}
		if (tk+Tau > T && tk < T)
			tk = T-Tau;

		Shag++;
	}while(tk < T);

	cout<<"Iterations quantity is "<< Shag <<endl;
	out<<endl<<"Iterations quantity is "<< Shag;
	out.close();

	delete [] yk;
}
