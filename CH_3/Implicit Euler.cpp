#include "head.h"


double f1(double* uk1, double* uk, double Tau, double l1, double l2, double l3)
{
return uk1[0]-uk[0]-Tau*(  uk[0] * (2 * l1 + 4 * l2) + uk[1] * 2 * (l1 - l2) + uk[2] * 2 * (l1 - l2) + (4*l1+2*l2)    )/6;
}

double f2(double *uk1, double *uk, double Tau, double l1, double l2, double l3)
{
return uk1[1]-uk[1]-Tau*(  uk[0] * (2 * l1 - 2 * l2) + uk[1] * (2*l1 + l2 + 3*l3) + uk[2] * (2 * l1 + l2 - 3*l3) + (4*l1 - l2 - 9*l3)  )/6;
}

double f3(double *uk1, double *uk, double Tau, double l1, double l2, double l3)
{
return uk1[2]-uk[2]-Tau*(  uk[0] * (2 * l1 -2 * l2) + uk[1] * (2*l1 + l2 - 3*l3) + uk[2] * (2*l1 + l2 + 3*l3) + (4*l1 - l2 + 9*l3)  )/6;
}

typedef double(*pf)(double*, double*, double, double, double, double);

double Differential(pf f, double* uk1, double* uk, double Tau, int n, double l1, double l2, double l3)
{
	double dx = 1e-9;
	double* D = new double[3];
	for (int i = 0; i < 3; i++)
	{
		
		if (i == n - 1)
		{

			D[i] = uk1[i] + dx;
		}
		else
		D[i] = uk1[i];
	}

	double F = f(uk1, uk, Tau, l1, l2, l3);
	double dF = f(D, uk, Tau, l1, l2, l3);

	delete[] D;
	return (dF - F) / dx;
}

double Nuton(double *uk1, double *uk, double Tau, int n, double l1, double l2, double l3)
{
	double e = 1e-9;
	double **A = mallocA(n);
	double* R = new double[n];
	double b1, b2 = 0;

	int i=0;
	do{
		i++;

		A[0][0]= Differential(f1, uk1, uk, Tau, 1, l1, l2, l3); // производна€ по u1
		A[0][1]= Differential(f1, uk1, uk, Tau, 2, l1, l2, l3); // производна€ по u2
		A[0][2]= Differential(f1, uk1, uk, Tau, 3, l1, l2, l3); // производна€ по u3
		A[0][3]=-f1(uk1,uk,Tau, l1, l2, l3);
		A[1][0]= Differential(f2, uk1, uk, Tau, 1, l1, l2, l3); // производна€ по u1
		A[1][1]= Differential(f2, uk1, uk, Tau, 2, l1, l2, l3);
		A[1][2]= Differential(f2, uk1, uk, Tau, 3, l1, l2, l3); // производна€ по u3
		A[1][3]=-f2(uk1,uk,Tau, l1, l2, l3);
		A[2][0]= Differential(f3, uk1, uk, Tau, 1, l1, l2, l3); // производна€ по u1
		A[2][1]= Differential(f3, uk1, uk, Tau, 2, l1, l2, l3); // производна€ по u2
		A[2][2]= Differential(f3, uk1, uk, Tau, 3, l1, l2, l3);
		A[2][3]=-f3(uk1,uk,Tau, l1, l2, l3);
		
		if(!Gaus(A, n, R))
		{
			cout<<"Error! Metod Gaussa\n";
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					cout << A[i][j] << setw(12);
				}
				cout << R[i] << endl;
			}
			return -1;
			break;
		}
		
		uk1[0]+=R[0];	uk1[1]+=R[1]; uk1[2]+=R[2];

		if(fabs(f1(uk1,uk,Tau,l1, l2, l3)) > fabs(f2(uk1,uk,Tau, l1, l2, l3)) && fabs(f1(uk1,uk,Tau, l1, l2, l3)) > fabs(f3(uk1,uk,Tau, l1, l2, l3))) //  ритерием завершени€ итерационного процесса из методички) eps1
			b1 = fabs(f1(uk1,uk,Tau, l1, l2, l3));
		else if(fabs(f2(uk1,uk,Tau, l1, l2, l3)) > fabs(f3(uk1,uk,Tau, l1, l2, l3)))
			b1 = fabs(f2(uk1,uk,Tau, l1, l2, l3));
		else
			b1 = fabs(f3(uk1,uk,Tau, l1, l2, l3));

		for(int k = 0; k < n; k++) //  ритерием завершени€ итерационного процесса из методички eps2
		{
			if(fabs(uk1[k]) < 1)
				b2 = fabs(R[k]);
			else
				if(fabs(uk1[k]) >= 1)
					b2 = fabs(R[k]/uk1[k]);
		}

	}while(b1>e && b2>e);

	delete[] R;
	freeA(A,n);
	return 0;
}


void EilerImplicit_(double* u, int n)
{
    double Tau, Tau_1, Tau1;
	double Eps = 1e-5, T = 1, TauMax = 0.1, TauMin = 0.01;
	double tk = 0, tk1 = 0;
	double l1, l2, l3;
	double* ek = new double[n];
	double* uk = new double[n];                //y k
	double* uk_1 = new double[n];              //y k-1
	double* uk1 = new double[n];               //y k+1
	Tau = Tau_1 = TauMin;

	cout << "¬ведите значени€ L1 L2 L3:\n";
	cin >> l1 >> l2 >> l3;
	for(int i = 0; i < n; i++) /// присваиваем наше u
		uk[i] = uk_1[i] = uk1[i] = u[i];


		///
	ofstream out;
	out.open("Eiler_Implicit.txt", ios::trunc);


	int shag = 0;
	do{

		Nuton(uk1, uk, Tau, n, l1, l2, l3); 
		
		if (Nuton(uk1, uk, Tau, n, l1, l2, l3) == -1)
		{
			break;
		}
		for(int k = 0; k < n; k++)  // ¬ычисл€ем  ek по формуле   3.16
			ek[k]=(-Tau/(Tau+Tau_1))*(uk1[k]-uk[k]-(Tau/Tau_1)*(uk[k]-uk_1[k]));

		bool F=false;
		for(int k = 0; k<n; k++)  //≈сли  ek > Eps, то выполн€ем и переходим в начало.
			if(fabs(ek[k]) > Eps && !F)
			{
				Tau/=2;
				tk1=tk;
				for(int j = 0; j < n; j++)
					uk1[j] = uk[j];
				F=true;
			}
		if(F)
			continue;


		double *tmp = new double[n];
		for(int f = 0; f < n; f++)  // ќпределить шаг   Tau1 по формулам (3.18)  » выбираетс€ минимальным среди всех значений  Tau1
		{
			if(fabs(ek[f]) > Eps)
				tmp[f] = Tau/2;
			if(Eps/4 < fabs(ek[f]) && fabs(ek[f]) <= Eps)
				tmp[f]=Tau;
			if(fabs(ek[f]) <= Eps/4)
				tmp[f]=2*Tau;
		}

		for(int t = 0; t < n-1; t++) 
		{
			if(tmp[t] > tmp[t+1])
				Tau1 = tmp[t+1];
			else
				Tau1 = tmp[t];
		}
		delete[] tmp;

		if(Tau1 > TauMax) // ≈сли Tau1 > TauMax, то присваиваем Tau1=TauMax
			Tau1 = TauMax;

		if((tk+Tau > T) && tk < T) // если tk перепрыгнула интервал T, то смещаем в “
		{
			tk = T;
			Tau = T-tk1;
		}

		tk1 = tk+Tau;

		for(int i = 0; i < n; i++) // записываем в файл, в консоли весь вывод не влезет
		{
			if(i == 0)
				out<<setiosflags(ios::fixed) << tk<< setw(15);
			out<< setiosflags(ios::fixed) <<setw(15)<<uk[i];
			if(i == n-1)
				out<< endl;
		}
		//

		//¬ыполн€ем сдвиг переменных и шагов интегрировани€
		for(int s = 0; s < n; s++)
		{
			uk_1[s] = uk[s];
			uk[s] = uk1[s];
		}
		Tau_1 = Tau;
		Tau = Tau1;
		tk = tk1;

		shag++;
	}while(tk < T);

	cout<<endl<<"Iterations quantity is "<< shag <<endl;

	out<<endl<<"Iterations	quantity is "<< shag;
	out.close();

	delete[] ek;
	delete[] uk;
	delete[] uk_1;
	delete[] uk1;
}

