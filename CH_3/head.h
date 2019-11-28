#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;


void EilerExplicit(double* y, int n);

void EilerImplicit_(double* y, int n);


bool Gaus(double** Ar, int n, double *X);
double** mallocA(int n);
void freeA(double** A, int n);
double** copyA(double** A, int n);
void print(double **T, int n);


