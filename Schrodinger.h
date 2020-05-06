#ifndef SCHRODINGER_H_INCLUDED
#define SCHRODINGER_H_INCLUDED

#include <complex>

using namespace std;

static const double PI= 3.14159265359; //Definicion de pi
static const double h=6.62607015e-34; //Constante de Planck
static const complex<double> i(0.0,1.0);//Definición de la unidad imaginaria


void generador (double s, double V[], double k0,int N,double lambda,complex<double> phi[],complex<double>a[]);
void Beta (double s, double V[],int N,complex <double> b[],complex<double>phi[],complex<double>a[]);
void Chi (int N,complex<double>a[],complex<double>b[],complex<double>chi[]);
void Phi (int N,complex<double>chi[],complex<double>phi[],double modulo[]);


#endif // SCHRODINGER_H_INCLUDED
