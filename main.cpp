/**********************************************************************/
/** Programa que implementa un algoritmo de resoluci�n de ecuaciones **/
/** diferenciales en derivadas parciales para resolver la ecuaci�n   **/
/** de Schr�dinger.                                                  **/
/** El programa pide el valor de N que ser�n las posiciones en las   **/
/** que se mida la funci�n de onda, el n�mero de ciclos n y la       **/
/** longitud de onda lambda.                                         **/
/** El programa escribe en un fichero los valores de la funci�n de   **/
/** onda en los instantes indicados y escribe en otro fichero la     **/
/** evoluci�n de la norma con el n�mero de ciclos. La cual debe      **/
/** mantenerse constante.                                            **/
/**********************************************************************/

#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include "Schrodinger.h"


using namespace std;

int main()
{
    double lambda,s,k0,norma;
    int ciclos,N;
    double *V,*modulo; //Potencial
    int j,t=0;//Indices el indice i queda reservado ya que nombro asi a la unidad imaginaria
    ofstream fich1,fich2,fich3;//Ficheros

    //Introducci�n de N
    cout<<"Introduzca N"<<endl;
    cin>>N;

    //Vectores necesarios para el c�lculo
    complex<double> phi[N];
    complex<double> a[N-1];
    complex<double> b[N-1];
    complex<double> chi[N];

    /** Recogida de parametros, inicializaci�n de memoria d�namica y de algunas variables iniciales**/

    //Inicializaci�n de la memoria d�namica

    V=new double [N];
    modulo=new double [N];

    //Introducci�n de la longitud de onda
    cout<<"Introduzca lambda"<<endl;
    cin>>lambda;

    fich1.open("F�ncion de onda.txt");
    fich2.open("Norma.txt");
    fich3.open("Probabilidad.txt");

    //Introducci�n del n�mero de ciclos
    do
    {
        cout<<"Introduzca el numero de ciclos a calcular entre 1 y "<<N/4 <<endl;
        cin>>ciclos;
    }
    while ((ciclos<1) || (ciclos >N/4));



    //Calculos necesarios para generar variables y condiciones
    k0=2.0*PI*ciclos/N;
    s=1.0/(4.0*k0*k0);

    //Generaci�n de las funciones de onda iniciales y las alfas
    generador(s,V,k0,N,lambda,phi,a);

    //Calculo la norma inicial
    norma=0.0;
    for(j=0;j<N;j++)
        modulo[j]=norm(phi[j]);
    for(j=0;j<N;j++)
        norma+=modulo[j];

    //Normalizo la funci�n de onda
    for(j=0;j<N;j++)
        phi[j]=phi[j]/sqrt(norma);

    //Esto es puramente estetico y es para tener el fichero norma con 1 en todos los pasos
    for(j=0;j<N;j++)
    {
        modulo[j]=norm(phi[j]);
        fich3<<j<<" "<<modulo[j]<<"\t";
    }
    fich3<<endl;
    norma=0.0;
    for(j=0;j<N;j++)
        norma+=modulo[j];
    fich2<<norma<<endl;

    //Escribo la funcion de onda inicial
    for(j=0;j<N;j++)
            fich1<<phi[j]<<"\t";
        fich1<<endl;



    /** Implementaci�n del algoritmo **/
    while (t<ciclos)
    {

        //Paso 2: Calculo de beta
        Beta(s,V,N,b,phi,a);

        //Paso 3: Calculo de chi
        Chi(N,a,b,chi);

        //Paso 4: Actualizaci�n de las funciones de onda
        Phi(N,chi,phi,modulo);

        //Guardo los modulos de las funciones de onda al cuadrado para su futura representaci�n
        for(j=0;j<N;j++)
            fich3<<j<<" "<<modulo[j]<<"\t";
        fich3<<endl;




        //A modo de comprobaci�n miro a ver si la norma se mantiene constante
        norma=0.0;

        for(j=0;j<N;j++)
            norma+=modulo[j];
        fich2<<norma<<endl;

        //Para ver la evoluci�n temporal de la funci�n de onda voy a escribirla en un fichero
        for(j=0;j<N;j++)
            fich1<<phi[j]<<"\t";
        fich1<<endl;

        t++;
    }

    /** Borrado de la memoria d�namica y fin del programa **/

    //Borrado de memoria d�namica
    delete[] V;
    delete[] modulo;

    //Cierre de ficheros
    fich1.close();
    fich2.close();
    fich3.close();


    return 0;
}
