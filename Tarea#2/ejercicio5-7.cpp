//============================================
//
// Pendulo doble compuesto con Runge-Kutta
// usando tamaño de paso adaptativo
//
//============================================

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "rkqs.hpp"
#include "odeint.hpp"

using namespace std;


void caidaLibreFriccion( double t, double *y, double *dydt );




int main()
{
  int nvar, nok, nbad;
  double t1, t2, eps, h, hmin, *ystart;


  /* memory space for variables */
  // numero de ecuaciones
  nvar = 2;

  // valor inicial de cada variable
  ystart = new double[ nvar ];

  /* other variables initialization */
  // tolerancia (error)
  eps  = 1e-10;
  h    = 1;
  hmin = 1e-5;
  nok  = 0;
  nbad = 0;

  

  /* initial condition */
  ystart[0]  = 1.0;
  ystart[1]  = 0.0;

  // tiempo inicial
  t1 =  0.0;

  // tiempo final
  t2 =  2000.0;
  

  odeint( ystart, nvar, t1, t2, eps, h, hmin, &nok, &nbad, &caidaLibreFriccion, &rkqs );

  cout << "nok = " << nok <<"\t nbad = " << nbad << endl;
  
  return 0;
}







void caidaLibreFriccion( double t, double *y, double *dydt )
{
  const double g = 9.8; // aceleracion gravedad
  const double m = 0.01; // masa
  const double k = .0001; // constante de friccion
    

  dydt[0] = y[1];
  dydt[1] = m*g - k*y[1]*y[1];

}
