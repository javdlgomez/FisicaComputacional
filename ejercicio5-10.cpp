///============================================
//
// Oscilacion de Van der Pol con Runge-Kutta
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


void VanderPol( double t, double *y, double *dydt );




int main()
{
  int nvar, nok, nbad;
  double t1, t2, eps, h, hmin, *ystart;


  /* memory space for variables */
  // numero de ecuaciones
  nvar = 3;

  // valor inicial de cada variable
  ystart = new double[ nvar ];

  /* other variables initialization */
  // tolerancia (error)
  eps  = .00001;
  h    = .001;
  hmin = 1e-5;
  nok  = 0;
  nbad = 0;

  

  /* initial condition */
  ystart[0]  = 0.5;
  ystart[1]  = 0.0;
  ystart[2] = -0.5;

  // tiempo inicial
  t1 =  0.0;

  // tiempo final
  t2 =  20.0;
  

  odeint( ystart, nvar, t1, t2, eps, h, hmin, &nok, &nbad, &VanderPol, &rkqs );

  cout << "nok = " << nok <<"\t nbad = " << nbad << endl;
  
  return 0;
}







void VanderPol( double t, double *y, double *dydt )

    
{
  dydt[0] = y[1];
  dydt[1] = y[2];
  dydt[2] = -y[0] - y[1]*(y[0]*y[0]-1);

}
