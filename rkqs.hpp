#ifndef RKQS_H
#define RKQS_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

//#include <math.h>
//#include <stdio.h>
//#include <stdlib.h>

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

using namespace std;

void rkqs( double y[], double dydx[], int n, double *x, double htry,
	   double eps, double yscal[], double *hdid, double *hnext,
	   void (*derivs)(double, double[], double[]) )
{
  void rkck( double y[], double dydx[], int n, double x, double h,
	     double yout[], double yerr[],
	     void (*derivs)(double, double[], double[] ));
  
  int i;
  double errmax, h, htemp, xnew, *yerr, *ytemp;

  yerr = new double[n];
  ytemp = new double[n];
  h = htry;

  for(;;){
    rkck( y, dydx, n, *x, h, ytemp, yerr, derivs );
    errmax = 0.0;
    for(i=0;i<n;i++) errmax = fmax( errmax, fabs(yerr[i]/yscal[i]) );
    errmax /= eps;
    if (errmax <= 1.0) break;
    htemp = SAFETY*h*pow(errmax, PSHRNK);
    h = ( h >= 0.0 ? fmax( htemp, 0.1*h) : fmin( htemp, 0.1*h ));
    xnew = (*x)+h;
    if (xnew == *x ) printf( "stepsize undewflow in rkqs\n" );
  }

  if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax, PGROW );
  else *hnext = 5.0*h;
  *x += (*hdid=h);
  for(i=0;i<n;i++) y[i]=ytemp[i];
  delete[] ytemp;
  delete[] yerr;
}


void rkck( double y[], double dydx[], int n, double x, double h,
	   double yout[], double yerr[],
	   void (*derivs)(double, double[], double[] ))
{
  int i;
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  double *ak2, *ak3, *ak4, *ak5, *ak6, *ytemp;

  ak2 = new double[n];
  ak3 = new double[n];
  ak4 = new double[n];
  ak5 = new double[n];
  ak6 = new double[n];
  ytemp = new double[n];

  
  for(i=0;i<n;i++)
    ytemp[i] =y[i] + b21*h*dydx[i];

  (*derivs)( x+a2*h, ytemp, ak2 );
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);

  (*derivs)(x+a3*h,ytemp,ak3);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);

  (*derivs)(x+a4*h,ytemp,ak4);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);

  (*derivs)(x+a5*h,ytemp,ak5);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);

  (*derivs)(x+a6*h,ytemp,ak6);
  for (i=0;i<n;i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);

  for (i=0;i<n;i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);

  delete[] ak6;
  delete[] ak5;
  delete[] ak4;
  delete[] ak3;
  delete[] ak2;
}

#endif
