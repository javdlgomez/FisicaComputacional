//============================================
//
// Metodo de RK4 caida con friccion
//
//============================================

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>


using namespace std;





void RK4( const double *y,
		     const int n_ec,
		     const double t,
		     const double h,
		     double *y_imas1,
	             void (*derivada)( const double *, const double, double * ) );
void salidaSolucion( const double t, const double *y, const int N );
void caidaLibreFriccion( const double *y, const double t, double *dydt );





int main()
{
  // Datos iniciales
  const double t0 = 0.0;
  const double h = .05;
  const int N = 10000; // numero de iteraciones
  const int out_cada = 1; // output cada out_cada iteraciones
  const int n_ec = 2; // numero de ecuaciones


  // reservar espacio para y
  double *y       = new double[ n_ec ];
  double *y_nueva = new double[ n_ec ];

  // inicializar cada variable segun las condiciones iniciales
  y[0]  = 1.0;
  y[1] = 0.0;
    

  // puntero a la funcion "derivada"
  void (*derivada)( const double *, const double, double * );
  derivada = caidaLibreFriccion;
  

  // inicializar y_nueva
  for( int i=0; i<n_ec; i++ ) y_nueva[i] = 0.0;
  
  double t = t0;
  

  salidaSolucion( t, y, n_ec );
  
  // ciclo de iteraciones
  for( int i=1; i<=N; i++ ){
    RK4( y, n_ec, t, h, y_nueva, derivada );

    y = y_nueva;
    t = t + h;

    if ( i%out_cada == 0){
      salidaSolucion( t, y, n_ec );
    }
    
  }
  
  return 0;
}







void salidaSolucion( const double t, const double *y, const int N )
{
  cout << fixed << setprecision(4) << t;

  for( int i=0; i<N; i++ )
    cout << scientific << setprecision(4) << "\t" << y[i];

  cout << endl;  
}








void RK4( const double *y,
		     const int n_ec,
		     const double t,
		     const double h,
		     double *y_imas1,
		     void (*derivada)( const double *, const double, double * ) )
{
  double *k0 = new double[ n_ec ];
  double *k1 = new double[ n_ec ];
  double *k2 = new double[ n_ec ];
  double *k3 = new double[ n_ec ];
  double *z  = new double[ n_ec ];

  (*derivada)( y, t, k0 );

  for( int i=0; i<n_ec; i++ )
    z[i] = y[i] + 0.5*k0[i]*h;

  (*derivada)( z, t+0.5*h, k1 );

  for( int i=0; i<n_ec; i++ )
    z[i] = y[i] + 0.5*k1[i]*h;

  (*derivada)( z, t+0.5*h, k2 );

  for( int i=0; i<n_ec; i++ )
    z[i] = y[i] + k2[i]*h;

  (*derivada)( z, t+h, k3 );

  for( int i=0; i<n_ec; i++ )
   y_imas1[i] = y[i] + h/6.0 * ( k0[i] + 2*k1[i] + 2*k2[i] + k3[i] );

  delete[] k0;
  delete[] k1;
  delete[] k2;
  delete[] k3;
  delete[] z;
}


void caidaLibreFriccion( const double *y, const double t, double *dydt )
{
  const double g = 9.8; // aceleracion gravedad
  const double m = 0.01; // masa
  const double k = .0001; // constante de friccion
  
  dydt[0] = y[1];
  dydt[1] = m*g - k*y[1]*y[1];
}






void derivada( const double y, const double t, double &dydt )
{
  const double g = 9.8; // aceleracion gravedad
  const double m = .01; // masa
  const double k = .0001; // constante de friccion
  
  dydt = m*g - k*y;
}
