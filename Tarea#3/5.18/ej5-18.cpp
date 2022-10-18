//============================================
//
// Metodo de Euler Mejorado para movimiento
// gravitacional en 2 dimensiones
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
void movGrav2D_3cuerpos(  const double *y, const double t, double *dydt );




int main()
{
  // Datos iniciales
  const double t0 = 0.0;
  const double h = 1;
  const int N = 200000; // numero de iteraciones
  const int out_cada = 100; // output cada out_cada iteraciones
  const int n_ec = 12; // numero de ecuaciones
  //estos parametros nos permiten jugar con las condiciones iniciales de la posicion del cohete
  const int m = 0;
  const double theta_init = m*M_PI_4;
  const double v0 = 50000; // velocidad inicial

  // Archivo que guarda la energia total


  // reservar espacio para y
  double *y       = new double[ n_ec ];
  double *y_nueva = new double[ n_ec ];

  // inicializar cada variable segun las condiciones iniciales
  y[0]  = 0.0; //tierra x
  y[1]  = 0.0; // tierra y
  y[2]  = 3.84e8; //luna x
  y[3]  = 0.0; //luna y
  y[4]  = (6800e3)*(cos(theta_init)); //cohete x
  y[5]  = (6800e3)*(sin(theta_init)); //cohete y
  y[6]  = 0.0; // tierra vx
  y[7]  = 0.0; //tierra cy
  y[8]  = 0.0;//luna vx
  y[9]  = 1033;//luna vy (promedio de rango de velocidades)
  y[10] = v0*cos(theta_init+.04); //cohete vx
  y[11] = v0*sin(theta_init+.04); //cohete vy

  //de esta forma podemos con parametros theta y v0 cambiar las condiciones iniciales
    

  // puntero a la funcion "derivada"
  void (*derivada)( const double *, const double, double * );
  derivada = movGrav2D_3cuerpos;
  

  // inicializar y_nueva
  for( int i=0; i<n_ec; i++ ) y_nueva[i] = 0.0;
  
  double t = t0;
  

  salidaSolucion( t, y, n_ec );
  
  // ciclo de iteraciones
  for( int i=1; i<=N; i++ ){
    RK4( y, n_ec, t, h, y_nueva, derivada );
    //RK4( y, n_ec, t, h, y_nueva, derivada );

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
  cout << fixed << setprecision(3) << t;

  for( int i=0; i<N; i++ )
    cout << scientific << setprecision(9) << "\t" << y[i];

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
  double dx = (y[4]-y[6])*(y[4]-y[6]);
  double dy = (y[5]-y[7])*(y[5]-y[7]);
  double dis = sqrt(dx+dy);


  (*derivada)( y, t, k0 );

  for( int i=0; i<n_ec; i++ )
    z[i] = y[i] + 0.5*k0[i]*h;

  (*derivada)( z, t+0.5*h, k1 );

  for( int i=0; i<n_ec; i++ )
    z[i] = y[i] + 0.5*k1[i]*h;

  (*derivada)( z, t+0.5*h, k2 );

  for( int i=0; i<n_ec; i++ )
    z[i] = y[i] + k2[i]*h;


    
    
    ;

  (*derivada)( z, t+h, k3 );

  for( int i=0; i<n_ec; i++ )

   y_imas1[i] = y[i] + h/6.0 * ( k0[i] + 2*k1[i] + 2*k2[i] + k3[i] );
    //if ( pow(((y[4]-y[6]^2+(y[5]-y[7]^2)),.5))){


    

  delete[] k0;
  delete[] k1;
  delete[] k2;
  delete[] k3;
  delete[] z;
}



void movGrav2D_3cuerpos(  const double *y, const double t, double *dydt )
{
  const double m1 = 5.97e24; //masa tierra
  const double m2 = 7.34e22; //masa luna es 0 para este problema
  const double m3 = 1000; //masa nave
  const double G  = 6.66e-11; // Constante de gravitacion universal

  // distancias
  const double r21_3 = pow( pow(y[2]-y[0],2) + pow(y[3]-y[1],2), 1.5 );
  const double r31_3 = pow( pow(y[4]-y[0],2) + pow(y[5]-y[1],2), 1.5 );
  const double r32_3 = pow( pow(y[4]-y[2],2) + pow(y[5]-y[3],2), 1.5 );

  dydt[0]  = y[6];
  dydt[1]  = y[7];
  dydt[2]  = y[8];
  dydt[3]  = y[9];
  dydt[4]  = y[10];
  dydt[5]  = y[11];
  dydt[6]  = -G*m2*( y[0]-y[2] ) / r21_3 - G*m3*( y[0]-y[4] ) / r31_3;
  dydt[7]  = -G*m2*( y[1]-y[3] ) / r21_3 - G*m3*( y[1]-y[5] ) / r31_3;
  dydt[8]  = -G*m1*( y[2]-y[0] ) / r21_3 - G*m3*( y[2]-y[4] ) / r32_3;
  dydt[9]  = -G*m1*( y[3]-y[1] ) / r21_3 - G*m3*( y[3]-y[5] ) / r32_3;
  dydt[10] = -G*m1*( y[4]-y[0] ) / r31_3 - G*m2*( y[4]-y[2] ) / r32_3;
  dydt[11] = -G*m1*( y[5]-y[1] ) / r31_3 - G*m2*( y[5]-y[3] ) / r32_3;

}
