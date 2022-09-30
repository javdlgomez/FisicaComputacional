//===================================
//
// Metodo de Euler Mejorado y Modificado
//
//===================================
//este ejercicio es muy parecido al anterior
//solo debemos agregar el paso del método modificado
//y definir los parametros que deseamos 
#include <iostream>
#include <cmath>


using namespace std;


double euler_mejorado( double y, double t, double h );
double euler_modificado( double y, double t, double h );
double derivada( double y, double t );

int main()
{
  // Datos iniciales
  const double y0 = 0.0;
  const double t0 = 0.0;
  const double h = 0.10;
  // 1 = N*h
  const int N = 10; // numero de iteraciones
  
  double y = y0;
  double t = t0;
  double y_nueva = 0.0;

  cout << t << "\t" << y << endl;
  
  // ciclo de iteraciones
  for( int i=1; i<=N; i++ ){

    y_nueva = euler_mejorado( y, t, h );

    y = y_nueva;
    t = t + h;

      cout << t << "\t" << y << endl;
  }
  
  return 0;
}

double euler_modificado( double y, double t, double h )
{
  //definimos los parámetros a utilizar
  double t_medio = t + h/2;
  double y_medio = y +h/2*derivada(y,t);
  //esta se la parte del método modificado

  double y_imas1 = y +h *derivada(y_medio,t_medio); 
  
  return y_imas1;
}

double euler_mejorado( double y, double t, double h )
{
  double y_tilde = y + h*derivada( y, t );
  
  //esta se la parte del método mejorado

  double y_imas1 = y + 0.5*h * ( derivada( y, t ) + derivada( y_tilde, t+h ) ); 
  
  return y_imas1;
}


double derivada( double y, double t )
{
  return y*y+1;
}
