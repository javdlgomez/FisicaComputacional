//===================================
//
// Metodo de Euler
//
//===================================

#include <iostream>

using namespace std;

//utilizamos la implementación más sencilla de este método vista en clase
double euler( double y, double t, double h );
double derivada( double y, double t );

int main()
{
  // Datos iniciales
  const double y0 = 0;
  const double t0 = 0;
  const double h = 0.1;

  //Para calcular N despejamos la ecuación
  // 1 = (N-2)h
  //el -2 aparece porque empezamos a contar en 2 y nos dicen que el dominio es menor a 1
  const int N = 9; // 
  
  double y = y0;
  double t = t0;
  double y_nueva = 0.0;

  cout << t << "\t" << y << endl;
  
  // ciclo de iteraciones
  for( int i=0; i<=N; i++ ){

    y_nueva = euler( y, t, h );

    y = y_nueva;
    t = t + h;

    cout << t << "\t" << y << endl;
  }
  
  return 0;
}


double euler( double y, double t, double h )
{
  return y + h*derivada( y, t );
}


double derivada( double y, double t )
{
  //aquí debemos tener cuidado para cambiar al problema
  //específico que deseamos solucionar
  return y*y+1;
}
