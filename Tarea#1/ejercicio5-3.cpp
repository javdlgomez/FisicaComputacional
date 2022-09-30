//===================================
//
// Euler modificado y conservacion de la energia
// compilacion: g++ -o programa ejercicio5-3.cpp
//===================================

//este ejercicio es muy parecido al anterior
//solo debemos agregar el paso del método modificado
//y definir los parametros que deseamos 
#include <iostream>
#include <cmath>


using namespace std;

//Hago 2 metodos de derivadas distintos y cada metodo de euler llama 
//su respectivo metodo de derivada

double euler_modificado1( double y, double t, double h ,double v);
double euler_modificado2( double y, double t, double h, double v );
double derivada1( double y, double t );
double derivada2( double y, double t );
//Ademas queremos medir la energia
double energia(double y, double t);


int main()
{
  // Datos iniciales
  const double y0 = 0.0;
  const double v0 = 1;
  const double t0 = 0.0;
  const double h = 0.10;
  const double E0 = 0.5;
  const int N = 100; // numero de iteraciones
  
  double y = y0;
  double t = t0;
  double v = v0;
  double E = E0;
  double v_nueva = 0.0;
  double y_nueva = 0.0;
  double E_nueva = 0.0;


  cout << t << "\t" << y << "\t" << v << "\t" << E << endl;
  
  // ciclo de iteraciones
  for( int i=1; i<=N; i++ ){

    v_nueva = euler_modificado2(y,t,h,v);
    y_nueva = euler_modificado1(y, t, h, v);


    y = y_nueva;
    v = v_nueva;

    E_nueva = energia(y,v);
    E = E_nueva;

    t = t + h;

      cout << t << "\t" << y << "\t" << v << "\t" << E << endl;
  }
  
  return 0;
}

//Aquí sería mucho mejor retornar una tupla
//en vez de hacer 2 veces los mismos calculos
//para nada mas retornar una variable
//esta implementacion no es buena pero funciona


double euler_modificado1( double y, double t, double h, double v )
{
  //definimos los parámetros a utilizar

  double t_medio = t + h/2;
  double v_medio = v +h/2*derivada2(y,t);
  // V = V0 -hY0
  double y_medio = y +h/2*derivada1(v,t);
  //Y = Y0+hV0

  //esta se la parte del método modificado

  double y_imas1 = y +h*v_medio;
  // Y = V0 -hV(t+h/2)

  return y_imas1;
}



double euler_modificado2( double y, double t, double h , double v)
{
  //definimos los parámetros a utilizar

  double t_medio = t + h/2;
  double v_medio = v +h/2*derivada2(y,t);
  // V = V0 -hY0
  double y_medio = y +h/2*derivada1(v,t);
  //Y = Y0+hV0

  //esta se la parte del método modificado

  double v_imas1 = v -h*y_medio; 
  // V = V0 -hY(t+h/2)
  
  return v_imas1;
}


double derivada1( double v, double t )
{
  return v;
}


double derivada2( double y, double t )
{
  return -y;
}

double energia(double y, double v)
{
  double energia = y*y/2+v*v/2;
  return energia;
}