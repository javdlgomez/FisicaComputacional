# Tarea#1

Solución de ejercicios utilizando distintas iteraciones del método de Euler para resolver EDOs en c++, solamente se utilizó la librería de stdio.

## Ejercicio 5.1: 

Nos piden resolver por medio del método de euler la siguiente ecuación diferencial:

$$ y'(x) = x^2+1 $$

Con Condiciones iniciales y de dominio
$$y(0) = 0, 0x<1 $$

Con distintos valores de h= 0.05, 0.10, 0.15 y 0.20 comparados con la solución analítica.

Código: 

    //===================================
    //
    // Metodo de Euler
    // compilacion: g++ -o programa ejercicio5-1.cpp
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

Solución:


Nos piden obtener distintas gráficas con los distintos pasos y graficarlas comparadas con la solución analítica.


![Alt text](https://github.com/javdlgomez/FisicaComputacional/blob/main/Tarea%231/graficas/d1.emf?raw=true "Title")
