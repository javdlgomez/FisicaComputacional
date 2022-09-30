# Tarea#1

Solución de ejercicios utilizando distintas iteraciones del método de Euler para resolver EDOs en c++, solamente se utilizó la librería de stdio.

## Ejercicio 5.1: 

Nos piden resolver por medio del método de euler la siguiente ecuación diferencial:

$$ y'(x) = x^2+1 $$

Con Condiciones iniciales y de dominio
$$y(0) = 0, 0x<1 $$

Con distintos valores de h1= 0.05, h2= 0.10,h3= 0.15 y h4= 0.20 comparados con la solución analítica.

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


Nos piden obtener distintas gráficas con los distintos pasos y comparadas con la solución analítica, en el repositorio se pueden visualizar el resto de las gráficas. Para producir los resultados vamos a dirigir el output de consola a un archivo, mediante .\programa > datoshcorrespondiente.txt escritos en UTF-8 y graficar los puntos encontrados en gnuplot. 

![Solución analítica vs Numéricas](https://github.com/javdlgomez/FisicaComputacional/blob/main/Tarea%231/graficas/djun.jpg)

Notemos que los valores más pequeños para iterar nuestro método numérico se asemejan más a la solución analítica.


## Ejercicio 5.2: 

Nos piden resolver el mismo problema del ejercicio 5.1 pero utilizando el método de Euler modificado, mejorado y comparar los resultados de los 3 métodos numéricos contra la solución analítica.

Código:

