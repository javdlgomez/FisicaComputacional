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

    //===================================
    //
    // Metodo de Euler Mejorado y Modificado
    // compilacion: g++ -o programa ejercicio5-2.cpp
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







Solución:


Nos piden obtener una gráfica que compare los datos obtenidos de los distintos métodos numéricos contra la solución analítica. Para producir los resultados vamos a dirigir el output de consola a un archivo, mediante .\programa > datos"metodocorrespondiente".txt escritos en UTF-8 y graficar los puntos encontrados en gnuplot. 
enor: Euler normal
emod: Euler modificado
em: Euler Mejorado


![Solución analítica vs Numéricas](https://github.com/javdlgomez/FisicaComputacional/blob/main/Tarea%231/graficas/522.jpg)

Notemos que los nuevos métodos de Euler tienen mayor presición para este problema en específico, esto era lo que se esperaba ya que son refinamientos del original.




