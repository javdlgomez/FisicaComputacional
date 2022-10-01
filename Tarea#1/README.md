# Tarea#1

Solución de ejercicios utilizando distintas iteraciones del método de Euler para resolver EDOs en c++, solamente se utilizó la librería de stdio.

## Ejercicio 5.1: 

Nos piden resolver por medio del método de euler la siguiente ecuación diferencial:

$$ y'(x) = x^2+1 $$

Con Condiciones iniciales y de dominio
$$y(0) = 0, 0≤x<1 $$

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

![image](https://user-images.githubusercontent.com/100542213/193376334-5ceb9547-282d-43fb-9ffc-b0f7694ec171.png)



Notemos que los valores más pequeños para iterar nuestro método numérico se asemejan más a la solución analítica.


## Ejercicio 5.2: 

Nos piden resolver el mismo problema del ejercicio 5.1 pero utilizando el método de Euler modificado, mejorado y comparar los resultados de los 3 métodos numéricos contra la solución analítica con h = 0.10.

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

Los nombres de los archivos tienen los siguientes significados

enor: Euler normal

emod: Euler modificado

em: Euler Mejorado


![image](https://user-images.githubusercontent.com/100542213/193376343-c34d52b2-8789-43af-ac7c-a71db6313300.png)

Notemos que los nuevos métodos de Euler tienen mayor presición para este problema en específico, esto era lo que se esperaba ya que son refinamientos del original.


## Ejercicio 5.3: 

Nos piden resolver un sistema de masa resorte por medio del Método de Euler Modificado con h = 0.10 y encontrar si se conserva la energía empleando este método numérico. Para ello primero debemos resolver el siguiente sistema de EDOs.

$$ y'(t) = v $$

$$ v'(t) = -x $$

Con Condiciones iniciales y de dominio
$$y(0) = 0, v(0) = 1,  t≥0 $$



Código: 

    //===================================
    //
    // Euler modificado y conservacion de la energia
    // compilacion: g++ -o programa ejercicio5-3.cpp
    //===================================


    #include <iostream>
    #include <cmath>


    using namespace std;

    //Hago 2 metodos de derivadas distintos y cada metodo de euler llama 
    //su respectivo metodo de derivada, esto no es una buena practica de programacion
    //pero lo hice asi para seguir con el mismo esquema del ejemplo pero
    //no hice bien la adaptacion

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
    
    
Solución:

Nos piden graficar las funciones de posición y velocidad contra el tiempo para analizar el movimiento del sistema y la energía contra el tiempo para encontrar si nuestro método numérico conserva el valor de esta. Para ello escribimos el output de consola en un archivo llamado datosej3 y estos los graficamos en gnuplot.

El orden de las columnas de nuestra base de datos es; t,y,v,E. Esto quiere decir que las gráficas son las siguientes:

datosej3 using 1:2; y vs t

datosej3 using 1:3; v vs t

datosej3 using 1:4; E vs t


![image](https://user-images.githubusercontent.com/100542213/193376353-6b4f14a9-b9a1-4362-8309-a2278337aa3c.png)

El valor de la energía parece conservarse, especialmente si lo comparamos con el método de Euler original. Además la solución del sistema se comporta como se espera que sea la solución analítica con respecto a la posición y velocidad.
    
    
    
