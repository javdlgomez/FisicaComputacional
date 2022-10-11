# Tarea#2

Solución de ejercicios utilizando distintas iteraciones del método de Runge-Kutta a orden 4 en c++.

## Ejercicio 5.5: 

Nos piden resolver por medio del método de RK4 la siguiente ecuación diferencial:

$$ p'(x) = mg-kv^2$

Con un tamaño de paso para obtener 4 cifras significativas en un intervalo de 10 segundos y comparar con la solución analítica a O(0).

Código: 

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

  
Solución:


Nos piden obtener distintas gráficas con los distintos pasos y comparadas con la solución analítica, en el repositorio se pueden visualizar el resto de las gráficas. Para producir los resultados vamos a dirigir el output de consola a un archivo, mediante .\programa > datoshcorrespondiente.txt escritos en UTF-8 y graficar los puntos encontrados en gnuplot. 

Nos piden realizar una comparación entre la velocidad de la solución analítica de un objeto en caída libre contra la solución numérica por medio de RK4 de un objeto cayendo tomando en cuenta la resistencia del aire a O(2). Para producir los resultados vamos a dirigir el output de consola a un archivo, mediante .\programa > datos#ejercicio.txt escritos en UTF-8 y graficar los puntos encontrados en gnuplot. 

![image](https://github.com/javdlgomez/FisicaComputacional/blob/main/Tarea%232/imagenes/ej5-5corto.emf)


El sistema con resistencia del aire llega a una velocidad terminal y la solución de caída libre no estpa acotado.


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
    
