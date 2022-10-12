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


Nos piden realizar una comparación entre la velocidad de la solución analítica de un objeto en caída libre contra la solución numérica por medio de RK4 de un objeto cayendo tomando en cuenta la resistencia del aire a O(2). Para producir los resultados vamos a dirigir el output de consola a un archivo, mediante .\programa > datos#ejercicio.txt escritos en UTF-8 y graficar los puntos encontrados en gnuplot. 

![image](https://github.com/javdlgomez/FisicaComputacional/blob/main/Tarea%232/imagenes/ej5-5lite.png)


El sistema con resistencia del aire llega a una velocidad terminal a diferencia de la solución de caída libre que no está acotada, se generaron imágenes que muestran de mejor manera el comportamiento pero son muy pesadas para desplegarlas.


## Ejercicio 5.7: 

Nos piden resolver el mismo problema del ejercicio 5.5 pero utilizando el método adaptativo que vimos en clase con una tolerancia de error de 10-5.

Código:

Adjunto únicamente el archivo cpp ya que los headers son equivalentes a los vistos en clase


    //============================================
    //
    // Pendulo doble compuesto con Runge-Kutta
    // usando tamaño de paso adaptativo
    //
    //============================================

    #include <iostream>
    #include <cmath>
    #include <iomanip>
    #include <fstream>

    #include "rkqs.hpp"
    #include "odeint.hpp"

    using namespace std;


    void caidaLibreFriccion( double t, double *y, double *dydt );




    int main()
    {
      int nvar, nok, nbad;
      double t1, t2, eps, h, hmin, *ystart;


      /* memory space for variables */
      // numero de ecuaciones
      nvar = 2;

      // valor inicial de cada variable
      ystart = new double[ nvar ];

      /* other variables initialization */
      // tolerancia (error)
      eps  = 1e-10;
      h    = 1;
      hmin = 1e-5;
      nok  = 0;
      nbad = 0;



      /* initial condition */
      ystart[0]  = 1.0;
      ystart[1]  = 0.0;

      // tiempo inicial
      t1 =  0.0;

      // tiempo final
      t2 =  2000.0;


      odeint( ystart, nvar, t1, t2, eps, h, hmin, &nok, &nbad, &caidaLibreFriccion, &rkqs );

      cout << "nok = " << nok <<"\t nbad = " << nbad << endl;

      return 0;
    }







    void caidaLibreFriccion( double t, double *y, double *dydt )
    {
      const double g = 9.8; // aceleracion gravedad
      const double m = 0.01; // masa
      const double k = .0001; // constante de friccion


      dydt[0] = y[1];
      dydt[1] = m*g - k*y[1]*y[1];

    }





Solución:


Empleamos el algoritmo proporcionado por los headers vistos en clase y con ello adaptamos la solución para el problema de caída con resistencia cuadrática del aire.
Los resultados de las gráficas son generados en un archivo solucion.dat y las iteraciones fallidas o aceptadas son impresas en consola.


Obtuvimos los siguientes resultados para un h = 0.5:

nok = 117        nbad = 0


Ahora no necesitamos tanto espacio en memoria para generar una gráfica que muestre lo deseado, se aumentó el ancho de línea de la solución analítica para que se pudiera visualizar.


![image](https://github.com/javdlgomez/FisicaComputacional/blob/main/Tarea%232/imagenes/ej5-7.png)


## Ejercicio 5.10: 

Nos piden resolver el oscilador de Van der Pol por medio del método de RK4 adaptativo, este sistema está descrito por la siguiente ecuación diferencial:

$$ x''(t) = -x-\epsilon(x^2-1)x' $$


Con Condiciones iniciales y de dominio
$$y(0) = 0.5, x'(0) = 0.0,  0≥t≥8\pi $$



Código: 

    ///============================================
    //
    // Oscilacion de Van der Pol con Runge-Kutta
    // usando tamaño de paso adaptativo
    //
    //============================================

    #include <iostream>
    #include <cmath>
    #include <iomanip>
    #include <fstream>

    #include "rkqs.hpp"
    #include "odeint.hpp"

    using namespace std;


    void VanderPol( double t, double *y, double *dydt );




    int main()
    {
      int nvar, nok, nbad;
      double t1, t2, eps, h, hmin, *ystart;


      /* memory space for variables */
      // numero de ecuaciones
      nvar = 3;

      // valor inicial de cada variable
      ystart = new double[ nvar ];

      /* other variables initialization */
      // tolerancia (error)
      eps  = .00001;
      h    = .001;
      hmin = 1e-5;
      nok  = 0;
      nbad = 0;



      /* initial condition */
      ystart[0]  = 0.5;
      ystart[1]  = 0.0;
      ystart[2] = -0.5;

      // tiempo inicial
      t1 =  0.0;

      // tiempo final
      t2 =  20.0;


      odeint( ystart, nvar, t1, t2, eps, h, hmin, &nok, &nbad, &VanderPol, &rkqs );

      cout << "nok = " << nok <<"\t nbad = " << nbad << endl;

      return 0;
    }







    void VanderPol( double t, double *y, double *dydt )


    {
      dydt[0] = y[1];
      dydt[1] = y[2];
      dydt[2] = -y[0] - y[1]*(y[0]*y[0]-1);

    }

Solución:

Nos piden graficar la posición contra el tiempo para analizar el movimiento del sistema. Para ello utilizamos los headers proporcionados en clase y el resultado es escrito en un archivo .dat.


![image](https://github.com/javdlgomez/FisicaComputacional/blob/main/Tarea%232/imagenes/ej5-10eps1.png)

Aquí tenemos un problema ya que el valor proporcionado de epsilón no produce una solución adecuada del sistema, para ello tomamos uno de orden 10-5.

![image](https://github.com/javdlgomez/FisicaComputacional/blob/main/Tarea%232/imagenes/ej5-10eps1e-5.png)


Ahora ya es posible apreciar el comportamiento oscilatorio de la solución en nuestra gráfica de posición contra el tiempo.
    
    
    
    
## Ejercicio 5.10:

Nos piden resolver el péndulo simple para cualquier amplitud utilizando el método adaptativo de RK4, para ello debemos resolver la siguiente ecuación diferencial


$$ \theta '' = -sin(\theta)


con las siguientes condiciones iniciales:


$$  \theta'(0) = 1,2,3,4,5,6 \ \ \ \theta''(0) = 0 $$



Código:

    Adjunto únicamente el archivo cpp ya que los headers son equivalentes a los vistos en clase

    //============================================
    //
    // Metodo de RK4 oscilador
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
    void oscilador( const double *y, const double t, double *dydt );





    int main()
    {
      // Datos iniciales
      const double t0 = -10.0;
      const double h = .5;
      const int N = 40; // numero de iteraciones
      const int out_cada = 1; // output cada out_cada iteraciones
      const int n_ec = 2; // numero de ecuaciones



      // reservar espacio para y
      double *y       = new double[ n_ec ];
      double *y_nueva = new double[ n_ec ];

      // inicializar cada variable segun las condiciones iniciales
      y[0]  = sqrt(6/4);
      y[1] = 0.0;


      // puntero a la funcion "derivada"
      void (*derivada)( const double *, const double, double * );
      derivada = oscilador;


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


    void oscilador( const double *y, const double t, double *dydt )
    {
      const double g = 9.8; // aceleracion gravedad
      const double m = 0.01; // masa
      const double k = .0001; // constante de friccion

      dydt[0] = y[1];
      dydt[1] = -sin(y[0]);
    }






    void derivada( const double y, const double t, double &dydt )
    {
      const double g = 9.8; // aceleracion gravedad
      const double m = .01; // masa
      const double k = .0001; // constante de friccion

      dydt = -sin(y);
    }
    
    
    
  Solución: 
  
  
Nos piden graficar los diagramas de fase para los niveles de energía pedidos, el libro nos pide considerar el comportamiento del sistema para distintas energías. Teorícamente esperamos que haya un cambio en el comportamiento del sistema cuando la energía cinéctica supera cierto margen. Pero aquí es donde nos encontramos con uno de los problemas de utilizar métodos numéricos y del por qué es importante la comuniación entre la teoría y la práctica.

Nuestro sistema no presentó el cambio para producir órbitas abiertas después de E=1, estuve buscando cuál podría ser la razón de esto ya que cambiando algunos parámetros pude encontrar comportamientos distintos, ya que no sucedió lo mismo utilizando nuestros scripts de Euler. Al parecer para obtener el resultado esperado utilizando este script es una combinación de colocar correctamente los parámetros en los pasos para que nuestro programa pueda diferenciar el seno de la aproximación en ángulos pequeños. En este caso no llegamos al resultado deseado pero descartamos un camino posible para tomar y que sería recomendado utilizar el método de RK4 adaptativo.

Para los parámetros dados que se encuentran anteriormente en el código nuestro resultado de diagramas de fase es el siguiente:


![image](https://github.com/javdlgomez/FisicaComputacional/blob/main/Tarea%232/imagenes/ejercicio5-12.PNG)

Este comportamiento es idéntico al esperado por la aproximación en movimiento armónico simple donde las órbitas nunca se abren. 
