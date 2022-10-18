# Tarea#3

Solución de el problema de 3 cuerpos sistema Tierra-Luna-Cohete utilizando el método de Runge-Kutta a orden 4 en c++ y adaptando las condiciones iniciales vistas en clase.

## Ejercicio 5.16: 

Nos piden resolver el problema para el caso de la masa Lunar = 0, lo que se reduce al problema de 2 cuerpos.





Código cpp: 


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
      const int N = 100000; // numero de iteraciones
      const int out_cada = 100; // output cada out_cada iteraciones
      const int n_ec = 12; // numero de ecuaciones
      //estos parametros nos permiten jugar con las condiciones iniciales de la posicion del cohete
      const int m = 3;
      const double theta_init = m*M_PI_4;

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
      y[10] = 650.0; //cohete vx
      y[11] =650.0; //cohete vy


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



    void movGrav2D_3cuerpos(  const double *y, const double t, double *dydt )
    {
      const double m1 = 5.97e24; //masa tierra
      const double m2 = 0; //masa luna es 0 para este problema
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



Código gp: 

         #graficador de utf8
        set terminal png notransparent
        set output "ej5-16-1vc.png"




        set title "3 Cuerpos no Luna poca vel inicial"
        set xlabel "x"
        set ylabel "y"

        plot "datosej5-16pc.txt" using 2:3 w lp lw 1 lc  "black"   title "Tierra", "datosej5-16pc.txt" using 6:7 w lp lw 1 lc "red" title "Cohete"



  
Solución:



Ajustamos el método de 3 cuerpos realizado en clase con los parámetros que nos piden y realizamos un output en consola tipo #ej.txt en utf8. Luego utilizamos gnuplot con las especificaciones anteriormente proporcionadas de nuestro graficador para obtener el comportamiento de la órbita del sistema.



![image](https://user-images.githubusercontent.com/100542213/196535874-555a5c39-4bc8-40a6-b3bb-ff6822d294fb.png)

Este es el comportamiento para una velocidad inicial orbital.


![image](https://user-images.githubusercontent.com/100542213/196535916-f7249fc1-0ba4-4780-a075-ff8aa21129c6.png)


Este es el comportamiento para una velocidad menor a la velocidad de escape del punto inicial.




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





![image](https://user-images.githubusercontent.com/100542213/195241516-adbb7ce7-7ac8-4a0d-a005-4a46720a420a.png)



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


![image](https://user-images.githubusercontent.com/100542213/195241547-1476dbae-7e8d-4af1-95a3-583a3c362173.png)


Aquí tenemos un problema ya que el valor proporcionado de epsilón no produce una solución adecuada del sistema, para ello tomamos uno de orden 10-5.

![image](https://user-images.githubusercontent.com/100542213/195241566-ef28c07c-a3fc-4455-867b-5b9018b99e67.png)



Ahora ya es posible apreciar el comportamiento oscilatorio de la solución en nuestra gráfica de posición contra el tiempo.
    
    
    
    
## Ejercicio 5.12:

Nos piden resolver el péndulo simple para cualquier amplitud utilizando el método adaptativo de RK4, para ello debemos resolver la siguiente ecuación diferencial



$$ \theta '' = -sin(\theta)$$


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


![image](https://user-images.githubusercontent.com/100542213/195241900-a2593172-2b20-4a06-8fa8-6d57862c7e89.png)


Este comportamiento es idéntico al esperado por la aproximación en movimiento armónico simple donde las órbitas nunca se abren. 


Jugando con los parámetros logramos reproducir el corte en E=1 en donde las elipses se empiezan a cortar hasta abrirse:


![image](https://user-images.githubusercontent.com/100542213/195241918-b510a36f-8c6b-47cc-8e5f-93fd2812bf64.png)




## Ejercicio 5.13: 

Nos piden animar el péndulo simple, para ello vamos a tomar los resultados de la curva de E=.25



Código: 

    #Generador de imagenes para la animacion




    #Si no hacemos esto no cargan las imágenes y solo las tira en consola

    set terminal pngcairo size 1024,1024 enhanced font "Verdana,12"


    #Ponemos nuestros ejes, en teoría están en metros pero así
    #se sobrecarga menos la animacion

    set xlabel 'x'
    set ylabel 'y'

    #Si no hacemos esto se comprime la imagen conforme la gráfica se vuelve horizontal o vertical
    set xrange [-1:1]
    set yrange [-1.5:.5]


    set title ""
    do for [i=0:19] {

        #colocamos la expresion regular para poder usar ffmpeg despues
        set output sprintf( "animacion/animacion%02d.png", i )

        #$2 significa que lo sacamos de la segunda columna de nuestros datos

        #esta es la cabeza del pendulo
        plot 'd5-12Ep25.txt' every ::i::i u (x1=sin($2)):(y1=-cos($2)) linetype rgb "red"  lw 100 title "" , \
         '' every ::i::i u (0.0):(0.0):(x1):(y1) with vectors nohead  title ""
        #aqui hacemos la cuerda del pendulo que va desde el origen hasta el punto del del dataset
        print i
    }


    # Para correrlo necesitamos ffmpeg en path o en la misma carpeta y lo corremos con
    # ffmpeg -r 10 -i animacion$02d.png -pix_fmt yuv420p animacion.mp4

  
Solución:

Generamos en Gnuplot por medio de un ciclo una imagen para cada punto obtenido en el dataset que están en intervalos de .5s, para ello necesitamos establecer los parámetros del output de la imagen. Además estilizamos un poco para que se vea una masa unida a una cuerda y tenemos un poco de cuidado con la generación de los nombres y en donde se encontraran almacenados. Una vez generadas las imágenes utilizamos ffmpeg para realizar una animación con 10 imágenes por segundo del péndulo simple, estos parámetros fueron elegidos por estilización aunque eso sacrifica la fidelidad física del experimento, esto se hizo para ejemplificar el proceso de animación ya que solo se deben obtener datos con mayor finura que avancen al mismo ritmo del tiempo para producir una simulación que intente ejemplicar la realidad.

![image](https://user-images.githubusercontent.com/100542213/195241958-d451d5db-0353-47db-acbc-0fe89937773d.png)



Este es un frame de ejemplo del resultado de la animación.

