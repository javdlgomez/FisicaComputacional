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


Este es el comportamiento para una velocidad inicial que no puede vencer la atracción gravitacional de la tierra.


También esperamos otra solución para una velocidad mayor a la orbital en la que el cohete escapa a la tierra.




## Ejercicio 5.17: 

Nos piden encontrar condiciones iniciales en las que el cohete logre acercarse lo más posible a la Luna (a menos de 300km del centro del satélite)

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
      const int N = 20000; // numero de iteraciones
      const int out_cada = 100; // output cada out_cada iteraciones
      const int n_ec = 12; // numero de ecuaciones
      //estos parametros nos permiten jugar con las condiciones iniciales de la posicion del cohete
      const int m = 0;
      const double theta_init = m*M_PI_4;
      const double v0 = 25000; // velocidad inicial

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
      y[10] = v0*cos(theta_init+.04); //cohete vx
      y[11] = v0*sin(theta_init+.04); //cohete vy

      //de esta forma podemos con parametros theta y v0 cambiar las condiciones iniciales


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
      double dx = (y[4]-y[6])*(y[4]-y[6]);
      double dy = (y[5]-y[7])*(y[5]-y[7]);
      double dis = sqrt(dx+dy);


      (*derivada)( y, t, k0 );

      for( int i=0; i<n_ec; i++ )
        z[i] = y[i] + 0.5*k0[i]*h;

      (*derivada)( z, t+0.5*h, k1 );

      for( int i=0; i<n_ec; i++ )
        z[i] = y[i] + 0.5*k1[i]*h;

      (*derivada)( z, t+0.5*h, k2 );

      for( int i=0; i<n_ec; i++ )
        z[i] = y[i] + k2[i]*h;




        ;

      (*derivada)( z, t+h, k3 );

      for( int i=0; i<n_ec; i++ )

       y_imas1[i] = y[i] + h/6.0 * ( k0[i] + 2*k1[i] + 2*k2[i] + k3[i] );
        //if ( pow(((y[4]-y[6]^2+(y[5]-y[7]^2)),.5))){




      delete[] k0;
      delete[] k1;
      delete[] k2;
      delete[] k3;
      delete[] z;
    }



    void movGrav2D_3cuerpos(  const double *y, const double t, double *dydt )
    {
      const double m1 = 5.97e24; //masa tierra
      const double m2 = 7.34e22; //masa luna es 0 para este problema
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


 Código de los gps:
 
 
 XvsY:
 
 
     #graficador de utf8
    set terminal png notransparent
    set output "ej5-17xy.png"




    set title "3 Cuerpos Aterrizaje XvsY "
    set xlabel "x(Mm)"
    set ylabel "y(Mm)"

    plot  "datosej5-17.txt" using ($6/1e6):($7/1e6) w lp lw .5 lc "red" title "Cohete",  "datosej5-17.txt" using ($4/1e6):($5/1e6) w lp lw .5 lc "blue" title "Luna"
    
    
xvst:

     #graficador de utf8
    set terminal png notransparent
    set output "ej5-17.png"




    set title "3 Cuerpos aterrizaje xvst "
    set xlabel "t(s)"
    set ylabel "x (m)"

    plot  "datosej5-17.txt" using 1:6 w lp lw .5 lc "red" title "Cohete",  "datosej5-17.txt" using 1:4 w lp lw .5 lc "blue" title "Luna"
    
    
    
yvst: 

    #graficador de utf8
    set terminal png notransparent
    set output "ej5-17y.png"




    set title "3 Cuerpos aterrizaje yvst"
    set xlabel "t(s)"
    set ylabel "y (m)"

    plot  "datosej5-17.txt" using 1:7 w lp lw .5 lc "red" title "Cohete",  "datosej5-17.txt" using 1:5 w lp lw .5 lc "blue" title "Luna"
 
Solución:

Realizamos el mismo procedimiento para el problema anterior ajustando las condiciones iniciales en términos de los parámetros v0 y theta_init. Existe un desfaz en las trigonométricas para producir el comportamiento de cercanía deseado. Luego el resultado producido en consola es guardado en un archivo #ej.txt y es graficado con los scripts presentados anteriormente en gnuplot.

![image](https://user-images.githubusercontent.com/100542213/196537552-ef60368c-b15c-4adb-9427-c5ca374a8a64.png)


![image](https://user-images.githubusercontent.com/100542213/196537575-4c07edf0-9849-4836-b0c4-adae8393ac69.png)


![image](https://user-images.githubusercontent.com/100542213/196537591-70355689-70aa-4738-b18e-79eb99a5c145.png)



## Ejercicio 5.10: 

Nos piden resolver el problema anterior pero ahora necesitamos que el cohete regrese a la tierra luego de visitar la Luna, los parámetros iniciales para esta situación son los más sensibles de las anteriores. La idea de nuestros parámetros es lanzar al cohete con una velocidad moderada a donde la Luna se moverá en el futuro, de esta manera el cohete seguirá en una órbita ligada al sistema y tendrá la velocidad suficiente para regresar. Luego estos resultados se visualizan en una animación.




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
      const int N = 1100000; // numero de iteraciones
      const int out_cada = 50000; // output cada out_cada iteraciones
      const int n_ec = 12; // numero de ecuaciones
      //estos parametros nos permiten jugar con las condiciones iniciales de la posicion del cohete
      const int m = 0;
      const double theta_init = m*M_PI_4;
      const double v0 = 14000; // velocidad inicial

      // Archivo que guarda la energia total


      // reservar espacio para y
      double *y       = new double[ n_ec ];
      double *y_nueva = new double[ n_ec ];

      // inicializar cada variable segun las condiciones iniciales
      y[0]  = 0.0; //tierra x
      y[1]  = 0.0; // tierra y
      y[2]  = 3.84e8; //luna x
      y[3]  = 0.0; //luna y
      y[4]  = (6800e3)*.8; //cohete x
      y[5]  = (6800e3)*sqrt(1-.8*.8); //cohete y
      y[6]  = 0.0; // tierra vx
      y[7]  = 0.0; //tierra cy
      y[8]  = 0.0;//luna vx
      y[9]  = 1033;//luna vy (promedio de rango de velocidades)
      y[10] =v0*.552; //cohete vx
      y[11] = v0*.53; //cohete vy

      //de esta forma podemos con parametros theta y v0 cambiar las condiciones iniciales


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
      double dx = (y[4]-y[6])*(y[4]-y[6]);
      double dy = (y[5]-y[7])*(y[5]-y[7]);
      double dis = sqrt(dx+dy);


      (*derivada)( y, t, k0 );

      for( int i=0; i<n_ec; i++ )
        z[i] = y[i] + 0.5*k0[i]*h;

      (*derivada)( z, t+0.5*h, k1 );

      for( int i=0; i<n_ec; i++ )
        z[i] = y[i] + 0.5*k1[i]*h;

      (*derivada)( z, t+0.5*h, k2 );

      for( int i=0; i<n_ec; i++ )
        z[i] = y[i] + k2[i]*h;




        ;

      (*derivada)( z, t+h, k3 );

      for( int i=0; i<n_ec; i++ )

       y_imas1[i] = y[i] + h/6.0 * ( k0[i] + 2*k1[i] + 2*k2[i] + k3[i] );
        //if ( pow(((y[4]-y[6]^2+(y[5]-y[7]^2)),.5))){




      delete[] k0;
      delete[] k1;
      delete[] k2;
      delete[] k3;
      delete[] z;
    }



    void movGrav2D_3cuerpos(  const double *y, const double t, double *dydt )
    {
      const double m1 = 5.97e24; //masa tierra
      const double m2 = 7.34e22; //masa luna es 0 para este problema
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


XvsY:

    #graficador de utf8
    set terminal png notransparent
    set output "ej5-18.png"




    set title "3 Cuerpos Regreso Órbita "
    set xlabel "x (m)"
    set ylabel "y (m)"

    plot  "datosej5-18.txt" using 6:7 w lp lw .5 lc "red" title "Cohete",  "datosej5-18.txt" using 4:5 w lp lw .5 lc "blue" title "Luna" ,  "datosej5-18.txt" using 2:3 w lp lw .5 lc "green" title "Tierra"


xvst:


    #graficador de utf8
    set terminal png notransparent
    set output "ej5-18xt.png"




    set title "3 Cuerpos Regreso xvst "
    set xlabel "t(s)"
    set ylabel "x (m)"

    plot  "datosej5-18.txt" using 1:6 w lp lw .5 lc "red" title "Cohete",  "datosej5-18.txt" using 1:4 w lp lw .5 lc "blue" title "Luna","datosej5-18.txt" using 1:2 w lp lw .5 lc "green" title "Tierra"


yvst:


    #graficador de utf8
    set terminal png notransparent
    set output "ej5-18y.png"




    set title "3 Cuerpos Regreso yvst"
    set xlabel "t(s)"
    set ylabel "y (m)"

    plot  "datosej5-18.txt" using 1:7 w lp lw .5 lc "red" title "Cohete",  "datosej5-18.txt" using 1:5 w lp lw .5 lc "blue" title "Luna", "datosej5-18.txt" using 1:3 w lp lw .5 lc "green" title "Tierra"


creador de frames para animación:


    #Generador de imagenes para la animacion




    #Si no hacemos esto no cargan las imágenes y solo las tira en consola

    set terminal pngcairo size 1024,1024 enhanced font "Verdana,12"


    #Ponemos nuestros ejes, en teoría están en metros pero así
    #se sobrecarga menos la animacion

    set xlabel 'x (1e8m)'
    set ylabel 'y (1e8m)'

    #Si no hacemos esto se comprime la imagen conforme la gráfica se vuelve horizontal o vertical
    set xrange [-5:5]
    set yrange [-5:5]




    set title ""
    do for [i=0:46] {

        #colocamos la expresion regular para poder usar ffmpeg despues
        set output sprintf( "animacion/animacion%02d.png", i )

        #$2 significa que lo sacamos de la segunda columna de nuestros datos

        #esta es la cabeza del pendulo
        plot 'datosej5-18.txt' every ::i::i u ($6/1e8):($7/1e8) linetype rgb "black"  lw 5 title "cohete" , 'datosej5-18.txt' every ::i::i u ($2/1e8):($3/1e8) linetype rgb "green"  lw 40 title "tierra", 'datosej5-18.txt' every ::i::i u ($4/1e8):($5/1e8) linetype rgb "blue"  lw 20 title "luna"
        #aqui hacemos la cuerda del pendulo que va desde el origen hasta el punto del del dataset

    }


    # Para correrlo necesitamos ffmpeg en path o en la misma carpeta y lo corremos con
    # ffmpeg -r 10 -i animacion$02d.png -pix_fmt yuv420p animacion.mp4






Solución:

Realizamos el mismo procedimiento de los ejercicios anteriores, en este caso tenemos mucho cuidad con las condiciones iniciales. Luego el resultado producido en consola es guardado en un archivo #ej.txt y es graficado con los scripts presentados anteriormente en gnuplot. Para producir la animación se dividie el intervalo de la simulación en 47 imágenes de XvsY de los 3 cuerpos y por medio de ffmpeg se animan en 1s.

![image](https://user-images.githubusercontent.com/100542213/196541393-2a3a01b1-e958-41bb-95e8-6077d49369c5.png)


![image](https://user-images.githubusercontent.com/100542213/196541435-a33484be-b08d-4830-8819-c67b3e6da86a.png)



![image](https://user-images.githubusercontent.com/100542213/196541452-7c1330cd-e6cb-4448-ac7d-adfbaa282d06.png)



![image](https://user-images.githubusercontent.com/100542213/196541499-e132bae2-a0b9-489d-8383-db07930109d7.png)

    
   


Este es un frame de ejemplo del resultado de la animación.

