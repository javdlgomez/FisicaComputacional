
# Tarea#4

Solución de la ecucación de Calor en una red 2D.



## Ejercicio 7.7: 

Nos piden resolver y luego como extra modelar la solución de las isotermas de la ecuación de calor estacionaria en 2D por medio del método de SOR. Esto en un cuadrado de 1m de largo con una red de 9x9, para las condiciones de frontera nos piden una temperatura de 100 en los bordes izquierdo e inferior y una temperatura de 0 en los bordes derecho y superior.



Código cpp: 

    //============================================
    //
    // Ecuacion de calor en estado estacionario
    // por el metodo SOR
    //
    //=============================================


    #include <cmath>
    #include <iostream>
    #include <fstream>

    using namespace std;



    void output( ostream &of, double **u, double *x, double *y, int Nx, int Ny );
    double p1( double x );
    double p2( double x );
    double q1( double y );
    double q2( double y );




    int main()
    {
      int ITMAX = 10000; // maximo numero de iteraciones
      double eps = 1e-8; // tolerancia de error
      int Nx = 9;
      int Ny = 9;
      double Lx = 1.0;
      double Ly = 1.0;
      double dx = Lx/Nx;
      double dy = dx;
      double alpha = 0.2; // factor para acelerar convergencia en SOR
      ofstream of( "solucion.dat", ios::out);


      // reservar memoria
      double *x        = new double[Nx+1];
      double *y        = new double[Nx+1];
      double **u       = new double*[Nx+1];
      double **u_nueva = new double*[Nx+1];

      for( int i=0; i<Nx+1; i++ ){
        u_nueva[i] = new double[Ny+1];
        u[i]       = new double[Ny+1];
      }


      // coordenadas
      for( int i=0; i<Nx+1; i++ ) x[i] = i*dx;
      for( int j=0; j<Ny+1; j++ ) y[j] = j*dy;


      // inicializar temperatura
      for( int i=0; i<Nx+1; i++ ){
        for( int j=0; j<Ny+1; j++ ){
          u[i][j]       = 0.0;
          u_nueva[i][j] = 0.0;
        }
      }


      // condiciones de frontera
      // lado inferior
      for( int i=0; i<Nx+1; i++ ) u[i][0]  = p1( x[i] );

      // lado superior
      for( int i=0; i<Nx+1; i++ ) u[i][Ny] = p2( x[i] );

      // lado izquierdo
      for( int j=0; j<Ny+1; j++ ) u[0][j]  = q1( y[j] );

      // lado derecho
      for( int j=0; j<Ny+1; j++ ) u[Nx][j] = q2( y[j] );


      // ciclo principal de SOR
      bool seguimos = true; // condicion de salida
      int k = 0; // numero de iteraciones


      cout << "Iniciando SOR" << endl;

      while( seguimos ){
        if ( k > ITMAX ){
          cerr << "Se alcanzó el número máximo de iteraciones para SOR" << endl;
          exit(1);
        }

        seguimos = false;

        for( int i=1; i<Nx; i++ ){
          for( int j=1; j<Ny; j++ ){
        u_nueva[i][j] = 0.25 * ( u[i+1][j] + u[i-1][j] + u[i][j-1] + u[i][j+1] );

        // verificamos si seguimos o no
        if ( fabs( u_nueva[i][j] - u[i][j] ) > eps )
          seguimos = true;

        // cambiamos iteraciones
        u[i][j] = u_nueva[i][j] + alpha*( u_nueva[i][j] - u[i][j] );	
          }
        }

        // terminamos la iteracion
        k++;
      }


      cout << "SOR finalizó en " << k << " iteraciones" << endl;

      // escribir solucion
      output( of, u, x, y, Nx, Ny );


      return 0;
    }




    void output( ostream &of, double **u, double *x, double *y, int Nx, int Ny )
    {
      for( int i=0; i<Nx+1; i++ ){
        for( int j=0; j<Ny+1; j++ )
          of << x[i] << "\t" << y[j] << "\t" << u[i][j] << endl;

        of << endl;
      }
    }




    double p1( double x )
    {
      return 100;
    }

    double p2( double x )
    {
      return 0.0;
    }


    double q1( double y )
    {
      return 100;
    }

    double q2( double y)
    {
      return 0.0;
    }

    
Código gp:




Código gp:



       set terminal qt
       set yrange [0:1]
       set xrange [0:1]
       set ylabel "y(m)"
       set xlabel "x(m)"


       sp  "solucion.dat" title "Mapa de Calor Solucion" with image


       save "hm.gp"



Solución:


Resolvemos la ecuación por el método de diferencias finitas y el resultado lo guardamos en un archivo solucion.dat codigicado en utf8, y los datos de la posición del mismo son graficados por medio de un script gpp para verificar el resultado. Realizamos 2 iteraciones distintas con las condiciones iniciales pedidas (signos iguales/signos distintos).

![image](https://user-images.githubusercontent.com/100542213/198173794-952c4c0b-afd1-430e-b448-103980ce97fd.png)

Interferencia Constructiva alpha y beta del mismo signo.

![image](https://user-images.githubusercontent.com/100542213/198173805-ae3658cb-9a1e-4b80-aac7-9d7055bfa7b3.png)

Interferencia Destructiva alpha y beta de signos distintos.




