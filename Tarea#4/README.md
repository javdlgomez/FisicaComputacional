# Tarea#4

Solución distintos casos de la ecuación de onda con el método de diferencias finitas.

## Ejercicio 7.1: 

Nos piden resolver la función de onda con extremos fijos con las condiciones iniciales proporcionadas por el libro.





Código cpp: 


    //===================================================
    //
    // Ecuacion de onda con diferencias finitas
    //
    //===================================================


    #include <iostream>
    #include <fstream>
    #include <cmath>

    using namespace std;

    double f_cond_ini( double x );
    double g_cond_ini( double x );
    double w_cond_frontera( double t ,double vel);
    double z_cond_frontera( double t ,double vel);
    void output( ostream &of, double *u, double *x, double t, int N );


    int main()
    {
      int N = 100; //numero de puntos en x
      //este aparece de que h es 1cm por lo tanto necesitamos 100 puntos de grid

      int out_cada = 1; //output cada no. de iteraciones
      double L = 1.0; //longitud del dominio en x
      double dx = L/N;
      double vel = sqrt(10/(.001/1)); // velocidad de la onda
      //para esta nos dan la tensión y densidad lineal de masa
      double dt = 0.0001;  //cond inicial dada
      double alfa = dt*vel/dx; //despejamos la ecuación
      //este es la raiz cuadrada del epsilon del libro

      int Niter = 200; // numero de iteraciones en el tiempo
      //aquí tomé 200 para que la animación se viera lo suficiente

      double tiempo = 0.0; // lleva la cuenta del tiempo
      ofstream outfile;
      outfile.open( "solucion.dat", ios::out );

      // variables para u
      double *u_nueva = new double[N+1]; // u_{i,j+1}
      double *u       = new double[N+1]; // u_{i,j}
      double *u_vieja = new double[N+1]; // u_{i,j-1}
      double *x       = new double[N+1]; // coordenada x


      // coordenada x
      for( int i=0; i<N+1; i++ )
        x[i] = i*dx;

      // condiciones iniciales u_{i0}
      for( int i=0; i<N+1; i++ )
        u_vieja[i] = f_cond_ini( x[i] );

      // condiciones iniciales u_{i1}
      for( int i=0; i<N+1; i++ )
        u[i] = u_vieja[i] + g_cond_ini( x[i] ) * dt;


      // condicion de frontera
      u[0] = w_cond_frontera( 0.0 ,vel);
      u[N] = z_cond_frontera( 0.0 ,vel);


      tiempo += dt;

      // ciclo principal
      for( int j=0; j<=Niter; j++ ){
        for( int i=1; i<N; i++ )
          u_nueva[i] = 2.*(1.-alfa*alfa) * u[i] + alfa*alfa*(u[i-1] + u[i+1]) - u_vieja[i];

        // condicion de frontera
        u_nueva[0] = w_cond_frontera( tiempo + dt ,vel);
        u_nueva[N] = z_cond_frontera( tiempo + dt ,vel);

        // cambiar instantes de tiempo
        for(int i=0; i<N+1; i++ ){
          u_vieja[i] = u[i];
          u[i]       = u_nueva[i];
        }

        tiempo += dt;

        // output
        if ( j % out_cada == 0 )
          output( outfile, u, x, tiempo, N );

      }



      return 0;
    }



    void output( ostream &of, double *u, double *x, double t, int N )
    {
      for( int i=0; i<N+1; i++ )
        of << t << "\t" << x[i] << "\t" << u[i] << endl;

      of << endl << endl;
    }



    double f_cond_ini( double x )
    {
      double L = 1.0; // longitud de la cuerda
      //return sin(4*2.*M_PI*x);
      //return exp(-100*pow(x-L/2,2));
      return 0.0;
    }


    double g_cond_ini( double x )
    {
      double L = 1.0; // longitud de la cuerda
      //return 10*exp(-100*pow(x-L/2,2));
      return 0.0;
    }


    double w_cond_frontera( double t ,double vel)
    {
      return exp(-100*((0-vel*t)-.5)*((0-vel*t)-.5));
    }


    double z_cond_frontera( double t ,double vel)
    {
      return exp(-100*((1-vel*t)-.5)*((1-vel*t)-.5));
    }




Código gp: 

    #graficador tipo animación de una onda

    set yrange [-2.5:2.5]
    dt=0.0001
    set ylabel "y(m)"
    set xlabel "x(m)"
    do for [it=0:200] {
        set title sprintf( "t = %f (s)", it*dt )
        plot 'solucion.dat' index it u 2:3 w l  title "Solucion 200 Iter"

        pause 0.05
    }


  
Solución:



Para resolver este ejercicio debemos interpretar los datos proporcionados por el libro en nuestra definición de código. Para ello tomamos h como la medida de un paso, esto quiere decir que $N = Lh = 100$. De la misma manera encontramos los demás elementos pedidos con las siguientes relaciones:

 $$ dx = \frac{L}{N}, \ \
 c = \sqrt{T/\mu}, \ \  
 \alpha = dt\frac{c}{dx}
 $$

Resolvemos la ecuación por el método de diferencias finitas y el resultado lo guardamos en un archivo solucion.dat codigicado en utf8, y los datos de la posición del mismo son graficados por medio de un script gpp para verificar el resultado.


![image](https://user-images.githubusercontent.com/100542213/197675150-fe30b10a-ed51-4feb-8167-a4836f4dcb6e.png)

Onda resultante de resolver la ecuación con los parámetros indicados, se tomaron 200 iteraciones para poder apreciar mejor el comportamiento.




## Ejercicio 7.2:

Nos piden realizar una animación de un sistema similar al anterior pero con una función de velocidad inicial dada.

Código cpp:


    //===================================================
    //
    // Ecuacion de onda con diferencias finitas
    //
    //===================================================


    #include <iostream>
    #include <fstream>
    #include <cmath>

    using namespace std;

    double f_cond_ini( double x );
    double g_cond_ini( double x , double vel);
    double w_cond_frontera( double t ,double vel);
    double z_cond_frontera( double t ,double vel);
    void output( ostream &of, double *u, double *x, double t, int N );


    int main()
    {
      int N = 100; //numero de puntos en x
      int out_cada = 1; //output cada no. de iteraciones
      double L = 1.0; //longitud del dominio en x
      double dx = L/N;
      double vel = sqrt(10/(.001/1)); // velocidad de la onda
      double dt = 0.0001;  
      double alfa = dt*vel/dx;
      int Niter = 20; // numero de iteraciones en el tiempo
     //aquí si nos piden que sea 20 y se entiende el comportamiento con esas 

      double tiempo = 0.0; // lleva la cuenta del tiempo
      ofstream outfile;
      outfile.open( "solucion.dat", ios::out );

      // variables para u
      double *u_nueva = new double[N+1]; // u_{i,j+1}
      double *u       = new double[N+1]; // u_{i,j}
      double *u_vieja = new double[N+1]; // u_{i,j-1}
      double *x       = new double[N+1]; // coordenada x


      // coordenada x
      for( int i=0; i<N+1; i++ )
        x[i] = i*dx;

      // condiciones iniciales u_{i0}
      for( int i=0; i<N+1; i++ )
        u_vieja[i] = f_cond_ini( x[i] );

      // condiciones iniciales u_{i1}
      for( int i=0; i<N+1; i++ )
        u[i] = u_vieja[i] + g_cond_ini( x[i] ,vel) * dt;


      // condicion de frontera
      u[0] = w_cond_frontera( 0.0 ,vel);
      u[N] = z_cond_frontera( 0.0 ,vel);


      tiempo += dt;

      // ciclo principal
      for( int j=0; j<=Niter; j++ ){
        for( int i=1; i<N; i++ )
          u_nueva[i] = 2.*(1.-alfa*alfa) * u[i] + alfa*alfa*(u[i-1] + u[i+1]) - u_vieja[i];

        // condicion de frontera
        u_nueva[0] = w_cond_frontera( tiempo + dt ,vel);
        u_nueva[N] = z_cond_frontera( tiempo + dt ,vel);

        // cambiar instantes de tiempo
        for(int i=0; i<N+1; i++ ){
          u_vieja[i] = u[i];
          u[i]       = u_nueva[i];
        }

        tiempo += dt;

        // output
        if ( j % out_cada == 0 )
          output( outfile, u, x, tiempo, N );

      }



      return 0;
    }



    void output( ostream &of, double *u, double *x, double t, int N )
    {
      for( int i=0; i<N+1; i++ )
        of << t << "\t" << x[i] << "\t" << u[i] << endl;

      of << endl << endl;
    }



    double f_cond_ini( double x )
    {
      double L = 1.0; // longitud de la cuerda
      //return sin(4*2.*M_PI*x);
      //return exp(-100*pow(x-L/2,2));
      return 0.0;
    }


    double g_cond_ini( double x ,double vel)
    {
      double L = 1.0; // longitud de la cuerda
      return -200*vel*(x-.5)*exp(-100*((x-0)-.5)*((x-0)-.5));
    }


    double w_cond_frontera( double t ,double vel)
    {
      return exp(-100*((0+vel*t)-.5)*((0+vel*t)-.5));
    }


    double z_cond_frontera( double t ,double vel)
    {
      return exp(-100*((1+vel*t)-.5)*((1+vel*t)-.5));
    }

Código gp:



    set yrange [-2.5:2.5]
    set xrange [0:.5]
    dt=0.0001
    set ylabel "y(m)"
    set xlabel "x(m)"

    do for [it=0:20] {
        set title sprintf( "t = %f (s)", it*dt )
        plot 'solucion.dat' index it u 2:3 w l  title "Solucion 20 Iter"

        pause .5
        #como son menos alargué un poco la animación
    }

 
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
    // Metodo de RK4 para movimiento
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



Frame #11 de la animación:


![image](https://user-images.githubusercontent.com/100542213/196541499-e132bae2-a0b9-489d-8383-db07930109d7.png)

    
