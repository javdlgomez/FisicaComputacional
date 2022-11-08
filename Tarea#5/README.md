
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

Nos piden realizar una animación de un sistema similar al anterior pero con una función de velocidad inicial dada y comprobar si esa onda viaja hacia la izquierda.

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



Resolvemos la ecuación por el método de diferencias finitas y el resultado lo guardamos en un archivo solucion.dat codigicado en utf8, y los datos de la posición del mismo son graficados por medio de un script gpp para verificar el resultado. 


![image](https://user-images.githubusercontent.com/100542213/197990655-674b4531-b4c4-4a24-8156-59cf2c0e6a1d.png)

![image](https://user-images.githubusercontent.com/100542213/197990686-4be94392-29a9-4d2c-9b56-243e7fe0c65f.png)

Ya que la onda se acerca al origen conforme avanza el tiempo podemos verificar que la velocidad de onda es hacia la izquierda.

## Ejercicio 7.3: 

Nos piden demostrar que se obtiene interfencia constructiva bajo los parámetros alpha y beta del mismo signo. Además nos dan ciertas condiciones iniciales para modelar.



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

    //podiamos definir C como global o añadirla a los métodos
    //ya que la solucion de la eq de onda es de tipo
    //f(x+ct)
    double f_cond_ini( double x );
    double g_cond_ini(double x, double vel);
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
      int Niter = 100; // numero de iteraciones en el tiempo
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
        u[i] = u_vieja[i] + g_cond_ini( x[i],vel ) * dt;


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

    //aquí hice una onda que no está fija en los extremos porque no
    //especifica completamente el el libro.
    //Como el fenómeno se cumple para ambos casos, en el siguiente
    //ejercicio como nos lo piden explícitamente tenemos el caso fijo

    double f_cond_ini( double x )
    {
      double L = 1.0; // longitud de la cuerda
      //return sin(4*2.*M_PI*x);
      //return exp(-100*pow(x-L/2,2));


      // por superposicion podemos hacer todo simultaneamente

      // return exp(-100*(x-0.75)*(x-0.75))+exp(-100*(x-0.25)*(x-0.25)); //constructiva 
      return exp(-100*(x-0.75)*(x-0.75))-exp(-100*(x-0.25)*(x-0.25)); //destructiva 
      //return 0.0;
    }


    double g_cond_ini( double x ,double vel)
    {
      double L = 1.0; // longitud de la cuerda
      //return 10*exp(-100*pow(x-L/2,2));

      // return -200*vel*(x-.75)*exp(-100*((x-.75)*(x-.75)))+200*vel*(x-.25)*exp(-100*((x-.25)*(x-.25))); //constructiva 
      return -200*vel*(x-.75)*exp(-100*((x-.75)*(x-.75)))-200*vel*(x-.25)*exp(-100*((x-.25)*(x-.25))); //destructiva 
      //return 0.0;
    }


    double w_cond_frontera( double t ,double vel)
    {
      //return  exp(-100*((0+vel*t)-0.75)*((0+vel*t)-0.75))+exp(-100*((0-vel*t)-0.25)*((0-vel*t)-0.25)); //constructiva 
      return  exp(-100*((0+vel*t)-0.75)*((0+vel*t)-0.75))-exp(-100*((0-vel*t)-0.25)*((0-vel*t)-0.25)); //destructiva
    }


    double z_cond_frontera( double t ,double vel)
    {

      // return  exp(-100*((1+vel*t)-0.75)*((1+vel*t)-0.75))+exp(-100*((1-vel*t)-0.25)*((1-vel*t)-0.25)); //constructiva
      return  exp(-100*((1+vel*t)-0.75)*((1+vel*t)-0.75))-exp(-100*((1-vel*t)-0.25)*((1-vel*t)-0.25)); //destructiva 
    }


Código gp:



    set yrange [-5:5]
    set xrange [0:1]
    dt=0.0001
    set ylabel "y(m)"
    set xlabel "x(m)"
    #cuidado al momento de graficar ya que cambie los nombres de las bases de datos para que se diferenciaran
    do for [it=0:100] {
        set title sprintf( "t = %f (s)", it*dt )
        plot 'solucion.dat' index it u 2:3 w l  title "Solucion Destructiva 100 Iter"

        pause 0.05
    }


Solución:


Resolvemos la ecuación por el método de diferencias finitas y el resultado lo guardamos en un archivo solucion.dat codigicado en utf8, y los datos de la posición del mismo son graficados por medio de un script gpp para verificar el resultado. Realizamos 2 iteraciones distintas con las condiciones iniciales pedidas (signos iguales/signos distintos).

![image](https://user-images.githubusercontent.com/100542213/198173794-952c4c0b-afd1-430e-b448-103980ce97fd.png)

Interferencia Constructiva alpha y beta del mismo signo.

![image](https://user-images.githubusercontent.com/100542213/198173805-ae3658cb-9a1e-4b80-aac7-9d7055bfa7b3.png)

Interferencia Destructiva alpha y beta de signos distintos.



## Ejercicio 7.4: 

Nos piden simular la ecuación diferencial anterior con extremos fijos para 500 iteraciones y reportar si se observa el fenómeno de inversión de fase.



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

    //podiamos definir C como global o añadirla a los métodos
    //ya que la solucion de la eq de onda es de tipo
    //f(x+ct)
    double f_cond_ini( double x );
    double g_cond_ini(double x, double vel);
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
      int Niter = 500; // numero de iteraciones en el tiempo
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
        u[i] = u_vieja[i] + g_cond_ini( x[i],vel ) * dt;


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


      // por superposicion podemos hacer todo simultaneamente

      return exp(-100*(x-0.75)*(x-0.75))+exp(-100*(x-0.25)*(x-0.25)); //constructiva 
      // return exp(-100*(x-0.75)*(x-0.75))-exp(-100*(x-0.25)*(x-0.25)); //destructiva 
      return 0.0;
    }


    double g_cond_ini( double x ,double vel)
    {
      double L = 1.0; // longitud de la cuerda
      //return 10*exp(-100*pow(x-L/2,2));

      return -200*vel*(x-.75)*exp(-100*((x-.75)*(x-.75)))+200*vel*(x-.25)*exp(-100*((x-.25)*(x-.25))); //constructiva 
      // return -200*vel*(x-.75)*exp(-100*((x-.75)*(x-.75)))-200*vel*(x-.25)*exp(-100*((x-.25)*(x-.25))); //destructiva 
      //return 0.0;
    }


    double w_cond_frontera( double t ,double vel)
    {
      //return  exp(-100*((0+vel*t)-0.75)*((0+vel*t)-0.75))+exp(-100*((0-vel*t)-0.25)*((0-vel*t)-0.25)); //constructiva 
      //return  exp(-100*((0+vel*t)-0.75)*((0+vel*t)-0.75))-exp(-100*((0-vel*t)-0.25)*((0-vel*t)-0.25)); //destructiva
      return 0; //extremos fijos 
    }


    double z_cond_frontera( double t ,double vel)
    {

      //return  exp(-100*((1+vel*t)-0.75)*((1+vel*t)-0.75))+exp(-100*((1-vel*t)-0.25)*((1-vel*t)-0.25)); //constructiva
      //return  exp(-100*((1+vel*t)-0.75)*((1+vel*t)-0.75))-exp(-100*((1-vel*t)-0.25)*((1-vel*t)-0.25)); //destructiva 
      return 0; //extremos fijos
    }


Código gp:



    set yrange [-5:5]
    set xrange [0:1]
    dt=0.0001
    set ylabel "y(m)"
    set xlabel "x(m)"
    #cuidado al momento de graficar ya que cambie los nombres de las bases de datos para que se diferenciaran
    do for [it=0:500] {
        set title sprintf( "t = %f (s)", it*dt )
        plot 'solucion.dat' index it u 2:3 w l  title "Solucion Constructiva 500 Iter"

        pause .001
    }

Solución:


Resolvemos la ecuación por el método de diferencias finitas y el resultado lo guardamos en un archivo solucion.dat codigicado en utf8, y los datos de la posición del mismo son graficados por medio de un script gpp para verificar el resultado. Para reproducir el fenómeno necesitamos colocar ciertas condiciones como nulas, el fenómeno ocurre para los casos constructivo y destructivo pero adjuntamos solamente los resultados del primero.

![image](https://user-images.githubusercontent.com/100542213/198174373-e6088ad3-d951-4296-93a1-23ebf408cc23.png)

Inversión de Fase caso constructivo.
