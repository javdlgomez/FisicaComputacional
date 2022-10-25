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
