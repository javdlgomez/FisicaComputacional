//===============================================
//
// Proyecto problema de n cuerpos gravitacional
// stcd++20
// Debe compilarse con GCC
//===============================================

#define _USE_MATH_DEFINES
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
// #include <sstream>


using namespace std;

// Variables globales
#define n_cuerpos (100)      // numero de cuerpos
#define n_ec (n_cuerpos * 4) // numero de ecuaciones
#define G (6.6743e-11)       // constante gravitacional

std::random_device rd;
std::default_random_engine generator(rd()); // rd() random seed
std::uniform_real_distribution<long double> random_menos1_a_1(-1.0, 1.0);

// Variables globales
long double *xp, *yp;                     // posicion en x, y
long double *vx, *vy;                     // velocidad en x, y
long double ax[n_cuerpos], ay[n_cuerpos]; // aceleracion en x, y
long double masa[n_cuerpos];              // masa de cada particula
long double E, P, L; // energia momento lineal y momento angular del sistema

// Inicializar masas
void init_masa() {
  // nuestra simulación toma todas las masas con el mismo valor
  for (int i = 0; i < n_cuerpos; i++) {
    masa[i] = 10e18;
  }
}

// Inicializar posiciones
void init_posicion() {
  // valores aleatorios entre [-1,1]UA en un cuadrado xy.

  for (int i = 0; i < n_cuerpos; i++) {
    xp[i] = random_menos1_a_1(generator) * 1.5e11;
    yp[i] = random_menos1_a_1(generator) * 1.5e11;
  }
}

// Inicializar velocidades
void init_velocidad() {
  // utilizamos la simplificación para una distribución uniforme de masas
  // puntuales a estas se les agrega luego un valor entre [-1/2,1/2] la magnitud
  // del valor obtenido por la distribución de masas.

  for (int i = 0; i < n_cuerpos; i++) {
    long double r = sqrt(pow(xp[i], 2) + pow(yp[i], 2));
    long double r_inverso = pow(r, -1);
    long double vel_compartida =
        sqrt(G * M_PI * n_cuerpos * masa[i] * r) * 3.0e-11;

    vx[i] =
        vel_compartida * (-1*yp[i] * r_inverso + random_menos1_a_1(generator));
    vy[i] = vel_compartida * (xp[i] * r_inverso + random_menos1_a_1(generator));
  }
}

void calc_aceleracion() {
  // consideramos solamente la fuerza de la gravedad.

  for (int i = 0; i < n_cuerpos; i++) {
    ax[i] = 0;
    ay[i] = 0;
    for (int j = 0; j < n_cuerpos; j++) {
      if (i == j)
        continue;
      long double delta_pX = xp[i] - xp[j];
      long double delta_pY = yp[i] - yp[j];
      long double G_masa_abs_ri_menos_rj_cubo =
          masa[j] * pow(pow(delta_pX, 2) + pow(delta_pY, 2), -1.5) * G;
      ax[i] -= delta_pX * G_masa_abs_ri_menos_rj_cubo;
      ay[i] -= delta_pY * G_masa_abs_ri_menos_rj_cubo;
    }
  }
}

void verificar_colisiones(long double t) {
  long double dist_lim = 1e8;

  for (int i = 0; i < n_cuerpos; i++) {
    if (masa[i] != 0.0) {
      for (int j = 0; j < i; j++) {
        if (masa[j] != 0.0) {
          long double dx = xp[i] - xp[j];
          long double dy = yp[i] - yp[j];
          long double distancia = sqrt(pow(dx, 2) + pow(dy, 2));
          if (distancia < dist_lim) {
            long double nueva_masa = masa[i] + masa[j];
            vx[i] = (masa[i] * vx[i] + masa[j] * vx[j]) / nueva_masa;
            vy[i] = (masa[i] * vy[i] + masa[j] * vy[j]) / nueva_masa;
            masa[i] = nueva_masa;

            // particula j sigue la misma trayectoria que particula i pero sin
            // masa
            xp[j] = (xp[i] + xp[j]) / 2;
            yp[j] = (yp[i] + yp[j]) / 2;
            vx[j] = vx[i];
            vy[j] = vy[i];
            ax[j] = ax[i];
            ay[j] = ay[i];
            masa[j] = 0.0;
            cout << "Colision " << i << " " << j << " en t = " << t << endl;
          }
        }
      }
    }
  }
}

void Energia_Pmomento_Lmomemnto() {

  // obtenemos E, P, L de todo el sistema no analizamos las partículas por
  // separado.

  E = 0.0;
  P = 0.0;
  L = 0.0;
  long double PCM = 0.0;
  long double mtot = 0.0;
  long double RCM = 0.0;
  long double Lloc = 0.0;
  for (int j = 0; j < n_cuerpos; j++) {
    long double vx2_mas_vy2 = pow(vx[j], 2) + pow(vy[j], 2);
    long double sqrt_vx2_mas_vy2 = sqrt(vx2_mas_vy2);
    E += masa[j] * vx2_mas_vy2 * 0.5;
    PCM += masa[j] * sqrt_vx2_mas_vy2;
    mtot += masa[j];
    RCM += masa[j] * sqrt(pow(xp[j], 2) + pow(yp[j], 2));
    Lloc = RCM * sqrt_vx2_mas_vy2;
  }
  RCM = RCM / mtot;
  E += pow(PCM, 2) * 0.5 / mtot;
  P = PCM;
  L = RCM * PCM + Lloc;
}

// Derivadas del sistemas de ecuaciones
void nCuerposGrav(long double y[n_ec], const long double t,
                  long double dydt[n_ec]) {
  xp = y;
  yp = y + n_cuerpos;
  vx = y + 2 * n_cuerpos;
  vy = y + 3 * n_cuerpos;
  // Calcular aceleraciones
  calc_aceleracion();

  for (int i = 0; i < n_cuerpos; i++) {
    dydt[i] = vx[i];
    dydt[i + n_cuerpos] = vy[i];
    dydt[i + 2 * n_cuerpos] = ax[i];
    dydt[i + 3 * n_cuerpos] = ay[i];
  }
}

void RK4(long double *y, const long double t, const long double h,
         long double *y_imas1,
         void (*derivada)(long double *, const long double, long double *)) {
  long double k0[n_ec];
  long double k1[n_ec];
  long double k2[n_ec];
  long double k3[n_ec];
  long double z[n_ec];

  (*derivada)(y, t, k0);

  for (int i = 0; i < n_ec; i++)
    z[i] = y[i] + 0.5 * k0[i] * h;

  (*derivada)(z, t + 0.5 * h, k1);

  for (int i = 0; i < n_ec; i++)
    z[i] = y[i] + 0.5 * k1[i] * h;

  (*derivada)(z, t + 0.5 * h, k2);

  for (int i = 0; i < n_ec; i++)
    z[i] = y[i] + k2[i] * h;

  (*derivada)(z, t + h, k3);

  for (int i = 0; i < n_ec; i++)
    y_imas1[i] = y[i] + h / 6.0 * (k0[i] + 2 * k1[i] + 2 * k2[i] + k3[i]);
}

void escribirArchivo(std::ofstream &of, std::stringstream &ss) {
  of << ss.str();
}

void salidaSolucion(const long double &t, std::stringstream &ss) {
  // escribimos en un archivo en la columna 0 el valor de t
  // en las siguientes n_cuerpos px y en las consiguientes py
  ss << t;
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << xp[i];
  }
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << yp[i];
  }
  ss << std::endl;
}

void salidaVelocidad(const long double &t, std::stringstream &ss) {
  // realizamos la misma estructura de base de datos para la posicion
  ss << t;
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << vx[i];
  }
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << vy[i];
  }
  ss << std::endl;
}

void salidaAceleracion(const long double &t, std::stringstream &ss) {
  // realizamos la misma estructura de base de datos para la posicion
  ss << t;
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << ax[i];
  }
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << ay[i];
  }
  ss << std::endl;
}

// energia, momentum lineal y angular
void salidaEnergia(const long double &t, std::stringstream &ss) {
  // escribimos en un archivo el valor del tiempo consiguiente
  // con el valor de la energía, momento lineal y momento angular.
  ss << t << "\t" << E << "\t" << P << "\t" << L << std::endl;
}

void salidaMasa(long double &t, std::stringstream &ss) {
  // escribimos en un archivo el valor del tiempo
  // seguido por el valor de la masa de todos los cuerpos.
  ss << t;
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << masa[i];
  }
  ss << std::endl;
}

int main() {
  // Parametros
  const int Niter = 52 * 5000;           // numero de iteraciones
  const long double h_step = 3600 * 24*7; // tamaño de paso
  const int out_cada = 10000;              // output cada out_cada iteraciones

  // Otras variables
  long double tiempo = 0.0;

  // reservar espacio para y
  long double y[n_ec];
  long double y_nueva[n_ec];

  // reservar espacio para posicion, velocidad, aceleracion y masa
  xp = y;
  yp = y + n_cuerpos;
  vx = y + 2 * n_cuerpos;
  vy = y + 3 * n_cuerpos;

  // puntero a la funcion "derivada"
  void (*derivada)(long double *, const long double, long double *);
  derivada = nCuerposGrav;

  // inicializar y_nueva
  for (int i = 0; i < n_ec; i++)
    y_nueva[i] = 0.0;

  // Inicializar masa
  init_masa();

  // Inicializar posicion
  init_posicion();

  // Inicializar velocidad
  init_velocidad();

  // Strings de salida
  std::stringstream ss_posicion;
  std::stringstream ss_velocidad;
  std::stringstream ss_aceleracion;
  std::stringstream ss_energia;
  std::stringstream ss_masa;

  // ciclo de iteraciones
  auto start = std::chrono::high_resolution_clock::now();
  for (int k = 0; k <= Niter; k++) {
    if (k % 1000 == 0) {
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration =
          std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      std::cout << "iteracion: " << k
                << ", tardo en calcularse: " << duration.count() << " ms"
                << std::endl;
      start = stop;
    }
    xp = y;
    yp = y + n_cuerpos;
    vx = y + 2 * n_cuerpos;
    vy = y + 3 * n_cuerpos;
    // energia
    Energia_Pmomento_Lmomemnto();
    // Verificar colisiones
    verificar_colisiones(tiempo);

    // euler_mejorado( y, n_ec, tiempo, h_step, y_nueva, derivada );
    RK4(y, tiempo, h_step, y_nueva, derivada);
    // salida
    if (k % out_cada == 0) {
      salidaSolucion(tiempo, ss_posicion);
      salidaVelocidad(tiempo, ss_velocidad);
      salidaAceleracion(tiempo, ss_aceleracion);
      salidaEnergia(tiempo, ss_energia);
      salidaMasa(tiempo, ss_masa);
    }

    // Intercambiar valores
    for (int i = 0; i < n_ec; i++)
      y[i] = y_nueva[i];

    // Incrementar el tiempo
    tiempo += h_step;
  }

  // Archivos de salida
  ofstream of_posicion("posicion.dat", ios::out);
  ofstream of_velocidad("velocidad.dat", ios::out);
  ofstream of_aceleracion("aceleracion.dat", ios::out);
  ofstream of_energia("energia_momLineal_momAngular.dat", ios::out);
  ofstream of_masa("masa.dat", ios::out);

  escribirArchivo(of_posicion, ss_posicion);
  escribirArchivo(of_velocidad, ss_velocidad);
  escribirArchivo(of_aceleracion, ss_aceleracion);
  escribirArchivo(of_energia, ss_energia);
  escribirArchivo(of_masa, ss_masa);

  return 0;
}
