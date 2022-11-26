//===============================================
//
// Tierra Sol Luna Gravitacional
// stdc++20
// Debe compilarse con GCC
//===============================================

#define _USE_MATH_DEFINES
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

// Variables globales
#define n_cuerpos (3) // numero de cuerpos
#define G (6.6743e-11)    // constante gravitacional

std::random_device rd;
std::default_random_engine generator(rd()); // rd() random seed
std::uniform_real_distribution<long double> random_menos1_a_1(-1.0, 1.0);

// Parametros
#define n_iteraciones (365L * 1000)  // numero de iteraciones
#define h_step (3600.0L * 0.24L)     // tamaño de paso
#define medio_h_step (h_step * 0.5L) // tamaño de paso / 2
#define out_cada (50)              // output cada out_cada iteraciones

#define distancia_colision (1.0e7L)

long double tiempo = 0.0;

// Estructuras
struct cuerpo {
  long double pX, pY; // posicion en x, y
  long double vX, vY; // velocidad en x, y
  long double aX, aY; // aceleracion en x, y
  long double masa;   // masa de cada particula
};

struct cuerpo cuerpos[n_cuerpos];

// Variables del sistema
long double E, P, L; // energia momento lineal y momento angular del sistema

// Inicializar masas
void init_masa() {
  // nuestra simulación toma todas las masas con el mismo valor
  cuerpos[0].masa = 5.9722e24L;
  cuerpos[1].masa = 7.342e22L;
  cuerpos[2].masa = 1.98847e30L;
  //for (int i = 3; i < n_cuerpos; i++) {
  //  cuerpos[i].masa = 1.0e17L;
  //}
}

// Inicializar posiciones
void init_posicion() {
  // valores aleatorios entre [-1,1]UA en un cuadrado xy.

  cuerpos[0].pX = -7.355e10L;
  cuerpos[0].pY = 0.0L;
  cuerpos[1].pX = -7.31656e10L;//w0w
  cuerpos[1].pY = 0.0L;
  cuerpos[2].pX = 7.355e10L;
  cuerpos[2].pY = 0.0L;
  //for (int i = 0; i < n_cuerpos; i++) {
  //  cuerpos[i].pX = random_menos1_a_1(generator) * 1.5e11;
  //  cuerpos[i].pY = random_menos1_a_1(generator) * 1.5e11;
  //}
}

// Inicializar velocidades
void init_velocidad() {
  // utilizamos la simplificación para una distribución uniforme de masas
  // puntuales a estas se les agrega luego un valor entre [-1/2,1/2] la magnitud
  // del valor obtenido por la distribución de masas.
  

  cuerpos[0].vX = 0.0L;
  cuerpos[0].vY = 30000L;
  cuerpos[1].vX = 0.0L;
  cuerpos[1].vY = 28900L;
  cuerpos[2].vX = 0.0L;
  cuerpos[2].vY = 0.0L;
  
  /*for (int i = 0; i < n_cuerpos; i++) {
    long double r = sqrt(pow(cuerpos[i].pX, 2) + pow(cuerpos[i].pY, 2));
    long double r_inverso = pow(r, -1);
    long double vel_compartida =
        sqrt(G * M_PI * n_cuerpos * cuerpos[i].masa * r) * 3.0e-11;

    cuerpos[i].vX = -vel_compartida * (cuerpos[i].pY * r_inverso + random_menos1_a_1(generator));
    cuerpos[i].vY = vel_compartida * (cuerpos[i].pX * r_inverso + random_menos1_a_1(generator));
    //(cuerpos[i].pX * r_inverso + random_menos1_a_1(generator) / 2);
  }*/
}

void calc_aceleracion() {
  // consideramos solamente la fuerza de la gravedad.

  for (int i = 0; i < n_cuerpos; i++) {
    cuerpos[i].aX = 0;
    cuerpos[i].aY = 0;
    for (int j = 0; j < n_cuerpos; j++) {
      if (i == j)
        continue;
      long double delta_pX = cuerpos[i].pX - cuerpos[j].pX;
      long double delta_pY = cuerpos[i].pY - cuerpos[j].pY;
      long double G_masa_abs_ri_menos_rj_cubo =
          cuerpos[j].masa * pow(pow(delta_pX, 2) + pow(delta_pY, 2), -1.5) * G;
      cuerpos[i].aX -= delta_pX * G_masa_abs_ri_menos_rj_cubo;
      cuerpos[i].aY -= delta_pY * G_masa_abs_ri_menos_rj_cubo;
    }
  }
}

void verificar_colisiones(long double &t) {
  for (int i = 0; i < n_cuerpos; i++) {
    if (cuerpos[i].masa != 0.0L) {
      for (int j = 0; j < i; j++) {
        if (cuerpos[j].masa != 0.0L) {
          long double deltaX = cuerpos[i].pX - cuerpos[j].pX;
          long double deltaY = cuerpos[i].pY - cuerpos[j].pY;
          long double distancia = sqrt(pow(deltaX, 2) + pow(deltaY, 2));
          if (distancia < distancia_colision) {
            long double nueva_masa = cuerpos[i].masa + cuerpos[j].masa;
            cuerpos[i].vX = (cuerpos[i].masa * cuerpos[i].vX +
                             cuerpos[j].masa * cuerpos[j].vX) /
                            nueva_masa;
            cuerpos[i].vY = (cuerpos[i].masa * cuerpos[i].vY +
                             cuerpos[j].masa * cuerpos[j].vY) /
                            nueva_masa;
            cuerpos[i].masa = nueva_masa;

            // particula j sigue la misma trayectoria que particula i pero sin
            // masa
            cuerpos[j].pX = (cuerpos[i].pX + cuerpos[j].pX) / 2;
            cuerpos[j].pY = (cuerpos[i].pY + cuerpos[j].pY) / 2;
            cuerpos[j].vX = cuerpos[i].vX;
            cuerpos[j].vY = cuerpos[i].vY;
            cuerpos[j].aX = cuerpos[i].aX;
            cuerpos[j].aY = cuerpos[i].aY;
            cuerpos[j].masa = 0.0L;
            std::cout << "Colision " << i << " " << j << " en t = " << t
                      << std::endl;
          }
        }
      }
    }
  }
}

void Energia_Pmomento_Lmomemnto() {

  // obtenemos E, P, L de todo el sistema no analizamos las partículas por
  // separado.

  E = 0.0L;
  P = 0.0L;
  L = 0.0L;
  long double PCM = 0.0L;
  long double mtot = 0.0L;
  long double RCM = 0.0L;
  long double Lloc = 0.0L;
  for (int j = 0; j < n_cuerpos; j++) {
    long double vx2_mas_vy2 = pow(cuerpos[j].vX, 2) + pow(cuerpos[j].vY, 2);
    long double sqrt_vx2_mas_vy2 = sqrt(vx2_mas_vy2);
    E += cuerpos[j].masa * vx2_mas_vy2 * 0.5;
    PCM += cuerpos[j].masa * sqrt_vx2_mas_vy2;
    mtot += cuerpos[j].masa;
    RCM +=
        cuerpos[j].masa * sqrt(pow(cuerpos[j].pX, 2) + pow(cuerpos[j].pY, 2));
    Lloc = RCM * sqrt_vx2_mas_vy2;
  }
  RCM = RCM / mtot;
  E += pow(PCM, 2) * 0.5 / mtot;
  P = PCM;
  L = RCM * PCM + Lloc;
}

void escribirArchivo(std::ofstream &of, std::stringstream &ss) {
  of << ss.str();
}

void salidaSolucion(const long double &t, std::stringstream &ss) {
  // escribimos en un archivo en la columna 0 el valor de t
  // en las siguientes n_cuerpos px y en las consiguientes py
  ss << t;
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << cuerpos[i].pX;
  }
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << cuerpos[i].pY;
  }
  ss << std::endl;
}

void salidaVelocidad(const long double &t, std::stringstream &ss) {
  // realizamos la misma estructura de base de datos para la posicion
  ss << t;
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << cuerpos[i].vX;
  }
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << cuerpos[i].vY;
  }
  ss << std::endl;
}

void salidaAceleracion(const long double &t, std::stringstream &ss) {
  // realizamos la misma estructura de base de datos para la posicion
  ss << t;
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << cuerpos[i].aX;
  }
  for (int i = 0; i < n_cuerpos; ++i) {
    ss << "\t" << cuerpos[i].pY;
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
    ss << "\t" << cuerpos[i].masa;
  }
  ss << std::endl;
}

// Derivadas del sistemas de ecuaciones
void nCuerposGrav(struct cuerpo (&y)[n_cuerpos],
                  struct cuerpo (&dydt)[n_cuerpos]) {
  for (int i = 0; i < n_cuerpos; i++) {
    cuerpos[i] = y[i];
  }
  // Calcular aceleraciones
  calc_aceleracion();

  for (int i = 0; i < n_cuerpos; i++) {
    dydt[i].pX = cuerpos[i].vX;
    dydt[i].pY = cuerpos[i].vY;
    dydt[i].vX = cuerpos[i].aX;
    dydt[i].vY = cuerpos[i].aY;
  }
}

void RK4() {
  struct cuerpo k0[n_cuerpos];
  struct cuerpo k1[n_cuerpos];
  struct cuerpo k2[n_cuerpos];
  struct cuerpo k3[n_cuerpos];

  struct cuerpo z[n_cuerpos];
  for (int i = 0; i < n_cuerpos; i++) {
    z[i] = cuerpos[i];
  }

  nCuerposGrav(cuerpos, k0);

  for (int i = 0; i < n_cuerpos; i++) {
    z[i].pX = cuerpos[i].pX + k0[i].pX * h_step;
    z[i].pY = cuerpos[i].pY + k0[i].pY * h_step;
    z[i].vX = cuerpos[i].vX + k0[i].vX * h_step;
    z[i].vY = cuerpos[i].vY + k0[i].vY * h_step;
  }

  nCuerposGrav(z, k1);

  for (int i = 0; i < n_cuerpos; i++) {
    z[i].pX = cuerpos[i].pX + k1[i].pX * h_step;
    z[i].pY = cuerpos[i].pY + k1[i].pY * h_step;
    z[i].vX = cuerpos[i].vX + k1[i].vX * h_step;
    z[i].vY = cuerpos[i].vY + k1[i].vY * h_step;
  }

  nCuerposGrav(z, k2);

  for (int i = 0; i < n_cuerpos; i++) {
    z[i].pX = cuerpos[i].pX + k2[i].pX * h_step;
    z[i].pY = cuerpos[i].pY + k2[i].pY * h_step;
    z[i].vX = cuerpos[i].vX + k2[i].vX * h_step;
    z[i].vY = cuerpos[i].vY + k2[i].vY * h_step;
  }

  nCuerposGrav(z, k3);

  for (int i = 0; i < n_cuerpos; i++) {
    cuerpos[i].pX += h_step *
                     (k0[i].pX + k1[i].pX * 2.0 + k2[i].pX * 2.0 + k3[i].pX) *
                     pow(6.0, -1);
    cuerpos[i].pY += h_step *
                     (k0[i].pY + k1[i].pY * 2.0 + k2[i].pY * 2.0 + k3[i].pY) *
                     pow(6.0, -1);
    cuerpos[i].vX += h_step *
                     (k0[i].vX + k1[i].vX * 2.0 + k2[i].vX * 2.0 + k3[i].vX) *
                     pow(6.0, -1);
    cuerpos[i].vY += h_step *
                     (k0[i].vY + k1[i].vY * 2.0 + k2[i].vY * 2.0 + k3[i].vY) *
                     pow(6.0, -1);
  }
}

int main() {
  // Strings de salida
  std::stringstream ss_posicion;
  std::stringstream ss_velocidad;
  std::stringstream ss_aceleracion;
  std::stringstream ss_energia;
  std::stringstream ss_masa;

  // Inicializar masa
  init_masa();

  // Inicializar posicion
  init_posicion();

  // Inicializar velocidad
  init_velocidad();

  // ciclo de iteraciones
  auto start = std::chrono::high_resolution_clock::now();
  for (long k = 0; k <= n_iteraciones; k++) {
    if (k % 1000 == 0) {
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration =
          std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
      std::cout << "iteracion: " << k
                << ", tardo en calcularse: " << duration.count() << " ms"
                << std::endl;
      start = stop;
    }
    // energia
    Energia_Pmomento_Lmomemnto();
    // Verificar colisiones
    verificar_colisiones(tiempo);

    // euler_mejorado( y, n_ecuaciones, tiempo, h_step, y_nueva, derivada );
    RK4();
    // salida
    if (k % out_cada == 0) {
      salidaSolucion(tiempo, ss_posicion);
      salidaVelocidad(tiempo, ss_velocidad);
      salidaAceleracion(tiempo, ss_aceleracion);
      salidaEnergia(tiempo, ss_energia);
      salidaMasa(tiempo, ss_masa);
    }

    // Incrementar el tiempo
    tiempo += h_step;
  }

  // Archivos de salida
  std::ofstream of_posicion("posicion.dat", std::ios::out);
  std::ofstream of_velocidad("velocidad.dat", std::ios::out);
  std::ofstream of_aceleracion("aceleracion.dat", std::ios::out);
  std::ofstream of_energia("energia_momLineal_momAngular.dat", std::ios::out);
  std::ofstream of_masa("masa.dat", std::ios::out);

  escribirArchivo(of_posicion, ss_posicion);
  escribirArchivo(of_velocidad, ss_velocidad);
  escribirArchivo(of_aceleracion, ss_aceleracion);
  escribirArchivo(of_energia, ss_energia);
  escribirArchivo(of_masa, ss_masa);

  return 0;
}
