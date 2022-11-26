# Proyecto final
## Introducción

El problema de los n cuerpos se refiere de forma usual a la solución del sistema de movimiento de n masas las cuales interactuan entre sí por medio de la fuerza gravitacional en 3D. Este problema no está limitado a esas restricciones y pueden añadirse otros factores o estudiar el mismo comportamiento con otras fuerzas que se comportan de manera similar como la fuerza eléctrica. 
Este sistema deja de tener solución analítica en el caso general para n>2, ya que la cantidad de grados de libertad es menor a la cantidad de ecuaciones disponibles para resolver el sistema. A pesar de esto existen casos particulares de los cuales se conocen las soluciones analíticas para más dimensiones, en especial n=3. La solucíón con condiciones arbitrarias de este sistema para los casos con n>2 puede ser encontrada mediante métodos numéricos, aproximaciones o casos reducidos de geoemtrías particulares.

Nosotros estudiamos el caso particular donde se tienen 100 masas puntuales y simétricas distribuidas aleatoriamente en un cuadrado 2D, afectadas únicamente por la fuerza gravitacional. Las velocidades iniciales se obtienen mediante una simplificación asumiendo una distribución de masa uniforme y se les es agregado un factor aleatorio proporcional a la misma. Y además si 2 de estas masas llegasen a acercarse de forma que la distancia entre ellas fuese menor a un threshold dado, estas participarían en una colisión inelástica perfecta.

Para modelar este sistema se realizó una simulación en c++, utilizando el método de RK4 para resolver el sistema de ecuaciones durante 5000 mil años de evolución. Se realizaron modelos con distintos intervalos de aumento en el tiempo y se encontró la evolución de las posiciones, velocidades, aceleraciones, energía, momento y momento angular total del sistema. 




## Métodos 
### Generalidades:
Para resolver el problema se realizaron dos scripts distintos de simulaciones computacionales en c++ bajo el estandar stdc++20 y fue compilado utilizando flags de optimización como Ofast, se realizaron pruebas para comprobar que no ocurrieran problemas de precisión numérica. Estas simulaciones realizan la aproximación de la solución del sistema de movimiento empleando el método de RK4, este funciona tomando las condiciones iniciales dadas por nuestra solución particular y obtiene las variables físicas en el siguiente intervalo de tiempo. Este algoritmo se repite tomando como los valores obtenidos las nuevas condiciones iniciales y se realiza hasta que se llegue al tiempo de duración deseado.


### Condiciones iniciales:
Para obtener las posiciones iniciales se utilizó el generador de números pseudo aleatorios por defecto de la libreria random de c++ y se distribuyeron las masas en un cuadrado centrado en el origen de lado 2UA.

Para obtener las velocidades iniciales se utilizó la expresión obtenida en :


$$ r_i = \sqrt{x_i^2+y_i^2}  $$

$$   v_{i0} = \frac{\sqrt{G \pi N m_i  r_i}}{L}(-\frac{y_{i0}}{r_i},\frac{x_{i0}{r_i}}) $$


Además se agregó un término aleatorio que oscila entre $\pm \frac{\sqrt{G\pi N m_i r_i}}{L}$.

Para encontrar la aceleraciones iniciales se despejó el sistema obtenido por ley de Newton y se llegó al siguiente resultado como fue dado en  :


$$ |r_i-r_j| = \sqrt{(x_i-x_j)^2+(y_i-y_j)^2} $$

 $$\ddot {{x_i }}= \sum_{i\neq j} - \frac{Gm_j}{|r_i-r_j|^3}(x_i-x_j) $$

$$ \ddot { {y_i }}= \sum_{i\neq j} - \frac{Gm_j}{|r_i-r_j|^3}(y_i-y_j) $$


Finalmente se coloca el valor de todas las masas a $m_i = 10^18 kg$.
Es importante agregar que nuestra simulación se generalizó para aceptar valores distintos de masas puntuales.

### Cálculo de las variables físicas:


Para encontrar la energía total se utiliza la ecuación de la energía para un sistema de partículas:

$$ E = \frac{MR^2}{2}    +\sum \frac{m_iv_i^2}{2}  $$

Donde M es la masa total y R es la posición del centro de masa del sistema.

Para encontrar el momento total se utilizó la ecuación del momento para un sistema de partículas:

$$ P = \sum m_iv_i $$

Para encontrar el momento agnular total se utilizó la ecuación del momento angular para un sistema de partículas:

$$ L = R \times P + \sum r_i \times m_i v_i $$


Para calcular la evolución de las variables de movimiento utilizamos la implementación del algoritmo de RK4 y repetimos el proceso.

### Cálculo de las colisiones:

Se utiliza un threshold de colisión bastante pequeño comparado con el tamaño de la geometría simulación pero basado físicamente en un radio de planeta ligeramente mayor al radio de la tierra $$ th = 1e^7$$. Por este motivo las colisiones son poco comunes y puede existir un traslape en las animaciones donde parece que ocurre una colisión porque el tamaño animado es más grande que el threshold, o sea su tamaño físico dado.
Al ocurrir una colisión se resuelve el sistema de ecuaciones para una colisión elástica perfecta con respecto de las variables físicas finales y se unen los 2 cuerpos sumando sus masas asignándole el valor a una de ellas y descartando la masa de la otra.

### Validez de los modelos:

Para verificar la validez de nuestra implementación se realizó una simulación del sistema de 3 cuerpos Tierra Luna Sol utilizando como condiciones iniciales los datos de la geometría y velocidad conocidas de los mismos. Además se miden las variables de Energía, Momento y Momentum angular con respecto del tiempo y al obtener que el comportamiento constante de las mismas es de acuerdo con la teoría entonces podemos tener cierto grado de seguridad de la fidelidad del mismo.



### Diferencias entre las simulaciones:
La primera implementación se utiliza como base una optimización de el script proporcionado por, este realiza el proceso del método de RK4 directamente en memoria sin protección del complilador y está construido con un paradigma de programación funcional.
Esta versión obtiene resultados más precisos, pero conlleva mayor tiempo de ejecución.

La segunda implementación utiliza un enfoque mixto empleando propiedades del OOP e intenta proteger un poco más a posibles errores durante la ejecución del programa. Aún así versión diverge con mayor facilidad a la solución esperada bajo las mismas condiciones, pero tiene un menor tiempo de ejecución.



## Resultados


### Funciones Simulación Script Original:

#### Posiciones iniciales:


     std::random_device rd;
     std::default_random_engine generator(rd()); // rd() random seed
     std::uniform_real_distribution<long double> random_menos1_a_1(-1.0, 1.0);
     // Inicializar posiciones
     void init_posicion() {
     // valores aleatorios entre [-1,1]UA en un cuadrado xy.

     for (int i = 0; i < n_cuerpos; i++) {
       xp[i] = random_menos1_a_1(generator) * 1.5e11;
       yp[i] = random_menos1_a_1(generator) * 1.5e11;
      }
      }

                                   
#### Velocidades iniciales: 


                                   
                                   
     // Inicializar velocidades
     void init_velocidad() {
     // utilizamos la simplificación para una distribución uniforme de masas
     // puntuales a estas se les agrega luego un valor entre [-1,1] la magnitud
     // del valor obtenido por la distribución de masas.

     for (int i = 0; i < n_cuerpos; i++) {
       long double r = sqrt(pow(xp[i], 2) + pow(yp[i], 2));
       long double r_inverso = pow(r, -1);
       long double vel_compartida =
           sqrt(G * M_PI * n_cuerpos * masa[i] * r) / (3.0e11);

       vx[i] =
           vel_compartida * (-1*yp[i] * r_inverso + random_menos1_a_1(generator));
       vy[i] = vel_compartida * (xp[i] * r_inverso + random_menos1_a_1(generator));
     }
     }                                   

 
#### Cálculo de las aceleraciones: 

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
 
 #### Verificación de Colisiones, se agregó animación: 
 
 
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

 
 #### Cálculo de E, P y L totales:
 
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
     //L tiene un problema y es que no está considerando las posiciones relativas
     // 
    }


#### Salida de Archivos:

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

Se omiten el resto de funciones de escritura por su similitud.
 
 
 ### Funciones Simulación Script Nuevo:
 
 #### Estructura de Cuerpos para la Simulación:
 

      struct cuerpo cuerpos[n_cuerpos];

     // Variables del sistema
     long double E, P, L; // energia momento lineal y momento angular del sistema

     // Inicializar masas
     void init_masa() {
       // nuestra simulación toma todas las masas con el mismo valor
       for (int i = 0; i < n_cuerpos; i++) {
         cuerpos[i].masa = 10.0e18L;
       }
     }



 
 #### Posiciones iniciales:
 

     std::random_device rd;
    std::default_random_engine generator(rd()); // rd() random seed
    std::uniform_real_distribution<long double> random_menos1_a_1(-1.0, 1.0);


     // Inicializar posiciones
    void init_posicion() {
      // valores aleatorios entre [-1,1]UA en un cuadrado xy.

      for (int i = 0; i < n_cuerpos; i++) {
        cuerpos[i].pX = random_menos1_a_1(generator) * 1.5e11;
        cuerpos[i].pY = random_menos1_a_1(generator) * 1.5e11;
      }
    }
    
    
   
 
#### Velocidades iniciales: 
 

     // Inicializar velocidades
    void init_velocidad() {
      // utilizamos la simplificación para una distribución uniforme de masas
      // puntuales a estas se les agrega luego un valor entre [-1/2,1/2] la magnitud
      // del valor obtenido por la distribución de masas.

      for (int i = 0; i < n_cuerpos; i++) {
        long double r = sqrt(pow(cuerpos[i].pX, 2) + pow(cuerpos[i].pY, 2));
        long double r_inverso = pow(r, -1);
        long double vel_compartida =
            sqrt(G * M_PI * n_cuerpos * cuerpos[i].masa * r) /(3e11);

        cuerpos[i].vX = -vel_compartida *
                        (cuerpos[i].pY * r_inverso + random_menos1_a_1(generator));
        cuerpos[i].vY = vel_compartida *
                        (cuerpos[i].pX * r_inverso + random_menos1_a_1(generator));
      }
    }
 
 
 #### Aceleraciones iniciales: 
 
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


#### Verificación de colisiones:


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


## Discusión de Resultados 
## Conclusiones 
## Referencias
