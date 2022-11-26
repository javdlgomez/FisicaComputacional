# Proyecto final
## Introducción

El problema de los n cuerpos se refiere de forma usual a la solución del sistema de movimiento de n masas las cuales interactuan entre sí por medio de la fuerza gravitacional en 3D. Este problema no está limitado a esas restricciones y pueden añadirse otros factores o estudiar el mismo comportamiento con otras fuerzas que se comportan de manera similar como la fuerza eléctrica. 
Este sistema deja de tener solución analítica en el caso general para n>2, ya que la cantidad de grados de libertad es menor a la cantidad de ecuaciones disponibles para resolver el sistema. A pesar de esto existen casos particulares de los cuales se conocen las soluciones analíticas para más dimensiones, en especial n=3. La solucíón con condiciones arbitrarias de este sistema para los casos con n>2 puede ser encontrada mediante métodos numéricos, aproximaciones o casos reducidos de geoemtrías particulares.

Nosotros estudiamos el caso particular donde se tienen 100 masas puntuales y simétricas distribuidas aleatoriamente en un cuadrado 2D, afectadas únicamente por la fuerza gravitacional. Las velocidades iniciales se obtienen mediante una simplificación asumiendo una distribución de masa uniforme y se les es agregado un factor aleatorio proporcional a la misma. Y además si 2 de estas masas llegasen a acercarse de forma que la distancia entre ellas fuese menor a un threshold dado, estas participarían en una colisión inelástica perfecta.

Para modelar este sistema se realizó una simulación en c++, utilizando el método de RK4 para resolver el sistema de ecuaciones durante 5000 mil años de evolución. Se realizaron modelos con distintos intervalos de aumento en el tiempo y se encontró la evolución de las posiciones, velocidades, aceleraciones, energía, momento y momento angular total del sistema. 




## Métodos 
### Generalidades:
Para resolver el problema se realizaron dos simulaciones computacionales en c++ bajo el estandar stdc++20. Estas simulaciones realizan la aproximación de la solución del sistema de movimiento empleando el método de RK4, este funciona tomando las condiciones iniciales dadas por nuestra solución particular y obtiene las variables físicas en el siguiente intervalo de tiempo. Este algoritmo se repite tomando como los valores obtenidos las nuevas condiciones iniciales y se realiza hasta que se llegue al tiempo de duración deseado.


### Condiciones iniciales:
Para obtener las posiciones iniciales se utilizó el generador de números pseudo aleatorios por defecto de la libreria random de c++ y se distribuyeron las masas en un cuadrado centrado en el origen de lado 1UA.

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


### Validez de los modelos:

Para verificar la validez de nuestra implementación se realizó una simulación del sistema de 3 cuerpos Tierra Luna Sol utilizando como condiciones iniciales los datos de la geometría y velocidad conocidas de los mismos. Además se miden las variables de Energía, Momento y Momentu angular con respecto del tiempo y al obtener que el comportamiento de las mismas es de acuerdo con la teoría entonces podemos tener cierto grado de seguridad de la fidelidad del mismo.



### Diferencias entre las simulaciones:
La primera implementación se utiliza como base una optimización de el script proporcionado por, este es ligeramente peligroso ya que realiza cambios directamente en memoria sin protección del complilador y está construido con un paradigma de programación funcional.
Esta versión obtiene resultados más precisos, pero conlleva mayor tiempo de ejecución.

La segunda implementación utiliza un enfoque mixto empleando propiedades del OOP e intenta proteger un poco más a posibles errores durante la ejecución del programa. Esta versión diverge con mayor facilidad a la solución esperada bajo las mismas condiciones, pero tiene un menor tiempo de ejecución.



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

## Discusión de Resultados 
## Conclusiones 
## Referencias
