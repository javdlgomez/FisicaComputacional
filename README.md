# FisicaComputacional
ejercicios del curso fisica computacional

# Proyecto final
## Introducción

El problema de los n cuerpos se refiere de forma usual a la solución del sistema de movimiento de n masas las cuales interactuan entre sí por medio de la fuerza gravitacional en 3D. Este problema no está limitado a esas restricciones y pueden añadirse otros factores o estudiar el mismo comportamiento con otras fuerzas que se comportan de manera similar como la fuerza eléctrica. 
Este sistema deja de tener solución analítica en el caso general para n>2, ya que la cantidad de grados de libertad es menor a la cantidad de ecuaciones disponibles para resolver el sistema. A pesar de esto existen casos particulares de los cuales se conocen las soluciones analíticas para más dimensiones, en especial n=3. La solucíón con condiciones arbitrarias de este sistema para los casos con n>2 puede ser encontrada mediante métodos numéricos, aproximaciones o casos reducidos de geoemtrías particulares.

Nosotros estudiamos el caso particular donde se tienen 100 masas puntuales y simétricas distribuidas aleatoriamente en un cuadrado 2D, afectadas únicamente por la fuerza gravitacional. Las velocidades iniciales se obtienen mediante una simplificación asumiendo una distribución de masa uniforme y se les es agregado un factor aleatorio proporcional a la misma. Y además si 2 de estas masas llegasen a acercarse de forma que la distancia entre ellas fuese menor a un threshold dado, estas participarían en una colisión inelástica perfecta.

Para modelar este sistema se realizó una simulación en c++, utilizando el método de RK4 para resolver el sistema de ecuaciones durante 5000 mil años de evolución. Se realizaron modelos con distintos intervalos de aumento en el tiempo y se encontró la evolución de las posiciones, velocidades, aceleraciones, energía, momento y momento angular total del sistema. 




## Métodos 
### Generalidades:
Para resolver el problema se realizó una simulación computacional en c++ bajo el estandar stdc++20. Esta simulación realiza la aproximación de la solución del sistema de movimiento empleando el método de RK4, este funciona tomando las condiciones iniciales dadas por nuestra solución particular y obtiene las variables físicas en el siguiente intervalo de tiempo. Este algoritmo se repite tomando como los valores obtenidos las nuevas condiciones iniciales y se realiza hasta que se llegue al tiempo de duración deseado.
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

## Cálculo de las variables físicas:


Para encontrar la energía total se utiliza la ecuación de la energía para un sistema de partículas:

$$ E = \frac{MR^2}{2}    +\sum \frac{m_iv_i^2}{2}  $$

Para encontrar el momento total se utilizó la ecuación del momento para un sistema de partículas:

$$ P = \sum m_iv_i $$

Para encontrar el momento agnular total se utilizó la ecuación del momento angular para un sistema de partículas:

$$ L = r \cross P + \sum r_i \cross m_i v_i $$





## Resultados
## Discusión de Resultados 
## Conclusiones 
## Referencias


## Ejercicio 5.1: 

Nos piden resolver por medio del método de euler la siguiente ecuación diferencial:

$$ y'(x) = x^2+1 $$

Con Condiciones iniciales y de dominio
$$y(0) = 0, 0x<1 $$

Con distintos valores de h= 0.05, 0.10, 0.15 y 0.20 comparados con la solución analítica.

Código: 
    //===================================
    //
    // Metodo de Euler
    // compilacion: g++ -o programa ejercicio5-1.cpp
    //
    //===================================
    
    #include <iostream>
    
    using namespace std;
    
    //utilizamos la implementación más sencilla de este método vista en clase
    double euler( double y, double t, double h );
    double derivada( double y, double t );
    
    int main()
    {
      // Datos iniciales
      const double y0 = 0;
      const double t0 = 0;
      const double h = 0.1;
    
      //Para calcular N despejamos la ecuación
      // 1 = (N-2)h
      //el -2 aparece porque empezamos a contar en 2 y nos dicen que el dominio es menor a 1
      const int N = 9; // 
      
      double y = y0;
      double t = t0;
      double y_nueva = 0.0;
    
      cout << t << "\t" << y << endl;
      
      // ciclo de iteraciones
      for( int i=0; i<=N; i++ ){
    
        y_nueva = euler( y, t, h );
    
        y = y_nueva;
        t = t + h;
    
        cout << t << "\t" << y << endl;
      }
      
      return 0;
    }
    
    
    double euler( double y, double t, double h )
    {
      return y + h*derivada( y, t );
    }
    
    
    double derivada( double y, double t )
    {
      //aquí debemos tener cuidado para cambiar al problema
      //específico que deseamos solucionar
      return y*y+1;
    }
