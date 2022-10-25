

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