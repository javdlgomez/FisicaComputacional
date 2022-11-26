set yrange [-3e11:3e11]
set xrange [-3e11:3e11]
dt=3600
set ylabel "y(m)"
set xlabel "x(m)"
set title "Problema 100 Cuerpos"
#cuidado al momento de graficar ya que cambie los nombres de las bases de datos para que se diferenciaran
#plot  for [i=2:101] 'posicion.dat'  using i:i+100 notitle
#pause -1
do for [it=0:15] {
   set title sprintf( "t = %f (s)", it*dt )
    plot  for [i=2:101] 'posicion.dat'  every ::it::it using i:i+100 notitle
    pause 0.01
}




























