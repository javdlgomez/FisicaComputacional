set yrange [-2e11:2e11]
set xrange [-2e11:2e11]
dt=3600*24
set ylabel "y(m)"
set xlabel "x(m)"
set title "Problema 100 Cuerpos Programa Alterado"
#cuidado al momento de graficar ya que cambie los nombres de las bases de datos para que se diferenciaran
plot  for [i=2:101] 'posicion.dat'  using i:i+100 notitle
pause -1
#do for [it=0:365] {
 #   set title sprintf( "t = %f (s)", it*dt )
  #  plot  for [i=2:101] 'posicion.dat'  every ::it::it using i:i+100 notitle
  #  pause 0.1
#}




























