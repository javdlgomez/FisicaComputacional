set yrange [-1.5e11:1.5e11]
set xrange [-1.5e11:1.5e11]
dt=3600*24
set ylabel "y(m)"
set xlabel "x(m)"
set title "Problema 3 Cuerpos Programa 42 Anos"
#cuidado al momento de graficar ya que cambie los nombres de las bases de datos para que se diferenciaran
plot 'posicion.dat'  using 2:5 t "tierra", 'posicion.dat'  using 3:6 t "Luna", 'posicion.dat'  using 4:7 t "Sol"
pause -1
#do for [it=0:365] {
 #   set title sprintf( "t = %f (s)", it*dt )
  #  plot  for [i=2:101] 'posicion.dat'  every ::it::it using i:i+100 
  #  pause 0.1
#}




























