

set yrange [-0.5:0.5]
set xrange [-35:35]
dt = 5
set xlabel "x(1e-1nm)"
do for [i=0:30] {
    set title sprintf( "t = %f (fs)", i*dt )
  #  plot 'solucion.dat' index i u 2:3 w lp, '' index i u 2:4 w lp 
    plot 'solucion.dat' index i u 2:($3**2+$4**2) w lp title "Amplitud de Onda",  step(x) = x>0 ? .2 : 0

    pause -1}