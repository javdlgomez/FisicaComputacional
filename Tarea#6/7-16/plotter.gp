

set yrange [-0.5:0.5]
set xlabel "x(1e-1nm)"
set xrange [-10:10]
dt = 5
do for [i=0:80] {
    set title sprintf( "t = %f (fs)", i*dt )
  #  plot 'solucion.dat' index i u 2:3 w lp, '' index i u 2:4 w lp 
    plot 'solucion.dat' index i u 2:($3**2+$4**2) w lp title "Amplitud de Onda"

    pause .09
}