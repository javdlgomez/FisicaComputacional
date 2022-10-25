

set yrange [-2.5:2.5]
dt=0.0001
set ylabel "y(m)"
set xlabel "x(m)"
do for [it=0:200] {
    set title sprintf( "t = %f (s)", it*dt )
    plot 'solucion.dat' index it u 2:3 w l  title "Solucion 200 Iter"

    pause 0.05
}