

set yrange [-2.5:2.5]
set xrange [0:.5]
dt=0.0001
set ylabel "y(m)"
set xlabel "x(m)"

do for [it=0:20] {
    set title sprintf( "t = %f (s)", it*dt )
    plot 'solucion.dat' index it u 2:3 w l  title "Solucion 20 Iter"

    pause .5
    #como son menos alargué un poco la animación
}