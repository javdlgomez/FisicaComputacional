
set terminal qt
set yrange [0:1]
set xrange [0:1]
set ylabel "y(m)"
set xlabel "x(m)"


sp  "solucion.dat" title "Mapa de Calor Solucion" with image


save "hm.gp"

