#graficador de utf8
set terminal png notransparent
set output "ej5-17y.png"




set title "3 Cuerpos aterrizaje yvst"
set xlabel "t(s)"
set ylabel "y (m)"

plot  "datosej5-17.txt" using 1:7 w lp lw .5 lc "red" title "Cohete",  "datosej5-17.txt" using 1:5 w lp lw .5 lc "blue" title "Luna"