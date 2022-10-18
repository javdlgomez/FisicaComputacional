#graficador de utf8
set terminal png notransparent
set output "ej5-18.png"




set title "3 Cuerpos Regreso Órbita "
set xlabel "x (m)"
set ylabel "y (m)"

plot  "datosej5-18.txt" using 6:7 w lp lw .5 lc "red" title "Cohete",  "datosej5-18.txt" using 4:5 w lp lw .5 lc "blue" title "Luna" ,  "datosej5-18.txt" using 2:3 w lp lw .5 lc "green" title "Tierra"