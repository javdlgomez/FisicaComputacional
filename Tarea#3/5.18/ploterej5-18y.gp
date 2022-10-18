#graficador de utf8
set terminal png notransparent
set output "ej5-18y.png"




set title "3 Cuerpos aterrizaje yvst"
set xlabel "t(s)"
set ylabel "y (m)"

plot  "datosej5-18.txt" using 1:7 w lp lw .5 lc "red" title "Cohete",  "datosej5-18.txt" using 1:5 w lp lw .5 lc "blue" title "Luna", "datosej5-18.txt" using 1:3 w lp lw .5 lc "black" title "Tierra"