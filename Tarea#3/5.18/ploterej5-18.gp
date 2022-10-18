#graficador de utf8
set terminal png notransparent
set output "ej5-18xt.png"




set title "3 Cuerpos aterrizaje xvst "
set xlabel "t(s)"
set ylabel "x (m)"

plot  "datosej5-18.txt" using 1:6 w lp lw .5 lc "red" title "Cohete",  "datosej5-18.txt" using 1:4 w lp lw .5 lc "blue" title "Luna","datosej5-18.txt" using 1:2 w lp lw .5 lc "black" title "Tierra"