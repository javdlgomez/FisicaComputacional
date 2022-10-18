#graficador de utf8
set terminal png notransparent
set output "ej5-17xy.png"




set title "3 Cuerpos Aterrizaje XvsY "
set xlabel "x(Mm)"
set ylabel "y(Mm)"

plot  "datosej5-17.txt" using ($6/1e6):($7/1e6) w lp lw .5 lc "red" title "Cohete",  "datosej5-17.txt" using ($4/1e6):($5/1e6) w lp lw .5 lc "blue" title "Luna"