#graficador de utf8
set terminal png notransparent
set output "ej5-16-1vc.png"




set title "3 Cuerpos no Luna poca vel inicial"
set xlabel "x"
set ylabel "y"

plot "datosej5-16pc.txt" using 2:3 w lp lw 1 lc  "black"   title "Tierra", "datosej5-16pc.txt" using 6:7 w lp lw 1 lc "red" title "Cohete"