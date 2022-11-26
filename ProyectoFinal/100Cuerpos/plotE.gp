set yrange [5e22:4e22]
set xrange [0:1.5768e11]
set ylabel "E(J)"
set xlabel "t(s)"
    set title "Evst"
    plot  "energia_momLineal_momAngular.dat" using 1:2  w lp notitle
    pause -1 
