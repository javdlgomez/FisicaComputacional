set yrange [3e22:4e22]
set xrange [0:1.6e11]
set ylabel "E(J)"
set xlabel "t(s)"
    set title "Evst"
    plot  "energia_momLineal_momAngular.dat" using 1:2  w lp notitle
    pause -1 
