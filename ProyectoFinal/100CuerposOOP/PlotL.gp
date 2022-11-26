set yrange [1.35e33:5e34]
set xrange [0:1.6e11]
set ylabel "L(kgms)"
set xlabel "t(s)"
    set title "Momento Angular"
    plot  "energia_momLineal_momAngular.dat" using 1:4  w lp notitle
    pause -1 
