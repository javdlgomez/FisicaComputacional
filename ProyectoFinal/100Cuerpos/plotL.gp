set yrange [1e33:8e33]
set xrange [0:1.5768e11]
set ylabel "L(kgms)"
set xlabel "t(s)"
    set title "Lvst"
    plot  "energia_momLineal_momAngular.dat" using 1:4  w lp notitle
    pause -1 
