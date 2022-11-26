set yrange [5e21:6e21]
set xrange [0:1.6e11]
set ylabel "P(kgs)"
set xlabel "t(s)"
set title "Pvst"
    plot  "energia_momLineal_momAngular.dat" using 1:3  w lp notitle
    pause -1 
