#Generador de imagenes para la animacion




#Si no hacemos esto no cargan las imágenes y solo las tira en consola

set terminal pngcairo size 1024,1024 enhanced font "Verdana,12"


#Ponemos nuestros ejes, en teoría están en metros pero así
#se sobrecarga menos la animacion

set xlabel 'x (1e8m)'
set ylabel 'y (1e8m)'

#Si no hacemos esto se comprime la imagen conforme la gráfica se vuelve horizontal o vertical
set xrange [-5:5]
set yrange [-5:5]




set title ""
do for [i=0:46] {

    #colocamos la expresion regular para poder usar ffmpeg despues
    set output sprintf( "animacion/animacion%02d.png", i )

    #$2 significa que lo sacamos de la segunda columna de nuestros datos

    #esta es la cabeza del pendulo
    plot 'datosej5-18.txt' every ::i::i u ($6/1e8):($7/1e8) linetype rgb "black"  lw 5 title "cohete" , 'datosej5-18.txt' every ::i::i u ($2/1e8):($3/1e8) linetype rgb "green"  lw 40 title "tierra", 'datosej5-18.txt' every ::i::i u ($4/1e8):($5/1e8) linetype rgb "blue"  lw 20 title "luna"
    #aqui hacemos la cuerda del pendulo que va desde el origen hasta el punto del del dataset
  
}


# Para correrlo necesitamos ffmpeg en path o en la misma carpeta y lo corremos con
# ffmpeg -r 10 -i animacion$02d.png -pix_fmt yuv420p animacion.mp4