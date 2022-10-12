#Generador de imagenes para la animacion




#Si no hacemos esto no cargan las imágenes y solo las tira en consola

set terminal pngcairo size 1024,1024 enhanced font "Verdana,12"


#Ponemos nuestros ejes, en teoría están en metros pero así
#se sobrecarga menos la animacion

set xlabel 'x'
set ylabel 'y'

#Si no hacemos esto se comprime la imagen conforme la gráfica se vuelve horizontal o vertical
set xrange [-1:1]
set yrange [-1.5:.5]


set title ""
do for [i=0:19] {

    #colocamos la expresion regular para poder usar ffmpeg despues
    set output sprintf( "animacion/animacion%02d.png", i )

    #$2 significa que lo sacamos de la segunda columna de nuestros datos

    #esta es la cabeza del pendulo
    plot 'd5-12Ep25.txt' every ::i::i u (x1=sin($2)):(y1=-cos($2)) linetype rgb "red"  lw 100 title "" , \
	 '' every ::i::i u (0.0):(0.0):(x1):(y1) with vectors nohead  title ""
    #aqui hacemos la cuerda del pendulo que va desde el origen hasta el punto del del dataset
    print i
}


# Para correrlo necesitamos ffmpeg en path o en la misma carpeta y lo corremos con
# ffmpeg -r 10 -i animacion$02d.png -pix_fmt yuv420p animacion.mp4