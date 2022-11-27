#Generador de imagenes para la animacion

#Si no hacemos esto no cargan las imágenes y solo las tira en consola

set terminal pngcairo size 1024,1024 enhanced font "Verdana,12"


#Ponemos nuestros ejes, en teoría están en metros pero así
#se sobrecarga menos la animacion

set xlabel 'x (m)'
set ylabel 'y (m)'

#Si no hacemos esto se comprime la imagen conforme la gráfica se vuelve horizontal o vertical
set xrange [-1.5e11:1.5e11]
set yrange [-1.5e11:1.5e11]

dt = 1

set title ""
do for [i=0:300] {
    set title sprintf( "t = %f (d)", it*dt*6050 )
    #colocamos la expresion regular para poder usar ffmpeg despues
    set output sprintf( "animacion/animacion%04d.png", i )

    #$2 significa que lo sacamos de la segunda columna de nuestros datos

    #esta es la cabeza del pendulo
    plot  for [i=2:101] 'posicion.dat'  every ::it::it using i:i+100 notitle
  
}



# Para correrlo necesitamos ffmpeg en path o en la misma carpeta y lo corremos con
# ffmpeg -r 10 -i animacion$02d.png -pix_fmt yuv420p animacion.mp4