# Est-Computacional
Trabajo de investigación para la asignatura de Estadística Computacional
El programa WHAM.py calcula un potencial de fuerza media (PMF, que debe ser graficado aparte) a partir de una serie de arhivos de coordenadas.
Solo está preparado para coordenadas no cíclicas, y requiere de un archivo de metadatos que indique los archivos obtenidos de simulaciones además de la coordenada y la constante para el potencial harmónico (en unidades consecuentes).,
Se adjuntan datos de ejemplo y un archivo de metadatos de ejemplo, con la salida que dan estos datos.

Para ejecutar usar el comando

./WHAM.py rmin rmax nbins tol T metadatafile freefile

donde rmin y rmax son los puntos mínimos y máximos del histograma, nbins el número de casillas, tol la tolerancia para determinar convergencia de la energía libre, T la temperatura de las simulaciones, metadatafile el archivo de metadatos y freefile el archivo de salida.
El archivo de ejemplo se produjo con el siguiente comando

./WHAM.py 0.3 1.5 30 0.01 298.15 metadata.dat freefile.dat

Por Mateo Barría
