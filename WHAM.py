#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
import sys

_, rmin, rmax, nbins, tol, T, metadatafile, freefile = sys.argv
#Comando debería tener el siguiente formato
#WHAM.py rmin rmax nbins tol T metadatafile freefile

rmin = float(rmin)
rmax = float(rmax)
nbins = int(nbins)
tol = float(tol)
T = float(T)

Na = 6.02214076e23 # mol⁻¹
k_SI = 1.38064852e-23 # E m² kg s⁻² K⁻¹ o J K⁻¹
k_kJ = k_SI / 1000 # 1.38064852e-23 * 10^-2 kJ K⁻¹
R_kJ = k_kJ * Na # kJ K⁻¹ mol⁻¹
beta = (R_kJ * T)**-1

data = [] #Datos totales
coordinates = [] #Ventanas
constants = [] #Constantes armónicas

#Leer Datos
with open(metadatafile) as m:
    for line in m:
        file, xi, k = line.split()
        data.append(np.loadtxt(file))
        coordinates.append(float(xi))
        constants.append(float(k))
    S = len(coordinates) #Contar cantidad de ventanas
    print('Cantidad de ventanas: %i' % S)
    
#Convertir de Listas a Arrays
data = np.asarray(data)
coordinates = np.asarray(coordinates)
constants = np.asarray(constants)

flat = data[:,:,1].flatten() #Aplanar datos para hacer el histograma
#Filtramos los datos para solo contar valores dentro de la coordenada de reacción como fue definida
filtered = flat[np.where(np.logical_and(flat>=rmin, flat<=rmax))] 
#Aprovechamos que los histogramas de matplotlib.pyplot devuelven las casillas y los conteos
counts, bins = plt.hist(filtered, bins = nbins)[0:2]
#counts[j] = sum_i^s n_{ij}
plt.close() #Cerramos el histograma ya que no lo necesitaremos más (no lo mostraremos)

N = [] #Ahora contamos los datos totales por simulación que caen dentro de la coordenada
for i in range(S):
    flat = data[i,:,1]
    filtered = flat[np.where(np.logical_and(flat>=rmin, flat<=rmax))]
    N.append(len(filtered))

#Calculamos los puntos al centro de las casillas
reaction = [(bins[i+1] + bins[i])/2 for i in range(len(bins) - 1)]

#Potencial de sesgo
def V(ξ,ξ0,k):
    return 1/2 * k * (ξ - ξ0)**2

#Factor de sesgo
def C(ξj,ξ0i,ki):
    return np.exp(-beta *V(ξj,ξ0i,ki))

#Calculamos el factor de sesgo para cada simulación, en cada casilla
c = []
for i in range(S): 
    fila = []
    for j in range(nbins):
        fila.append(C(reaction[j],coordinates[i],constants[i]))
    c.append(fila)
c = np.asarray(c)

#Energía libre en función de xi
def F(p):
    return -R_kJ * T * np.log(p)

#Energía libre de cada ventana
def A(f):
    return beta**-1 * np.log(f)

#Probabilidad des-sesgada, primera ecuación de WHAM
def pj(j,f):
    num = counts[j]
    den = 0
    for i in range(S):
        den += N[i] * f[i] * c[i,j]
    return num/den
        
#Factor de normalización, segunda ecuación de WHAM
def fi(i,p):
    s = 0
    for j in range(nbins):
        s += c[i,j] * p[j]
    return s**-1

#Estimaciones iniciales de las ecuaciones de WHAM
#Primero asumimos que todos los F son iguales a 1
f = np.asarray([1] * S)
#Y a partir de ellos calculamos los p
p = np.asarray([pj(j,f) for j in range(nbins)])

for l in range(10000): #Máximo de mil iteraciones
    f = np.asarray([fi(i,p) for i in range(S)]) 
    newp = np.asarray([pj(j,f) for j in range(nbins)])
    if False not in (abs(F(p) - F(newp)) < tol): 
        #Romper el loop si TODAS las energías libres convergen
        print('Convergencia en %i iteraciones\n' % l)
        p = newp
        break
    else:
        p = newp
else: #Si el loop se completa, ejecutar esta linea
    print('No converge en %i iteraciones, usando última estimación\n' % l)
    
#Calcular energías libres finales y desplazarlas para que la de xi sea igual a 0
# en su mínimo mientras que la de las ventanas sea igual a 0 en su primer punto
Free = abs(F(p) - min(F(p)))
Helmholtz = A(f) - A(f[0])

#Guardar todo a un archivo externo que fácilmente permita graficar el PMF
with open(freefile,'w') as file:
    head = '#Coor\tFree\tProb'
    print(head)
    file.write(head + '\n')
    for j in range(nbins):
        line = '%.6f\t%.6f\t%.6f' % (reaction[j] , Free[j] , p[j])
        print(line)
        file.write(line + '\n')
    #Guardar además como comentarios las energías libres de las ventanas, porsiacaso
    head2 = '#Window\tFree'
    print(head2)
    file.write(head2 + '\n')
    for i in range(S):
        line = '#%i\t%.6f' % (i, Helmholtz[i])
        print(line)
        file.write(line + '\n')
