import numpy as np
import math
import matplotlib.pyplot as plt
import random
from scipy.spatial import Voronoi,voronoi_plot_2d
from colour import Color
def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = np.ptp(vor.points, axis=0).max()

    # Diccionario de aristas para tomando como clave cada punto
    all_ridges = {}

    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices): 
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Se decostruyen las regiones infinitas
    for p1, region in enumerate(vor.point_region): #pasa por todos los puntos asociados a una region
        vertices = vor.regions[region] #para cada region hace una lista de los vertices que la componen
        
        if all(v >= 0 for v in vertices ): 
            # finite region
            new_regions.append(vertices)
            continue

        #  non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0] 

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Encontrar los vértices que la lejanos

            t = vor.points[p2] - vor.points[p1] # tangent

            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal
            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n

            far_point = vor.vertices[v2] + direction * radius
    
            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]
        # finish
        new_regions.append(new_region.tolist())
   
    return new_regions, np.asarray(new_vertices)
""" Funcion que cambia para cada vertice contenido en la celda voronoi, sus vecinos anteriores con los nuevos, se utiliza después de agregar puntos"""
def cambia_vertices(lista_valores,lista_cambia,indicev_filtrado,verticesb2):
    regiones_finales_vertices=[]
    subregion_final=[]
    for i in lista_cambia:
        for j in i:
            for k in lista_valores:
                if j==k[0]:
                    subregion_final.append(verticesb2[k[1]])
            if j in indicev_filtrado:
                subregion_final.append(vertices2[j])
        regiones_finales_vertices.append(subregion_final)
        subregion_final=[]
    #print(regiones_finales_vertices)
    return regiones_finales_vertices
"""Función que calcula el area de un poligono con el input de las coordenadas de los vértices de este"""
def area_agujetas_gauss(coordenadas):
    area=0.0
    for i in range (len(coordenadas)):
        j=(i+1)%(len(coordenadas))
        area=area+ (coordenadas[i][0]*coordenadas[j][1])
        area=area-(coordenadas[j][0]*coordenadas[i][1])
    areafinal=.5*abs(area)
    return areafinal

# puntos iniciales
np.random.seed(1234)
points = np.random.rand(15, 2)
alturas_iniciales=[40.0,500.0,20.0,300.0,190.0,1990.0,250.0,2000.0,100.0,220.0,1700.0,30.0,200.0,400.0,1750.0]
# Voronoi tesselation
vor = Voronoi(points)

# plot
regions, vertices = voronoi_finite_polygons_2d(vor)

#puntos_extra

contador=0
print("Ingresa la cantidad de puntos a interpolar, un entero postivo")
puntos_interpolar=int(input())

while contador<puntos_interpolar:
    b=[]
    m=random.uniform(-0.24,1.058)
    n=random.uniform(-.084,.98)
    b.append(m)
    b.append(n)
    #print(b)
    """
    print("ingresa la coordenada x del punto ",contador," a interpolar, en el rango (-0.24,1.058)")
    x=float(input())
    while x>1.058 or x<-0.24:
        x=float(input())
    b.append(x)
    print("ingresa la coordenada y del punto",contador," a interpolar, en el rango(-.084,.98)")
    y=float(input())
    while y>.98 or y<(-0.084):
        y=float(input())
    b.append(y)
    """
    #b=[.63,.33]
    points=np.vstack([points,b])
    vorb = Voronoi(points)

    # vamos a dejar solo la region interesada
    # plot lineas
    regionsb, verticesb = voronoi_finite_polygons_2d(vorb)
    vertices2=vertices.tolist()
    verticesb2=verticesb.tolist()

    v1=[]
    v2=[]
    ln=[]
    lnb=[]
    indicesv=[]
    indicesvb=[]
    vecinos=[]
    vecinos_todos=[]
    regiones_filtradas=[]
    #quitamos los vertices repetidos de la lista orginal y los agregamos a v1
    for d in vertices2: 
        if d not in v1:
            v1.append(d)
    #quitamos los vertices repetidos de la lista de vertices después de agregar b y los agregamos a v2
    for d2 in verticesb2: 
        if d2 not in v2:
            v2.append(d2)
    #creamos dos listas, ln e indicesv que contienen los vertices y sus respectivos indices de la lista de vertices original que no están en v2        
    for i in range(len(v1)): 
        if v1[i] not in v2:
            ln.append(v1[i])
            indicesv.append(i)
    #creamos dos listas, lnb e indicesvb que contienen los vertices y sus respectivos indices de la lista de vertices nuevos que no están en v1        
    for j in range(len(v2)):
        if v2[j] not in v1:
            lnb.append(v2[j])
            indicesvb.append(j)
    #encontramos las regiones que tienen esos vertices del diagrama original
    for i in indicesv:
        for j in regions:
            if i in j and (j not in regiones_filtradas):
               regiones_filtradas.append(j)
    #utlizando la lista de regiones anterior, encontramos los vecinos vecinos de los vertices originales restantes
    for i in indicesv:    
        vecinos.append(i)
        for j in range(len(regiones_filtradas)):
            for k in range(len(regiones_filtradas[j])):
                if i==regiones_filtradas[j][k]:
                  if 0<k<len(regiones_filtradas[j])-1:
                      if regiones_filtradas[j][k+1] not in vecinos:
                          vecinos.append(regiones_filtradas[j][k+1])
                      if regiones_filtradas[j][k-1] not in vecinos:
                          vecinos.append(regiones_filtradas[j][k-1])
                  elif k==0 and k<len(regiones_filtradas[j])-1:
                      if regiones_filtradas[j][k+1] not in vecinos:
                          vecinos.append(regiones_filtradas[j][k+1])
                      if regiones_filtradas[j][len(regiones_filtradas[j])-1] not in vecinos:
                          vecinos.append(regiones_filtradas[j][len(regiones_filtradas[j])-1])
                  elif k>0 and k==len(regiones_filtradas[j])-1:
                      if regiones_filtradas[j][0] not in vecinos:
                          vecinos.append(regiones_filtradas[j][0])
                      if regiones_filtradas[j][k-1] not in vecinos:
                          vecinos.append(regiones_filtradas[j][k-1])
        vecinos_todos.append(vecinos)
        vecinos=[]
    #de los vertices nuevos, vamos a ver cuál vertice original está unido a cada vértice nuevo
    distancia=[]
    verticenuevo_verticeo=[]
    for i in indicesvb:
        #verticenuevo_verticeo.append(i)
        for j in indicesv:
            distancia.append(math.sqrt((verticesb2[i][1]-vertices2[j][1])**2+(verticesb2[i][0]-vertices2[j][0])**2))
        a=min(distancia)
        minindice=distancia.index(a)
        verticenuevo_verticeo.append(indicesv[minindice])
        distancia=[]
    #vamos a reorganizar la lista de tal manera que los vertcies orginales esten al inicio
    #primero quitamos los vertices orginales que no estén contenidos en la figura (vertices que se construyeron para obtener poligonos finitos)
    verticeo_verticenuevo=[]
    vecinoindice=[]
    maximos=[]
    indicesv_filtrado=[]
    for j in vor.ridge_vertices:
        maximos.append(max(j))
    for v in indicesv:
        if v<max(maximos):
            indicesv_filtrado.append(v)
    #print("indices filtrados")
    for i in indicesv_filtrado:
        vecinoindice.append(i) 
        for j in range(len(verticenuevo_verticeo)):
            if(i==verticenuevo_verticeo[j]): #buscamos cada elemento de indicesv en la lista que contiene el indice al que está unido cada vértice nuevo
               vecinoindice.append( indicesvb[j])
        verticeo_verticenuevo.append(vecinoindice)
        vecinoindice=[]
    """
    A continuación vamos a crear para los vertice orginales, dos listas
    1) lista que contiene los angulos que forma cada vertice con sus vecinos nuevos (vertices nuevos)
    2)lista que contiene los angulos que forma cada vertice con sus vecinos orginales (vertices originales)
    """
    #angulos que forma los vertices orginales con los nuevos
    angulos_vertices_ocn=[]
    angulos=[]
    for i in verticeo_verticenuevo:
        angulos.append(i[0])
        for j in range(1,len(i)):
            angulos.append(math.degrees(math.atan((verticesb2[i[j]][1]-vertices2[i[0]][1])/((verticesb2[i[j]][0]-vertices2[i[0]][0])+.000000000000001))))
        angulos_vertices_ocn.append(angulos)
        angulos=[]
    #print("angulos vertices originales con vertices nuevos")
    #print(angulos_vertices_ocn)

    #angulos que forman los vertices orginales con sus vecinos originales

    angulos_vertices_oco=[]
    for j in vecinos_todos:
        angulos.append(j[0])
        for k in range(1,len(j)):
            angulos.append(math.degrees(math.atan((vertices2[j[k]][1]-vertices2[j[0]][1])/(vertices2[j[k]][0]-vertices2[j[0]][0]))))
        angulos_vertices_oco.append(angulos)
        angulos=[]
    #print("angulos vertices originales con vertices oroginales")
    #print(angulos_vertices_oco)

    #vamos a ver si hay dos angulos que se repitan para cada vertice original en ambas listas
    #tenemos  dos listas angulos_vertices_ocn y angulos_vertices_oco si de los orginales con nuevos estan en originales con orginales, remplazamos los vertices orginales por los nuevos
    #empezamos encontrando que vértices orginales voy a cambiar y por cuál vertice nuevo
    cambia=[]
    pares=[]
    for i in angulos_vertices_ocn:
        for j in range(1,len(i)):
            #print(round(i[j],6),"--")
            for k in range(1,len(angulos_vertices_oco[angulos_vertices_ocn.index(i)])):
                #print(round(angulos_vertices_oco[angulos_vertices_ocn.index(i)][k],6))
                if round(i[j],6)==round(angulos_vertices_oco[angulos_vertices_ocn.index(i)][k],6):
                    pares.append(vecinos_todos[angulos_vertices_ocn.index(i)][k])
                    pares.append(verticeo_verticenuevo[angulos_vertices_ocn.index(i)][j])
                    cambia.append(pares)
                pares=[]

    #print("cambia")
    #print(cambia)
    regiones_finales=[]
    #print("coordenadas de sub regiones finales")
    regiones_filtradas_grandes=regiones_filtradas
    #Obtengo los vértices que conforman las regiones chiquitas
    regiones_finales = cambia_vertices(cambia, regiones_filtradas, indicesv_filtrado,verticesb2)
    #cambio las regiones_filtradas las regiones grandes a listas de vertices
    regiones_finales_grandes=[]
    subregiones_finales_grandes=[]
    for i in regiones_filtradas_grandes:
        for j in i:
            subregiones_finales_grandes.append(vertices2[j])
        regiones_finales_grandes.append(subregiones_finales_grandes)
        subregiones_finales_grandes=[]
            
    # calculamos el area de esas regiones
    #print(regiones_finales)
    #print(regiones_finales_grandes)

    area_regiones_finales_grandes=[]
    area_regiones_finales=[]
    for i in regiones_finales_grandes:
        area=area_agujetas_gauss(i)
        area_regiones_finales_grandes.append(area)
    for j in regiones_finales:
        areaf=area_agujetas_gauss(j)
        area_regiones_finales.append(areaf)
    #print(area_regiones_finales_grandes)
    #print(area_regiones_finales)

    #calculo porcentaje
    porcentaje=[]
    for i in range(len(area_regiones_finales)):
        porcentaje.append(area_regiones_finales[i]/area_regiones_finales_grandes[i])
    #print(porcentaje)
    #calcula que punto pertenece a cada area
    puntos_finales=[]
    for i in regiones_filtradas:
        if i in regions:
            puntos_finales.append(regions.index(i))
    suma_de_areas=[]
    #calcula la altura del punto por el porcentaje del area chiquita de toda la region
    for i in range(len(porcentaje)):
        suma_de_areas.append(porcentaje[i]*alturas_iniciales[i])
    #se suman todas las areas    
    valor=0
    for i in suma_de_areas:
        valor=valor+i
    #print("la altura aproximada en ese punto es, ")
    #print(valor)
    valor=round(valor,3)
    alturas_iniciales.append(valor)
    contador=contador+1
# colorize
alturas_ordenadas=[]
alturas_ordenadas=sorted(alturas_iniciales,key = lambda x:float(x))
##print("alturas ordenadas")
##print(len(alturas_ordenadas))

#print("el mapa incial es")
#fig1=voronoi_plot_2d(vor)
##print("el mapa interpolado es")

#se crea una lista de colores
yellow = Color("yellow")
colors = list(yellow.range_to(Color("blue"),len(alturas_ordenadas)))

for region in regionsb:
    i=regionsb.index(region)
    polygonb = verticesb[region]
    figure= plt.fill(*zip(*polygonb), alpha=.8,lw=1, ec='k')             
    plt.setp(figure, facecolor=f'{colors[alturas_ordenadas.index(alturas_iniciales[i])]}')
for i in range(len(points)):
    x=points[i][0]
    y = points[i][1]
    plt.plot(x, y, 'bo')
    #plt.text(x * (1 + 0.01), y * (1 + 0.01) , alturas_iniciales[i], fontsize=5.5)

#plt.plot(points[:,0], points[:,1], 'ko')
plt.xlim(vorb.min_bound[0] - 0.1, vorb.max_bound[0] + 0.1)
plt.ylim(vorb.min_bound[1] - 0.1, vorb.max_bound[1] + 0.1)


plt.show()
