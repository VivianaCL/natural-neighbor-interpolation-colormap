import numpy as np
import math
import matplotlib.pyplot as plt
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
        radius = vor.points.ptp().max()

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
np.random.seed(1234)
points = np.random.rand(15, 2)
alturas_iniciales=[40.0,500.0,20.0,300.0,190.0,1990.0,250.0,2000.0,100.0,220.0,1700.0,30.0,200.0,400.0,1750.0]
alturas_ordenadas=[]
alturas_ordenadas=sorted(alturas_iniciales,key = lambda x:float(x))

yellow = Color("yellow")
colors = list(yellow.range_to(Color("blue"),len(alturas_ordenadas)))
vor = Voronoi(points)
regionsb, verticesb = voronoi_finite_polygons_2d(vor)
for region in regionsb:
    i=regionsb.index(region)
    polygonb = verticesb[region]
    figure= plt.fill(*zip(*polygonb), alpha=.8,lw=1, ec='k')             
    plt.setp(figure, facecolor=f'{colors[alturas_ordenadas.index(alturas_iniciales[i])]}')
for i in range(len(points)):
    x=points[i][0]
    y = points[i][1]
    plt.plot(x, y, 'bo')
    plt.text(x * (1 + 0.01), y * (1 + 0.01) , alturas_iniciales[i], fontsize=5.5)
print(points)
#plt.plot(points[:,0], points[:,1], 'ko')
plt.xlim(vor.min_bound[0] - 0.1, vor.max_bound[0] + 0.1)
plt.ylim(vor.min_bound[1] - 0.1, vor.max_bound[1] + 0.1)
plt.show()
