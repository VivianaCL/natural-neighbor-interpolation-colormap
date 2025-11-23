# Construcción de Mapas Topográficos: Interpolación de Vecino Natural

Este repositorio contiene la implementación en Python de algoritmos de **Geometría Computacional** (Voronoi y Delaunay) para la generación de mapas topográficos a partir de datos irregulares. 

El proyecto se centra en la **Interpolación de Vecino Natural (Natural Neighbor Interpolation)**, una técnica que estima la altura de puntos desconocidos basándose en la ponderación de áreas de las celdas de Voronoi.

## Descripción

El objetivo es procesar un conjunto disperso de coordenadas con alturas conocidas y generar una superficie continua y suave. Para la visualización, se utiliza una **escala de colores gradiente**: desde tonos amarillos (alturas bajas) hasta tonos azules (alturas altas), facilitando la identificación de patrones geoespaciales[cite: 10, 17].

Se utilizan:
* **Diagramas de Voronoi:** Para la subdivisión del plano.
* **Fórmula de la "Agujeta" de Gauss (Shoelace Formula):** Para el cálculo exacto de áreas de polígonos irregulares.
* **Interpolación Ponderada:** El valor de un nuevo punto se calcula según el porcentaje de área que su celda "roba" a las celdas vecinas originales.

##  Contenido del Repositorio

El proyecto consta de 3 scripts principales:

### 1. `Voronoi Puntos Iniciales Coloreados.py`
Este script toma el conjunto de datos base (15 puntos iniciales con alturas conocidas) y genera el Diagrama de Voronoi correspondiente.
* **Función:** Visualiza la topografía inicial aplicando la escala de colores (Amarillo $\to$ Azul) a las celdas existentes según la altura de su punto generador].

### 2. `Interpolación de puntos con coordenadas ingresadas por el usuario.py`
Herramienta interactiva que permite calcular la altura específica en una coordenada arbitraria.
* **Uso:** El usuario ingresa las coordenadas $(x, y)$ de un punto $q$.
* **Proceso:** El algoritmo inserta $q$ en el diagrama, calcula su nueva celda de Voronoi, determina las subregiones de traslape con los vecinos y devuelve la altura interpolada precisa.
* **Visualización:** Muestra el mapa actualizado con el nuevo punto integrado.

### 3. `Interpolación de puntos con coordenadas aleatorias.py`
Simulación masiva para densificar el mapa topográfico.
* **Funcionamiento:** Genera automáticamente una cantidad definida de puntos aleatorios (ej. 500, 1000, 2000) dentro del rango del mapa utilizando una distribución uniforme.
##  Tecnologías y Librerías

El código fue desarrollado en Python utilizando las siguientes librerías:

* **`scipy.spatial`**: Para el cálculo de `Voronoi` y estructuras geométricas.
* **`numpy`**: Para manejo de arreglos y cálculos vectoriales.
* **`matplotlib.pyplot`**: Para la graficación de los diagramas y mapas.
* **`colour`**: Para la generación de la escala de colores gradiente.
* **`math`** y **`random`**: Para cálculos auxiliares y generación de puntos de prueba.

