# Natural Neighbor Interpolation for Topographic Map Generation

This repository contains a Python implementation of **Computational Geometry** algorithms (Voronoi and Delaunay) to generate smooth topographic surfaces from irregular elevation data.

The core method used is **Natural Neighbor Interpolation**, which estimates the height of any new point based on the relative areas of Voronoi cells.

## Overview

The goal of this project is to transform scattered points with known elevations into a continuous topographic map. A **yellow → blue color gradient** is applied to visualize low-to-high elevation patterns.

The project relies on:
- **Voronoi diagrams** for domain partitioning  
- **Gauss’s Shoelace Formula** for exact polygon area computation  
- **Area-weighted interpolation** based on Voronoi cell overlap  

## Repository Structure

### 1. `Voronoi_InitialPoints_Colored.py`
Generates the initial Voronoi diagram of the base dataset and visualizes it using the height-based color gradient.

### 2. `Interpolation_UserInput.py`
Interactive tool where the user provides coordinates \((x, y)\).  
The algorithm inserts the new point, recomputes local Voronoi regions, and returns the interpolated elevation.

### 3. `Interpolation_RandomPoints.py`
Runs large-scale simulations by generating random points to densify the terrain surface.


For more detailed explanations, mathematical derivations, and experiments, refer to the full project  `./docs/Construyendo Mapas Topograficos.pdf`
