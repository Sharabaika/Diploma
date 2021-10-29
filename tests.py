import numpy as np
from MeshReader import readPoints
from Plotter import PlotMesh

def main():
    grid = open("ez_mesh/mesh.dat", "r")

    central_mesh = open("ez_mesh/central_region.dat", "r")
    central_border = open("ez_mesh/central_border.dat", "r")

    inner_mesh = open("ez_mesh/inner_region.dat", "r")
    inner_border = open("ez_mesh/inner_border.dat", "r")

    outer_mesh = open("ez_mesh/outer_region.dat", "r")
    outer_border = open("ez_mesh/outer_border.dat", "r")
    
    points, triangles, segment_index, neighbors = readPoints(grid,
        (outer_mesh, 6), (outer_border, 5),
        (inner_mesh, 4), (inner_border, 3),
        (central_mesh, 2), (central_border, 1))
    # PlotMesh(points, triangles, segment_index, True, True, True)
    print(neighbors)

    

if __name__ == "__main__":
    main()