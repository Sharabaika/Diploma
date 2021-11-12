import numpy as np
from MeshReader import readPoints
from Plotter import PlotMesh

def main():
    grid = open("meshes/angles/mesh.dat", "r")

    central_mesh = open("meshes/angles/central_region.dat", "r")
    central_border = open("meshes/angles/central_border.dat", "r")

    inner_mesh = open("meshes/angles/inner_region.dat", "r")    
    inner_border = open("meshes/angles/inner_border.dat", "r")

    outer_mesh = open("meshes/angles/outer_region.dat", "r")    
    outer_border = open("meshes/angles/outer_border.dat", "r")

    nodes, triangles, segment_indices, neighbors = readPoints(grid,
        (outer_mesh, 6), (outer_border, 5),
        (inner_mesh, 4), (inner_border, 3),
        (central_mesh, 2), (central_border, 1))

    PlotMesh(nodes, triangles, segment_indices, False, True, False)

    

if __name__ == "__main__":
    main()