import numpy as np
from MeshReader import readPoints
from Plotter import PlotMesh

def main():
    grid = open("ez_mesh/full_mesh.dat", "r")
    inner_border = open("ez_mesh/inner_border.dat", "r")
    inner_grid = open("ez_mesh/inner_grid.dat", "r")
    center_border = open("ez_mesh/center_border.dat", "r")
    outer_grid = open("ez_mesh/outer_grid.dat", "r")
    outer_border = open("ez_mesh/outer_border.dat", "r")
    points, triangles, segment_index, neighbors = readPoints(grid, (inner_grid, 1), (outer_grid, 2), (inner_border, 3), (center_border, 4),  (outer_border, 5))
    PlotMesh(points, triangles, segment_index, False, True)

if __name__ == "__main__":
    main()