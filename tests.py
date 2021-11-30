import numpy as np
from MeshFileHandling.MeshReader import ReadRaw, ReadSaved
from MeshFileHandling.MeshWriter import SaveMesh

import sys
sys.path.insert(0, '/home/amninder/Desktop/Folder_2')

def main():
    grid = open("meshes/square/grid.dat", "r")
    left = open("meshes/square/left_bound.dat", "r")
    right = open("meshes/square/right_bound.dat", "r")
    upper = open("meshes/square/upper_bound.dat", "r")
    bottom = open("meshes/square/bottom_bound.dat", "r")


    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadRaw(grid, [2000], 
        (left, [10, 1000]), (right, [11, 1000]),
         (bottom, [12, 1000]), (upper, [13, 1000]))

    SaveMesh("SavedMeshes", "square_saved", nodes, triangles, segment_indices, trig_neighbors, node_neighbours)
    
    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadSaved("SavedMeshes/square_saved.dat")

    from MeshHandling.Plotter import PlotMesh 
    PlotMesh(nodes, triangles, segment_indices, False, True, True)


if __name__ == "__main__":
    main()