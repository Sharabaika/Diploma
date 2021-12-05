import numpy as np
from MeshFileHandling.MeshReader import ReadRaw, ReadSaved
from MeshFileHandling.MeshWriter import SaveMesh

import sys
sys.path.insert(0, '/home/amninder/Desktop/Folder_2')

def main():
    grid = open("meshes/circle/mesh.dat", "r")
    inner = open("meshes/circle/inner.dat", "r")
    outer = open("meshes/circle/outer.dat", "r")


    # nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadRaw(grid, [2000], 
    #     (inner, [10, 1000]), (outer, [11, 1000]))

    # SaveMesh("SavedMeshes", "circle_saved", nodes, triangles, segment_indices, trig_neighbors, node_neighbours)
    
    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadSaved("SavedMeshes/circle_saved.dat")

    from MeshHandling.Plotter import PlotMesh 
    PlotMesh(nodes, triangles, segment_indices, False, False, True)


if __name__ == "__main__":
    main()