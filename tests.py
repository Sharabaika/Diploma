import numpy as np
from MeshFileHandling.MeshReader import ReadRaw, ReadSaved
from MeshFileHandling.MeshWriter import SaveMesh

import sys
sys.path.insert(0, '/home/amninder/Desktop/Folder_2')

def main():
    mesh_name = "small_square"

    grid = open(f"meshes/{mesh_name}/grid.dat", "r")

    left = open(f"meshes/{mesh_name}/left_bound.dat", "r")
    right = open(f"meshes/{mesh_name}/right_bound.dat", "r")
    bottom = open(f"meshes/{mesh_name}/bottom_bound.dat", "r")
    upp = open(f"meshes/{mesh_name}/upper_bound.dat", "r")


    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadRaw(grid, [2000], 
        (left, [10, 1000]), (right, [11, 1000]), (bottom, [12, 1000]), (upp, [13, 1000]))

    SaveMesh("SavedMeshes", f"{mesh_name}", nodes, triangles, segment_indices, trig_neighbors, node_neighbours)
    
    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadSaved(f"SavedMeshes/{mesh_name}.dat")

    from MeshHandling.Plotter import PlotMesh 
    PlotMesh(nodes, triangles, segment_indices, False, False, True)


if __name__ == "__main__":
    main()