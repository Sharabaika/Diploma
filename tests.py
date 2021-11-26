import numpy as np
from MeshFileHandling.MeshReader import ReadRaw, ReadSaved
from MeshFileHandling.MeshWriter import SaveMesh

import sys
sys.path.insert(0, '/home/amninder/Desktop/Folder_2')

def main():
    # grid = open("ez_mesh/mesh.dat", "r")

    # central_mesh = open("ez_mesh/central_region.dat", "r")
    # central_border = open("ez_mesh/central_border.dat", "r")

    # inner_mesh = open("ez_mesh/inner_region.dat", "r")
    # inner_border = open("ez_mesh/inner_border.dat", "r")

    # outer_mesh = open("ez_mesh/outer_region.dat", "r")
    # outer_border = open("ez_mesh/outer_border.dat", "r")

    # nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadRaw(grid,
    #     (outer_mesh, 6), (outer_border, 5),
    #     (inner_mesh, 4), (inner_border, 3),
    #     (central_mesh, 2), (central_border, 1))

    # SaveMesh("SavedMeshes", "ez_mesh_saved", nodes, triangles, segment_indices, trig_neighbors, node_neighbours)
    
    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadSaved("SavedMeshes/ez_mesh_saved.dat")

    from MeshHandling.Plotter import PlotMesh 
    PlotMesh(nodes, triangles, segment_indices)


if __name__ == "__main__":
    main()