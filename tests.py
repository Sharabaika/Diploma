import matplotlib
import numpy as np
from Scripts.MeshReader import ReadRaw, ReadSaved
from Scripts.MeshWriter import SaveMesh
import pandas as pd
import matplotlib.tri as tri
import sys
from Scripts.Plotter import CoolPlots, DynamycsAnalysis, PlotNodes, ResultAnalysis

from Scripts.ResultFileHandling import ResultSaving
sys.path.insert(0, '/home/amninder/Desktop/Folder_2')

def main():
    mesh_name = "square_saved"

    # grid = open(f"meshes/{mesh_name}/grid.dat", "r")

    # left = open(f"meshes/{mesh_name}/left_bound.dat", "r")
    # right = open(f"meshes/{mesh_name}/right_bound.dat", "r")
    # bottom = open(f"meshes/{mesh_name}/bottom_bound.dat", "r")
    # upp = open(f"meshes/{mesh_name}/upper_bound.dat", "r")


    # nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadRaw(grid, [2000], 
    #     (left, [10, 1000]), (right, [11, 1000]), (bottom, [12, 1000]), (upp, [13, 1000]))

    # SaveMesh("SavedMeshes", f"{mesh_name}", nodes, triangles, segment_indices, trig_neighbors, node_neighbours)
    
    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadSaved(f"SavedMeshes/{mesh_name}.dat")

    # from MeshHandling.Plotter import PlotMesh 
    # PlotMesh(nodes, triangles, segment_indices, False, False, True)

    an = DynamycsAnalysis("SavedResults", "FinalTestV2")

    # an.PlotPsi()
    # an.PlotW()
    an.PlotErrors(xrange = (10,-1))



if __name__ == "__main__":
    main()