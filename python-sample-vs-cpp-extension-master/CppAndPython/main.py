import matplotlib
import numpy as np
from Scripts.MeshReader import ReadRaw, ReadSaved
from Scripts.MeshWriter import SaveMesh
import pandas as pd
import matplotlib.tri as tri
import sys
from Scripts.Plotter import CoolPlots, DynamycsAnalysis, PlotNodes, ResultAnalysis

from Scripts.ResultFileHandling import ResultSaving

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
    
    #nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadSaved(f"SavedMeshes/{mesh_name}.dat")

    # from MeshHandling.Plotter import PlotMesh 
    # PlotMesh(nodes, triangles, segment_indices, False, False, True)

    #from superfastcode import test_fun
    #print(test_fun((1,2)))

    print("aboba123")

    from superfastcode import SolveFluids

    #X = [0,0,1]
    #Y = [20,30, 40]
    #Triangles = [[10, 11], [100,111]]
    #Segments = [[1,2,3], [1,2,4]]
    #Trig_neighbours = [[1,2,3], [1,2,4]]
    #Node_neighbours = [[1,2,3], [1,2,4]]

    mesh_name = "square_saved"
    nodes, Triangles, Segments, Trig_neighbours, Node_neighbours = ReadSaved(f"SavedMeshes/{mesh_name}.dat")
    X = nodes[:,0]
    Y = nodes[:,1]

    Re = 300
    Vx = 10

    QPsi = 1.2
    QW = 0.5

    Max_error = 1e-5
    Max_cycles = 10000

    print(SolveFluids((X,Y, Triangles, Segments, Trig_neighbours, Node_neighbours, Re, Vx, QPsi, QW, Max_error, Max_cycles)))

    print("aboba")

    #an = DynamycsAnalysis("SavedResults", "ReworkedRe1000")

    #an.PlotPsi()
    #an.PlotW()
    #an.PlotErrors()



if __name__ == "__main__":
    main()