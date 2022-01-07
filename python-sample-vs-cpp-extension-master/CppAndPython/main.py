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

    mesh_name = "square_saved"
    nodes, Triangles, Segments, Trig_neighbours, Node_neighbours = ReadSaved(f"SavedMeshes/{mesh_name}.dat")
    X = nodes[:,0].tolist()
    Y = nodes[:,1].tolist()

    Re = 300
    Vx = -1

    QPsi = 1.2
    QW = 0.5

    Max_error = 1e-5
    Max_cycles = 10000

    triangulation = tri.Triangulation(X,Y, Triangles)

    Saver = ResultSaving("W", "Psi")
    Saver.AddParams(mesh_name = mesh_name, Re = Re, QPsi = QPsi, QW = QW)

    #import Scripts.Plotter as plotter
    #plotter.PlotMesh(nodes, Triangles, Segments, False, False, False)

    #from superfastcode import SolveFluids
    #Psi, W, DPsi, DW = SolveFluids((X,Y, Triangles, Segments, Trig_neighbours, Node_neighbours, Re, Vx, QPsi, QW, Max_error, Max_cycles))

    #print(len(DPsi))

    #Saver.logger.LogErrorsList(Psi = DPsi, W = DW)
    #Saver.SaveResults("SavedResults", "CPPTestFinal", W = W, Psi = Psi)


    an = DynamycsAnalysis("SavedResults", "CPPTestFinal")

    an.PlotPsi()
    an.PlotW()
    an.PlotErrors()



if __name__ == "__main__":
    main()