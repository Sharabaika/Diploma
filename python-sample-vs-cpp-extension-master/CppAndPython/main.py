import matplotlib
import numpy as np
from Scripts.MeshReader import ReadRaw, ReadSaved
from Scripts.MeshWriter import SaveMesh
import pandas as pd
import matplotlib.tri as tri
import sys
from Scripts.ResultAnalysis import CoolPlots, DynamycsAnalysis, PlotMesh, PlotNodes, ResultAnalysis
from math import atan2, exp, sqrt
import matplotlib.pyplot as plt

from Scripts.ResultFileHandling import ResultSaving


def PlotResults(result_name):
    results = DynamycsAnalysis("SavedResults", f"{result_name}")

    results.PlotPsi()
    

def SaveRawMesh(name):
    mesh_name = name

    grid = open(f"MeshProjects/{mesh_name}/grid.dat", "r")

    inner = open(f"MeshProjects/{mesh_name}/inner.dat", "r")
    outer = open(f"MeshProjects/{mesh_name}/outer.dat", "r")


    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadRaw(grid, [2000], 
        (inner, [10, 1000]), (outer, [11, 1000]))

    SaveMesh("SavedMeshes", f"{mesh_name}", nodes, triangles, segment_indices, trig_neighbors, node_neighbours)


def CPPStuff():
    from superfastcode import test_fun
    print(test_fun((1,2)))

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

    from superfastcode import SolveFluids
    Psi, W, DPsi, DW = SolveFluids((X,Y, Triangles, Segments, Trig_neighbours, Node_neighbours, Re, Vx, QPsi, QW, Max_error, Max_cycles))

    print(len(DPsi))

    Saver.logger.LogErrorsList(Psi = DPsi, W = DW)
    Saver.SaveResults("SavedResults", "CPPTestFinal", W = W, Psi = Psi)


def test():

    mesh_name = "small"

    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadSaved(f"SavedMeshes/{mesh_name}.dat")
    N_nodes = len(nodes)

    import matplotlib.pyplot as plt
    ax = plt.axes()
    ax.set_aspect('equal')

    for n_node in range(N_nodes):
        segment_index = segment_indices[n_node]

        if 1000 in segment_index:            
            # normal #
            # ------ #
            normalX, normalY = 0, 0

            border_neighbours = []
            for neighbour in node_neighbours[n_node]:
                if 1000 in segment_indices[neighbour]:
                    border_neighbours.append(neighbour)

            x, y = nodes[n_node]

            l, r = border_neighbours
            xl, yl = nodes[l]
            xr, yr = nodes[r]
            
            fl = atan2(yl, xl)
            fr = atan2(yr, xr)

            if fl >= fr:
                fl, fr = fr, fl
                xl, yl, xr, yr = xr, yr, xl, yl

            if fl*fr < 0 and abs(fr > 1):
                fl, fr = fr, fl
                xl, yl, xr, yr = xr, yr, xl, yl  

            dx, dy = xr-xl, yr-yl
            dr = sqrt(dx*dx + dy*dy)
            dx, dy = dx/dr, dy/dr

            normalX, normalY = dy, -dx
            ax.arrow(x, y, normalX / 100, normalY / 100, head_width=0.005, head_length=0.01, fc='k', ec='k')


    plt.show()


def Nulselt(result_name):
    results = DynamycsAnalysis("SavedResults", f"{result_name}")

    return results.CalculateLocalNulselt()

def PlotSavedMesh(name):
    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadSaved(f"SavedMeshes/{name}.dat")
    PlotMesh(nodes, triangles, segment_indices, False, True, True)
    
def main():
    # fi, nu = Nulselt("circle_convectionV2")
    # plt.plot(fi, nu)
    # plt.show()

    meshes = ["N120_n0_R1_dr0", "N120_n3_R1_dr0.3", "N120_n4_R1_dr0.3"]
    for mesh in meshes:
        # SaveRawMesh(mesh)
        PlotSavedMesh(mesh)

if __name__ == "__main__":
    # test()
    main()