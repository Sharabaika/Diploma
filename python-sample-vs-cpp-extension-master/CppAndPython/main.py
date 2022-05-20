from turtle import title
from unittest import result
import matplotlib
import numpy as np
from MagnetismSolver import Solve, SolveMagnetics
from Scripts.MeshReader import MeshAnalysis, ReadRaw, ReadSaved
from Scripts.MeshWriter import SaveMesh
import pandas as pd
import matplotlib.tri as tri
import sys
from Scripts.Plotters import MagneticsPlot, PlotMesh, SavePlot, ShowPlot
from Scripts.ResultAnalysis import DynamycsAnalysis, MagneticsAnalysis, NuseltTable
from math import atan2, exp, sqrt
import matplotlib.pyplot as plt

from Scripts.ResultFileHandling import ResultSaving
from Scripts.settings import MeshNames, ParamsSettings, ResultName
from Solver import solve_fast


def PlotResults(result_name):
    results = DynamycsAnalysis("SavedResults", f"{result_name}")

    results.PlotPsi()
    

def SaveRawMesh(path, mesh_name):
    grid = open(f"MeshProjects/{path}/grid.dat", "r")
    print("grid is open")

    conductor = open(f"MeshProjects/{path}/conductor_region.dat", "r")
    print("conductor is open")
    conductor_border = open(f"MeshProjects/{path}/conductor_border.dat", "r")
    print("conductor_border is open")

    medium = open(f"MeshProjects/{path}/medium_region.dat", "r")
    print("medium is open")
    medium_border = open(f"MeshProjects/{path}/medium_border.dat", "r")
    print("medium_border is open")

    void = open(f"MeshProjects/{path}/void_region.dat", "r")
    print("void is open")
    void_border = open(f"MeshProjects/{path}/void_border.dat", "r")
    print("void_border is open")

    nodes, triangles, segment_indices, trig_neighbors, node_neighbours, trianlge_indices = ReadRaw(grid, -1, 
        (conductor, 0), (void, 4), (medium, 2), (void_border, 5), (medium_border, 3), (conductor_border, 1))
    print("raw is red")


    SaveMesh("SavedMeshes", f"{mesh_name}", nodes, triangles, segment_indices, trig_neighbors, node_neighbours, trianlge_indices)
    print("mesh is saved")


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


def LocalNulselt(result_name):
    results = DynamycsAnalysis("SavedResults", f"{result_name}")
    return results.CalculateLocalNulselt()

def Nulselt(result_name):
    results = DynamycsAnalysis("SavedResults", f"{result_name}")
    return results.CalculateNulselt()

def PlotSavedMesh(name):
    nodes, triangles, segment_indices, trig_neighbors, node_neighbours, triangle_indeces = ReadSaved(f"SavedMeshes/{name}.dat")
    PlotMesh(nodes, triangles, segment_indices, False, True, False)
    
def PlotNuselt():
    def fit_funct(x, mult_coef, exp_coef):
        return mult_coef*x**exp_coef
        
    w, h = matplotlib.rcParams["figure.figsize"] 
    fig, ax = plt.subplots(figsize=(w*2.5, h))

    ram_range = ParamsSettings.ram_range_short
    x_fit = np.arange(100, 101000, 1000)
    xdata = ram_range

    table = NuseltTable.LoadFromCSV()

    for mesh in MeshNames.mesh_list:
        nus = []
        for ram in ram_range:    
            nu = table.GetNuselt(ResultName.MakeName(mesh, ram))
            nus.append(nu)
        ydata = np.array(nus)

        from scipy.optimize import curve_fit
        popt, pcov = curve_fit(fit_funct, xdata, ydata)

        ax.plot(x_fit, fit_funct(x_fit, *popt),
            label=f"{MeshNames.GetShortName(mesh)}: mult=%5.4f, pow=%5.4f" % tuple(popt))

        ax.scatter(ram_range, nus)
    ax.legend()
    fig.savefig("nus_plot.png", dpi = 1000)
    plt.show()

def CompairNus():
    w, h = matplotlib.rcParams["figure.figsize"] 
    fig, ax = plt.subplots(figsize=(w*2.5, h))

    ram = 100000
    table = NuseltTable.LoadFromCSV()

    nus = []
    x = []

    for mesh in MeshNames.mesh_list_n0:
        nu = table.GetNuselt(ResultName.MakeName(mesh, ram))
        nus.append(nu)
        an = MeshAnalysis(mesh)
        x.append(an.GetNodesNumber(2))

    ax.plot(x, nus)
    ax.scatter(x,nus)
    ax.set_xlabel("number of nodes in region 2")
    ax.set_ylabel("Nu")
    ax.set_title("Mesh size validation")
    # fig.savefig("nus_validation.png")
    plt.show()

def main():
    mesh_name_full = MeshNames.n4_600_dr_03_rot

    # SaveRawMesh("n4/N100-600-600-100_rotated", mesh_name_full)
    # PlotSavedMesh(mesh_name_full)


    # SolveMagnetics(mesh_name = mesh_name_full)
    
    chunks = ParamsSettings.ram_chunks(4)

    ram_range = chunks[0]

    print(F"STARTING CHUNK {ram_range}")

    last_result = ""
    for ram in ram_range:    
        result_name = ResultName.MakeName(mesh_name_full, ram)

        table = NuseltTable.LoadFromCSV()
        nus = table.GetNuselt(result_name)
        if nus is not np.NaN:
            continue

        initials =  last_result
        print(F"SOLVING {ram}")
        solve_fast(Ra = 0, Ram = ram, mesh_name = mesh_name_full, result_name = result_name, initials = initials)
        last_result = result_name

        table = NuseltTable.LoadFromCSV()
        nus = table.GetNuselt(result_name, True)
        print(f"Ram = {ram} nu = {nus}")

    # PlotNuselt()
    # CompairNus()



if __name__ == "__main__":
    # test()
    main()