from turtle import title
from unittest import result
import matplotlib
import numpy as np
from MagnetismSolver import SolveMagnetics
from Scripts.MeshReader import ReadRaw, ReadSaved
from Scripts.MeshWriter import SaveMesh
import pandas as pd
import matplotlib.tri as tri
import sys
from Scripts.Plotters import MagneticsPlot, PlotMesh, SavePlot, ShowPlot
from Scripts.ResultAnalysis import DynamycsAnalysis, MagneticsAnalysis
from math import atan2, exp, sqrt
import matplotlib.pyplot as plt

from Scripts.ResultFileHandling import ResultSaving
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
    PlotMesh(nodes, triangles, segment_indices, False, False, True)
    
def PlotNuselt(folder):
    mesh_name = "N120_n0_R1_dr0"
    nus0 = []
    ram0_range = np.array([1000, 5000, 20000, 30000, 40000, 50000, 60000, 70000, 80000])
    for ram in ram0_range:    
       result_name = f"{folder}_{ram}"
       nus0.append(Nulselt(result_name))

    from scipy.optimize import curve_fit

    def fit_funct(x, mult_coef, exp_coef):
        return mult_coef*x**exp_coef

    xdata = ram0_range
    ydata = np.array(nus0)

    popt, pcov = curve_fit(fit_funct, xdata, ydata)

    x_fit = np.arange(1000, 100000, 1000)
    plt.plot(x_fit, fit_funct(x_fit, *popt), 'r-',
         label='fit: mult=%5.3f, pow=%5.3f' % tuple(popt))

    plt.scatter(ram0_range, nus0, label="cpp")


    #plt.show()

def main():
    # ram_range = [70000, 80000, 90000]

    # last_result = f"Computations/n0/n0_N100-500-500-100_Ram_{60000}"
    # mesh_name = "n0_N100-500-500-100"
    # mesh_name_full = "Computational/n0_N100-500-500-100"
    # for ram in ram_range:    
    #     result_name = f"Computations/n0/n0_N100-500-500-100_Ram_{ram}"
    #     result_name_save = f"SavedResults/{result_name}"
    #     initials = last_result
    #     solve_fast(Ra = 0, Ram = ram, mesh_name = mesh_name_full, result_name = result_name_save, initials = initials)
    #     last_result = result_name
    
    # PlotSavedMesh("Computational/n0_N100-500-500-100")
    # PlotNuselt(f"Computations/n0/n0_N100-500-500-100_Ram")
    # PlotNuselt("N120_n0_R1_dr0/validation_final_cpp/ra_0_H_5_chi0_2_Pr_700_ram")
    # plt.legend()
    # plt.savefig("disaster")
    # plt.show()

    SolveMagnetics()
    

    # results = MagneticsAnalysis("SavedMagnetics", f"Computational/n0_N100-500-500-100/magnetics_H_5_chi0_2_mu_1000")
    # results = MagneticsAnalysis("SavedMagnetics", f"N120_n0_R1_dr0/magnetics_H_5_chi0_2_N120_n0_R1_dr0")

    # plotter = MagneticsPlot(results)
    # plotter.PlotMu()
    # ShowPlot()

    # plotter.PlotH()
    # ShowPlot()

if __name__ == "__main__":
    # test()
    main()