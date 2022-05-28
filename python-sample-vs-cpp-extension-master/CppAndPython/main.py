import math
from turtle import title
from unittest import result
import matplotlib
import numpy as np
from tenacity import retry
from MagnetismSolver import Solve, SolveMagnetics
from Scripts.MeshGenerator import RotateMeshByAngle
from Scripts.MeshReader import MeshAnalysis, ReadRaw, ReadSaved
from Scripts.MeshWriter import SaveMesh
import pandas as pd
import matplotlib.tri as tri
import sys
from Scripts.Plotters import ComapairH, CompairLocalNuselts, CompairMeshes, DynamycsPlot, MagneticsPlot, PlotMaxH, PlotMesh, PlotPsiT_Table, PlotPsiTable, SavePlot, ShowPlot
from Scripts.ResultAnalysis import DynamycsAnalysis, MagneticsAnalysis, NuseltTable
from math import atan2, exp, sqrt
import matplotlib.pyplot as plt

from Scripts.ResultFileHandling import ResultSaving
from Scripts.settings import MagneticsResultName, MeshNames, ParamsSettings, ResultName
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
    PlotMesh(nodes, triangles, segment_indices, False, False, False)
    
def CompairSavedMeshes(*meshes):
    loaded = [ReadSaved(f"SavedMeshes/{name}.dat") for name in meshes]
    CompairMeshes(*loaded)

def PlotNuselt():
    def fit_funct(x, mult_coef, exp_coef):
        return mult_coef*x**exp_coef
        
    w, h = matplotlib.rcParams["figure.figsize"] 
    fig, ax = plt.subplots(figsize=(w*2.5, h))

    ram_min = 275

    x_fit = np.arange(ram_min, 100001, 25)

    table = NuseltTable.LoadFromCSV()
    df = table.table.sort_values(by='Ram', ascending=False)

    # fits
    meshes = MeshNames.mesh_list
    meshes = [MeshNames.n0_600, MeshNames.n2_600_dr_03_rot, MeshNames.n2_600_dr_03]
    meshes = { 
        MeshNames.n0_600 : ["o", "Circle"]
        # , MeshNames.n2_600_dr_03_rot : ["_", "Horizontal"]
        # , MeshNames.n2_600_dr_03 : ["|", "Vertical"]
    }
    for (mesh, (marker, label)) in meshes.items():
        results = df.loc[df[NuseltTable.mesh_name_literal] == mesh]
        
        # all

        rams = results[[NuseltTable.ram_literal]].to_numpy().flatten()
        nus = results[[NuseltTable.nuselt_result_literal]].to_numpy().flatten()

        ax.scatter(rams, nus, marker = marker, label = label)

        # fit

        results =  results[results[NuseltTable.ram_literal] >=  ram_min]

        rams = results[[NuseltTable.ram_literal]].to_numpy().flatten()
        nus = results[[NuseltTable.nuselt_result_literal]].to_numpy().flatten()

        from scipy.optimize import curve_fit
        popt, pcov = curve_fit(fit_funct, rams, nus)
        
        # ax.plot(x_fit, fit_funct(x_fit, *popt),
        #     label=f"{label} fit: mult=%5.4f, pow=%5.4f" % tuple(popt))

    nus = nus[::-1]
    rams = rams[::-1]
    dir = (nus[1:]-nus[:-1])/(rams[1:] - rams[:-1])
    dir2 = (dir[1:]-dir[:-1])/(rams[1:-1]-rams[:-2])
    for i, v in enumerate(dir2):
        print(rams[i],v)


    ax.legend()
    ax.set_xlabel("Ram")
    ax.set_ylabel("Nu")
    # ax.set_xscale('log')
    # fig.savefig("SavedPlots/n0_nus_plot.png")
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

    for i, xi in enumerate(x):
        ax.annotate(f"({ x[i] },{ nus[i] :.2f})", (x[i], nus[i]), 
        size=15)

    ax.set_xlabel("количество точек в области 2", size=15)
    ax.set_ylabel("Nu", size=15)
    # ax.set_title("Mesh size validation")
    fig.savefig("SavedPlots/nus_validation.png")
    plt.show()

def walk():
    import os
    from os import path

    table = NuseltTable.LoadFromCSV()

    for current_path, subdirs, files in os.walk("SavedResults\\Computational"):
        if "nodes.csv" in files:
            print(current_path)
            result = DynamycsAnalysis("", os.path.join(current_path))
            mesh_name_full = result.GetMeshName()
            ram = result.GetParam("Ram")
            name = ResultName.MakeName(mesh_name_full, ram)
            table.GetNuselt(name)
 
def PLotNusPlotly():
    import plotly.express as px
    from scipy.optimize import curve_fit
    import plotly.graph_objects as go
    def fit_funct(x, mult_coef, exp_coef):
        return mult_coef*(x**exp_coef)

    table = NuseltTable.LoadFromCSV()
    df = table.table.sort_values(by='Ram', ascending=False)

    fig = px.scatter(df, x="Ram", y="nu", color='mesh_name_short', symbol="mesh_name_short")

    # fits
    meshes = set(df[NuseltTable.mesh_name_literal])
    x_fit = np.arange(100, 100001, 25)
    for mesh in meshes:
        if not mesh in MeshNames.mesh_list:
            continue

        results = df.loc[df[NuseltTable.mesh_name_literal] == mesh]

        if (len(results) < 2):
            continue
        
        rams = results[[NuseltTable.ram_literal]].to_numpy().flatten()
        nus = results[[NuseltTable.nuselt_result_literal]].to_numpy().flatten()

        popt, pcov = curve_fit(fit_funct, rams, nus)

        # print(fit_funct(x_fit, *popt))
        name = f"{mesh}, fit mult=%5.4f, pow=%5.4f" % tuple(popt)
        hovertemplate = f"{mesh},<br> fit mult=%5.4f,<br> pow=%5.4f <extra></extra>" % tuple(popt)
        fig.add_trace(go.Scatter(x=x_fit, y=fit_funct(x_fit, *popt), mode='lines', name=name, hovertemplate=hovertemplate))


    fig.write_html("SavedPlots/Nus.html")
    fig.show()

def RotateMeshAndSave(mesh_name_to_rotate, angle, result_mesh_name):
    from os import path
    if path.exists(f"SavedMeshes/{result_mesh_name}.dat"):
        return 

    mesh = ReadSaved(f"SavedMeshes/{mesh_name_to_rotate}.dat")
    new_mesh = RotateMeshByAngle(mesh, angle)
    SaveMesh("SavedMeshes", f"{result_mesh_name}", *new_mesh)

    import os
    try:
        os.makedirs(f"SavedResults/{result_mesh_name}")
    except FileExistsError:
        print("results directory exists")

    try:
        os.makedirs(f"SavedMagnetics/{result_mesh_name}/magnetics_H_5_chi0_2_mu_1000")
    except FileExistsError:
        print("magnetics directory exists")

def PlotNusVsAngle():
    import plotly.express as px
    import plotly.graph_objects as go

    meshes = [MeshNames.n2_rotated_format.format(angle=angle) for angle in range(0,91,9)]

    table = NuseltTable.LoadFromCSV()
    df = table.table.sort_values(by='Ram', ascending=False)
    df = df.loc[df[NuseltTable.mesh_name_literal].isin(meshes)]

    default_mesh = MeshNames.n2_rotated_format.format(angle = 0)

    angles = []
    fraction = []
    for index, row in df.iterrows():
        mesh_name = row[NuseltTable.mesh_name_literal]
        angle = float(mesh_name.split("_")[-1])
        angles.append(angle)

        ram = row["Ram"]

        result_name = ResultName.MakeName(default_mesh, ram)
        


    df['Angle'] = angles
    df = df.sort_values(by='Angle')

    

    fig = px.scatter(df, x="Angle", y="nu", group="Ram")
    fig.write_html("SavedPlots/n2_angles_Nus.html")
    fig.show()

def PlotNusVsAngleV2():
    meshes = [MeshNames.n2_rotated_format.format(angle=angle) for angle in range(0,91,9)]

    table = NuseltTable.LoadFromCSV()
    df = table.table.sort_values(by='Ram', ascending=False)
    df = df.loc[df[NuseltTable.mesh_name_literal].str.contains(MeshNames.n2_rotated_format.format(angle=""))]

    angles = []
    fraction = []
    for index, row in df.iterrows():
        mesh_name = row[NuseltTable.mesh_name_literal]
        angle = float(mesh_name.split("_")[-1])
        angles.append(angle)
        


    df['Angle'] = angles
    df = df.sort_values(by='Angle')

    fig, ax = plt.subplots()
    
    rams = {
        1000 : "v"
        , 50000 : "s"
        , 100000 : "^"
    }

    for ram, marker in rams.items():
        results = df[df['Ram'] == ram]  

        angles = results['Angle'].to_numpy()
        nus = results['nu'].to_numpy()
        nus = nus/nus[0]

        ax.plot(angles, nus)
        ax.scatter(angles, nus, label = f"Ram = {ram}", marker = marker)

        if ram == 1000:
            for i, nu in enumerate(nus):
                if (i > len(nus)-5):
                    ax.annotate(f"{nus[i] :.2f}", (angles[i], nus[i]*1.005))

    ax.legend()
    ax.set_ylabel('Nu/Nu0')
    ax.set_xlabel("alpha")

    plt.show()


def main():
    walk()
    PlotNusVsAngleV2()
    return
    # PlotNuselt()
    # return
    # PlotNusVsAngleV2()
    # return
    # mesh_name = MeshNames.n2_600_dr_03_rot
    # ram = 100000
    # name = ResultName.MakeName(mesh_name, ram)
    # result = DynamycsAnalysis("SavedResults", name)
    # plotter = DynamycsPlot(result)
    # plotter.PlotPsiT()

    # results = []
    # mesh_names = [MeshNames.n3_600_dr_03, MeshNames.n4_600_dr_03, MeshNames.n4_600_dr_03_rot, MeshNames.n5_600_dr_03]

    # for mesh in mesh_names:
    #    name = MagneticsResultName.MakeName(mesh)
    #    result = MagneticsAnalysis("SavedMagnetics", name)
    #    results.append(result)

    # ComapairH(*results)
    # return

    # results = []
    # mesh_names = [MeshNames.n2_rotated_format.format(angle=angle) for angle in range(0, 91, 9)]

    # for ram in [1000, 50000, 100000]:
    #     res = []
    #     for mesh_name in mesh_names:
    #         name = ResultName.MakeName(mesh_name, ram)
    #         result = DynamycsAnalysis("SavedResults", name)
    #         res.append(result)

    #     results.append(res)

    # CompairLocalNuselts(results)
    # return

    starting = 0
    angle = 13

    new_mesh_name = MeshNames.n2_rotated_format.format(angle=angle)

    print(f"STARTING MESH {new_mesh_name}")

    RotateMeshAndSave(MeshNames.n2_rotated_format.format(angle=starting), angle, new_mesh_name)

    print(f"MESH IS ROTATED")

    # SolveMagnetics(mesh_name = new_mesh_name)

    print("MAGNETICS ARE SOLVED")

    for ram in [100000]:
        result_name = ResultName.MakeName(new_mesh_name, ram)

        table = NuseltTable.LoadFromCSV()
        # nus = table.GetNuselt(result_name, True)
        # if not math.isnan(nus):
        #     print(f"{ram} ALREADY HAS NU {nus}")
        #     continue

        initials = ""
        solve_fast(Ra = 0, Ram = ram, mesh_name = new_mesh_name,initials = initials, result_name = result_name, QW = 0.045)

        table = NuseltTable.LoadFromCSV()
        nus = table.GetNuselt(result_name, True)

        print(f"Ram = {ram} nu = {nus}")


    ## table = NuseltTable.LoadFromCSV()
    ## table.RedoTable()

    # mesh_name_full = MeshNames.n0_600

    # ## SaveRawMesh("n5/N100-600-600-100", mesh_name_full)
    # ## PlotSavedMesh(mesh_name_full)
    # ## SolveMagnetics(mesh_name = mesh_name_full)

    # ram_range = [375, 380, 385, 390, 395]
    # # ram_range = [405, 410, 415, 420]

    # print(F"STARTING MESH {mesh_name_full}")

    # last_result = ""
    # for ram in ram_range:    
    #    result_name = ResultName.MakeName(mesh_name_full, ram)

    #    table = NuseltTable.LoadFromCSV()
    #    nus = table.GetNuselt(result_name, False)
    #    if not math.isnan(nus):
    #        print(f"{ram} ALREADY HAS NU {nus}")
    #        continue

    #    initials =  ""
    #    print(F"SOLVING {ram}")
    #    solve_fast(Ra = 0, Ram = ram, mesh_name = mesh_name_full, result_name = result_name, initials = initials)
    #    last_result = result_name

    #    table = NuseltTable.LoadFromCSV()
    #    nus = table.GetNuselt(result_name, True)
    #    print(f"Ram = {ram} nu = {nus}")

    # return

    # mesh_name_full = MeshNames.n2_600_dr_03
    # ram_range = [14000, 16000, 18000]

    # #mesh_name_full = MeshNames.n2_600_dr_03_rot
    # #ram_range = [14000, 16000, 18000]

    # print(F"STARTING MESH {mesh_name_full}")

    # last_result = ""
    # for ram in ram_range:    
    #     result_name = ResultName.MakeName(mesh_name_full, ram)

    #     table = NuseltTable.LoadFromCSV()
    #     nus = table.GetNuselt(result_name, True)
    #     if not math.isnan(nus):
    #         print(f"{ram} ALREADY HAS NU {nus}")
    #         continue

    #     initials =  ""
    #     print(F"SOLVING {ram}")
    #     solve_fast(Ra = 0, Ram = ram, mesh_name = mesh_name_full, result_name = result_name, initials = initials)
    #     last_result = result_name

    #     table = NuseltTable.LoadFromCSV()
    #     nus = table.GetNuselt(result_name, True)
    #     print(f"Ram = {ram} nu = {nus}")


    # PlotNuselt()
    # CompairNus()



if __name__ == "__main__":
    # test()
    main()