import matplotlib.pyplot as plt
import numpy as np

def PlotMesh(points, triangles, segment_idex, index_nodes = False, scatted_nodes = False, index_regions = False):
    x, y = points[:, 0], points[:, 1]
    plt.triplot(x, y, triangles, color='green')

    if scatted_nodes:
        plt.scatter(x, y, s=100, c=segment_idex)    

    if index_nodes:
        for point_index in range(len(x)):
            plt.text(x=x[point_index], y=y[point_index], s = point_index, color='red', fontsize=14)

    if index_regions:
        for point_index in range(len(x)):
            plt.text(x=x[point_index], y=y[point_index], s = segment_idex[point_index], color='red', fontsize=14)
    plt.show()

def PlotNodes(points, Fi):
    x, y = points[:, 0], points[:, 1]

    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    cf = ax.tricontourf(x,y,Fi)
    fig.colorbar(cf, ax=ax)
    plt.show()

def PlotScatter(points, z):
    x, y = points[:, 0], points[:, 1]

    plt.scatter(x, y, s=100, c=z) 
    plt.show()

def PlotElements(triang, z):
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    tpc = ax1.tripcolor(triang, z, shading='flat')
    fig1.colorbar(tpc)
    ax1.set_title('tripcolor of Delaunay triangulation, flat shading')
    plt.show()
