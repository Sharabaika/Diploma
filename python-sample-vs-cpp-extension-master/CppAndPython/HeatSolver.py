import numpy as np
from MeshReader import readPoints
from Plotter import PlotMesh
from Plotter import PlotT

grid = open("grid.dat", "r")
outer_border = open("border.dat", "r")
nodes, triangles, segment_indices, neighbors = readPoints(grid, (outer_border, 1))
# PlotMesh(nodes, triangles, segment_indices)

N_nodes = len(nodes)

T = np.zeros(N_nodes)

# PARAMS #
# ====== #
QT = 1.2

N_cyclies = 100
for n_cycle in range(N_cyclies):
    for n_node in range(N_nodes):
        a0T = 0
        anbT = 0
        for n_neigbor in neighbors[n_node]:
            neigbor = triangles[n_neigbor]
            n0, n1, n2 = neigbor
            if n_node == n1:
                n0, n1, n2 =  n_node, n2, n0
            elif n_node == n2:
                n0, n1, n2 = n_node, n0, n1

            x0, y0 = nodes[n0]
            x1, y1 = nodes[n1]
            x2, y2 = nodes[n2]
            x10, y01 = x1-x0, y0-y1
            x21, y12 = x2-x1, y1-y2
            x02, y20 = x0-x2, y2-y0
            Delta = x10*y20-x02*y01

            a0T += +0.5/Delta*( y12*y12 + x21*x21 )
            anbT += +0.5/Delta*( T[n1]*(y12*y20+x21*x02) + T[n2]*(y12*y01+x21*x10) )            
        TOld=T[n_node]
        segment_index = segment_indices[n_node]
        if segment_index == 0:
            T[n_node] = -anbT/a0T*QT + TOld*(1-QT)
        elif segment_index == 1:
            T[n_node] = 1

        if n_node == 53:
            T[n_node] = 0
                    
PlotT(nodes, T)