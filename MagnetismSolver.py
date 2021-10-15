import numpy as np
from MeshReader import readPoints
import Plotter
from math import sqrt

grid = open("grid.dat", "r")
outer_border = open("border.dat", "r")
nodes, triangles, segment_indices, neighbors = readPoints(grid, (outer_border, 1))
# PlotMesh(nodes, triangles, segment_indices)

N_nodes = len(nodes)
N_trigs = len(triangles)

# PARAMS #
# ====== #
QF = 1.2
chi0 = 3
H0 = 1

# Arrays #
# ====== #
H = np.zeros(N_trigs)
Mu = np.zeros(N_trigs) + chi0
Fi = np.zeros(N_nodes)    

# Init Arrays #
# ----------- #        
for n_trig, triangle in enumerate(triangles):
    isBorder = any(segment_indices[n_node] == 1 for n_node in triangle)
    Mu[n_trig] = 1 if isBorder else 1 + chi0

# Solve #
# ===== #
N_cyclies = 10
for n_cycle in range(N_cyclies):
    for n_node in range(N_nodes):

        # Fi #
        # -- #
        a0F = 0
        anbF = 0
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

            a0F += a0F + Mu[n_neigbor]*0.5/Delta*( y12*y12 + x21*x21 )
            anbF += anbF + Mu[n_neigbor]*0.5/Delta*(Fi[n1]*(y12*y20+x21*x02)+Fi[n2]*(y12*y01+x21*x10))  
        FiOld = Fi[n_node]    

        segment_index = segment_indices[n_node]
        if segment_index == 0:
            Fi[n_node] =- anbF/a0F*QF + FiOld*(1-QF)
        elif segment_index == 1 :
            Fi[n_node] = H0 * nodes[n_node][1]

    for n_triangle, triangle in enumerate(triangles):
        # continue
        n0, n1, n2 = triangle
        x0, y0 = nodes[n0]
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        x10, y01 = x1-x0, y0-y1
        x21, y12 = x2-x1, y1-y2
        x02, y20 = x0-x2, y2-y0
        Delta = x10*y20-x02*y01

        DeltaA = Fi[n0]*y12 + Fi[n1]*y20 + Fi[n2]*y01
        DeltaB = Fi[n0]*x21 + Fi[n1]*x02 + Fi[n2]*x10

        # H #
        # - #
        Hx = DeltaA/Delta
        Hy = DeltaB/Delta 
        H[n_triangle] = sqrt( Hx**2+Hy**2 );   

        # Mu #
        # -- #
        segment_index = segment_indices[n_node]
        if segment_index == 0:
            Mu[n_node] = 1 + chi0*H[n_triangle]
        elif segment_index ==1 :
            Mu[n_node] = 1
                    
Plotter.PlotFi(nodes, Fi)
