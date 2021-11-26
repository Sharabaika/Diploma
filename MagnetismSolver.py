import matplotlib
from matplotlib.pyplot import title
import numpy as np
from MeshFileHandling.MeshReader import ReadRaw
import MeshHandling.Plotter as Plotter
from math import exp, sqrt
import matplotlib.tri as tri

def most_common(lst):
    return max(set(lst), key=lst.count)

grid = open("ez_mesh/mesh.dat", "r")

central_mesh = open("ez_mesh/central_region.dat", "r")
central_border = open("ez_mesh/central_border.dat", "r")

inner_mesh = open("ez_mesh/inner_region.dat", "r")
inner_border = open("ez_mesh/inner_border.dat", "r")

outer_mesh = open("ez_mesh/outer_region.dat", "r")
outer_border = open("ez_mesh/outer_border.dat", "r")

nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadRaw(grid,
    (outer_mesh, 6), (outer_border, 5),
    (inner_mesh, 4), (inner_border, 3),
    (central_mesh, 2), (central_border, 1))
# PlotMesh(nodes, triangles, segment_indices)

x, y = nodes[:, 0], nodes[:, 1]

triangulation = tri.Triangulation(x, y, triangles) 


N_nodes = len(nodes)
N_trigs = len(triangles)

# PARAMS #
# ====== #
QF = 1.6
chi0 = 3.0
H0 = 10000.0
mu0 = 1000000.0

# Arrays #
# ====== #
H = np.zeros(N_trigs)
Mu = np.zeros(N_trigs)
Fi = np.zeros(N_nodes)  


# Init Arrays #
# =========== #   

# Field # 
# ----- #
for n_trig, triangle in enumerate(triangles):
    H[n_trig] = H0

    segment_index = most_common([ segment_indices[n_node] for n_node in triangle ])
    if segment_index in [2]:
        Mu[n_trig] = mu0
    elif segment_index in [1,3,4]:
        Mu[n_trig] = 1 + chi0/(1+chi0*H[n_trig])
    elif segment_index in [5,6]:
        Mu[n_trig] = 1

H_new = np.array(H)
Mu_new = np.array(Mu)
Fi_new = np.array(Fi)


# Solve #
# ===== #
N_cyclies = 25
for n_cycle in range(N_cyclies):
    if n_cycle % 25 == 0:
        print("cycle n{0}".format(n_cycle))

    for n_node in range(N_nodes):

        a0F = 0
        anbF = 0

        Psi_BorderIntegral_a0 = 0
        Psi_BorderIntegral_nb = 0

        Psi_AreaIntegral = 0

        W_BorderIntegral = 0
        W_BorderIntegral_k0 = 0

        W_AreaIntegral = 0

        for n_trig_neigbor in trig_neighbors[n_node]:
            n0, n1, n2 = triangles[n_trig_neigbor]
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

            # Fi #
            # -- #
            a0F = a0F + Mu[n_trig_neigbor]*0.5/Delta*( y12*y12 + x21*x21 )
            anbF = anbF + Mu[n_trig_neigbor]*0.5/Delta*(Fi[n1]*(y12*y20+x21*x02)+Fi[n2]*(y12*y01+x21*x10))  


        segment_index = segment_indices[n_node]
        if segment_index in [1,2,3,4,6]:
            if abs(anbF) > 0 and abs(a0F) >0:
                Fi_new[n_node] = -anbF/a0F*QF + Fi[n_node]*(1-QF)
        elif segment_index in [5] :
            Fi_new[n_node] = H0 * nodes[n_node][1]
        else:
            print("keke")

    for n_triangle, triangle in enumerate(triangles):
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
        H_new[n_triangle] = sqrt( Hx**2+Hy**2 );   

        # Mu #
        # -- #
        segment_index = most_common([ segment_indices[n] for n in triangle ])
        if segment_index in [2]:
            Mu_new[n_triangle] = mu0
        elif segment_index in [1,3,4]:
            # Mu_new[n_triangle] = 1 + chi0*H[n_triangle]/(1+chi0*H[n_triangle])
            Mu_new[n_triangle] = 1 

        elif segment_index in [5,6]:
            Mu_new[n_triangle] = 1
        else:
            print("aboba")

    Mu = Mu_new
    H = H_new
    Fi = Fi_new

Plotter.CoolPlots.PlotLevelNodes(nodes, Fi, xrange = (-1.5,1.5), yrange = (-1.5,1.5), nlevels = 30, manual = True,  title="Fi")
# Plotter.CoolPlots.PlotLevelTriangles(nodes, triangles, H,  xrange = (-1.5,1.5), yrange = (-1.5,1.5), title="H")
# Plotter.CoolPlots.PlotLevelTriangles(nodes, triangles, Mu,  xrange = (-1.5,1.5), yrange = (-1.5,1.5), title="Mu")
