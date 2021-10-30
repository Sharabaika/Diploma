import matplotlib
import numpy as np
from MeshReader import readPoints
import Plotter
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

nodes, triangles, segment_indices, neighbors = readPoints(grid,
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
chi0 = 3
H0 = 10000
mu0 = 1000000

Vc = 1           # Viscosity
Pr = 1

# Arrays #
# ====== #
H = np.zeros(N_trigs)
Mu = np.zeros(N_trigs)
Fi = np.zeros(N_nodes)  

Psi = np.zeros(N_nodes)
W = np.zeros(N_nodes)

# Init Arrays #
# ----------- #   

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


# Dynamics #
# -------- #
Psi_new = np.array(Psi)
W_new = np.array(W)


# Solve #
# ===== #
N_cyclies = 300
for n_cycle in range(N_cyclies):
    for n_node in range(N_nodes):

        a0F = 0
        anbF = 0

        Psi_BorderIntegral_a0 = 0
        Psi_BorderIntegral_nb = 0

        Psi_AreaIntegral = 0

        for n_trig_neigbor in neighbors[n_node]:
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

            # Psi #
            # --- #
            Delta_PsiA = Psi[n0]*y12 + Psi[n1]*y20 + Psi[n1]*y01
            Delta_PsiB = Psi[n0]*x21 + Psi[n1]*x02 + Psi[n1]*x10

            A_Psi = Delta_PsiA / Delta
            B_PSi = Delta_PsiB / Delta

            a0      = (x21*x21+y12*y12)/Delta
            a1_psi1 = Psi[n1] * (y20*y12+x02*x21)
            a2_psi2 = Psi[n2] * (y01*y12+x10*x21)

            Psi_BorderIntegral_a0 += a0
            Psi_BorderIntegral_nb += (a1_psi1+a2_psi2)/2

            Psi_AreaIntegral += (22.*W[n0]+7.*W[n1]+7.*W[n2])*Delta/216.0

            # W #
            # - # 
            U_x = B_PSi
            U_y = -A_Psi
            U_ell = sqrt(U_x*U_x + U_y*U_y)

            sina = U_y /U_ell
            cosa = U_x / U_ell

            xc = (x0+x1+x2)/3
            yc = (y0+y1+y2)/3

            def ToLocal(x, y):
                X =  (x-xc)*cosa + (y-yc)*sina
                Y = -(x-xc)*sina + (y-yc)*cosa
                return X,Y
            
            X0, Y0 = ToLocal(x0,y0)
            X1, Y1 = ToLocal(x1,y1)
            X2, Y2 = ToLocal(x2,y2)

            Y12 = Y1 - Y2
            Y20 = Y2 - Y0
            Y01 = Y0 - Y1
            
            X_max = max(X0, X1, X2)

            EW0 = exp(U_ell*(X0-X_max)/(Pr*Vc))
            EW1 = exp(U_ell*(X1-X_max)/(Pr*Vc))
            EW2 = exp(U_ell*(X2-X_max)/(Pr*Vc))

            Delta_W   = EW0*Y12 + EW1*Y20 + EW2*Y01
            Delta_W_A = W[n0]*Y12 + W[n1]*Y20 + W[n2]*Y01
            Delta_W_B = W[n0]*(EW2-EW1) + W[n1]*(EW0-EW2) + W[n2]*(EW1-EW0)
            Delta_W_C = + W[n0]*(EW1*Y2 - EW2*Y1)  \
                        + W[n1]*(EW2*Y0 - EW0*Y2)  \
                        + W[n2]*(EW0*Y1 - EW1*Y0)
            
            A_W = Delta_W_A / Delta_W
            B_W = Delta_W_B / Delta_W
            C_W = Delta_W_C / Delta_W

        Psi_new[n_node] = (-Psi_BorderIntegral_nb + Psi_AreaIntegral)/Psi_BorderIntegral_a0
        

        segment_index = segment_indices[n_node]
        if segment_index in [1,2,3,4,6]:
            if abs(anbF) > 0 and abs(a0F) >0:
                Fi_new[n_node] = -anbF/a0F*QF + Fi[n_node]*(1-QF)
        elif segment_index in [5] :
            Fi_new[n_node] = H0 * nodes[n_node][1]
        else:
            print("keke")

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
                    
Plotter.PlotElements(triangulation, H)
Plotter.PlotNodes(nodes, Fi)
# Plotter.PlotScatter(nodes, Fi)
Plotter.PlotElements(triangulation, Mu)