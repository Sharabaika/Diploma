import matplotlib
import numpy as np
from MeshReader import readPoints
import Plotter
from math import exp, sqrt
import matplotlib.tri as tri

grid = open("meshes/square/grid.dat", "r")

borders = open("meshes/square/bounds.dat", "r")
bottom = open("meshes/square/bottom_bound.dat", "r")

nodes, triangles, segment_indices, neighbors = readPoints(grid, (borders, 1), (bottom, 2))

N_nodes = len(nodes)
N_trigs = len(triangles)

# PARAMS #
# ====== #

N_cyclies = 1000

Pr = 40.0
Vc = 1           # Viscosity


# Arrays #
# ====== #
Psi = np.zeros(N_nodes)
W = np.zeros(N_nodes)


# Init Arrays #
# =========== # 

# Dynamics #
# -------- #
for n_node in range(N_nodes):
    Psi[n_node] = np.random.rand(1)*0.001
    W[n_node] = np.random.rand(1)*0.001


Psi_new = np.array(Psi)
W_new = np.array(W)


for n_cycle in range(N_cyclies):
    if n_cycle % 50 == 0:
        print("cycle n = {n}".format(n=n_cycle))
    for n_node in range(N_nodes):
        Psi_BorderIntegral_a0 = 0
        Psi_BorderIntegral_nb = 0
        Psi_AreaIntegral = 0

        W_BorderIntegral = 0
        W_BorderIntegral_k0 = 0
        W_AreaIntegral = 0

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

            # /2 ?
            Psi_BorderIntegral_nb += (a1_psi1+a2_psi2)

            Psi_AreaIntegral += (22.*W[n0]+7.*W[n1]+7.*W[n2])*Delta/216.0

            # W #
            # - # 
            U_x = B_PSi
            U_y = -A_Psi
            U_ell = sqrt(U_x*U_x + U_y*U_y)

            sina = 0
            cosa = 1

            if U_ell > 0:
                sina = U_y / U_ell
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

            X12, Y12 = X1 - X2, Y1 - Y2
            X20, Y20 = X2 - X0, Y2 - Y0
            X01, Y01 = X0 - X1, Y0 - Y1


            Ya, Yb = (Y1 - Y0)/2, (Y2 - Y0)/2
            Xa, Xb = (X1 - X0)/2, (X2 - X0)/2

            
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

            W0, W1, W2 = W[n0], W[n1], W[n2]

            aBw = U_ell*Y12*Y0/ 8.0 / Pr + Vc*X12/2
            aCw = (U_ell*Y12)/(2*Pr)
            
            W_BorderIntegral += W1*(aBw*(EW0-EW2)+ aCw*(EW2*Y0-EW0*Y2)) \
                             +  W2*(aBw*(EW1-EW0)+ aCw*(EW0*Y1-EW1*Y0))

            W_BorderIntegral_k0 += aBw*(EW2-EW1)+ aCw*(EW1*Y2-EW2*Y1)

        segment_index = segment_indices[n_node]
        if segment_index in [0]:
            Psi_new[n_node] = (-Psi_BorderIntegral_nb + Psi_AreaIntegral)/Psi_BorderIntegral_a0
            W_new[n_node] = (W_AreaIntegral-W_BorderIntegral)/W_BorderIntegral_k0
        elif segment_index in [1,2]:
            Psi_new[n_node] = 0
            W_new[n_node] = (W_AreaIntegral-W_BorderIntegral)/W_BorderIntegral_k0      
        


    Psi = Psi_new
    W = W_new


Plotter.PlotNodes(nodes, Psi)
Plotter.PlotNodes(nodes, W)