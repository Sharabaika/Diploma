from types import MemberDescriptorType
import matplotlib
import numpy as np
from Scripts.MeshReader import ReadRaw, ReadSaved
from math import exp, sqrt
import matplotlib.tri as tri
from Scripts.Plotter import PlotMesh, PlotScatter 
import Scripts.ResultFileHandling as files

from time import perf_counter

from pstats import Stats, SortKey

def main():
    Saver = files.ResultSaving("W", "Psi")

    # Regions
    HEATCONVECTION_REGION_INDEX = 2000  # solve fluid and heat equations

    CONVECTION_REGION_INDEX = 2000      # "water like" medium  
    CONVECTION_BORDER_INDEX = 1000      # solid border for fluid equation

    CONSTANT_T_REGION_INDEX = 2300      # region with constant temperature
    CONSTANT_T_FLOW_BORDER_INDEX = 2400 # border with defined q

    mesh_name = "square"

    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadSaved(f"SavedMeshes/{mesh_name}.dat")

    N_nodes = len(nodes)
    N_trigs = len(triangles)


    # PARAMS #
    # ====== #
    Re = 300
    Pr = 1/Re
    Vc = 1 

    Re_T = 1
    Pr_T = 1
    Vc_T = 1 

    N_CYCLIES_MAX = 1000
    MAX_DELTA_ERROR = 1e-5

    Vx = 0

    QPsi = 0.6
    QW = 0.5
    QT = 1

    Saver.AddParams(mesh_name = mesh_name, Re = Re, QPsi = QPsi, QW = QW, QT = QT)


    # Arrays #
    # ====== #
    Psi = np.zeros(N_nodes)
    W = np.zeros(N_nodes)
    T = np.zeros(N_nodes)

    for n_node in range(N_nodes):
        tags = segment_indices[n_node]
        is_wall = CONVECTION_BORDER_INDEX in tags
        x, y = nodes[n_node]
        # W[n_node] = np.random.rand(1)*0.001+0.001 
        Psi[n_node] = 0 if is_wall else 0.0001*np.sin(x*np.pi)*np.sin(y*np.pi)
        T [n_node] = np.random.rand(1)*0.001+0.001


    Psi_new = np.array(Psi)
    W_new = np.array(W)
    T_new = np.array(T)

    Error = 2*MAX_DELTA_ERROR

    n_cycle = 0

    while n_cycle < N_CYCLIES_MAX and Error>=Error:
        for n_node in range(N_nodes):
            segment_index = segment_indices[n_node]
            
            # USUAL MEDIUM #
            # ============ #
            if CONVECTION_REGION_INDEX in segment_index :
                Psi_BorderIntegral_a0 = 0
                Psi_BorderIntegral_nb = 0
                Psi_AreaIntegral = 0

                W_BorderIntegral = 0
                W_BorderIntegral_k0 = 0
                W_AreaIntegral = 0

                T_BorderIntegral = 0
                T_BorderIntegral_k0 = 0
                T_AreaIntegral = 0

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

                    # Psi #
                    # --- #
                    Psi0, Psi1, Psi2 = Psi[n0], Psi[n1], Psi[n2]
                    W0, W1, W2 = W[n0], W[n1], W[n2]
                    T0, T1, T2 = T[n0], T[n1], T[n2]

                    Delta_PsiA = Psi0*y12 + Psi1*y20 + Psi2*y01
                    Delta_PsiB = Psi0*x21 + Psi1*x02 + Psi2*x10

                    A_Psi = Delta_PsiA / Delta
                    B_PSi = Delta_PsiB / Delta

                    a0      = (x21*x21+y12*y12)/Delta
                    a1_psi1 = Psi1 * (y20*y12+x02*x21) / Delta
                    a2_psi2 = Psi2 * (y01*y12+x10*x21) / Delta

                    Psi_BorderIntegral_a0 += a0
                    Psi_BorderIntegral_nb += (a1_psi1+a2_psi2)

                    Psi_AreaIntegral += (22.0*W0+7.0*W1+7.0*W2)*Delta/216.0

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
                    
                    # X0, Y0 = (x0-xc)*cosa + (y0-yc)*sina, -(x0-xc)*sina + (y0-yc)*cosa
                    # X1, Y1 = (x1-xc)*cosa + (y1-yc)*sina, -(x1-xc)*sina + (y1-yc)*cosa
                    # X2, Y2 = (x2-xc)*cosa + (y2-yc)*sina, -(x2-xc)*sina + (y2-yc)*cosa

                    X0, Y0 = ToLocal(x0,y0)
                    X1, Y1 = ToLocal(x1,y1)
                    X2, Y2 = ToLocal(x2,y2)

                    X12, Y12 = X1 - X2, Y1 - Y2
                    X20, Y20 = X2 - X0, Y2 - Y0
                    X01, Y01 = X0 - X1, Y0 - Y1


                    Ya, Yb = (Y1 + Y0)*0.5, (Y2 + Y0)*0.5
                    Xa, Xb = (X1 + X0)*0.5, (X2 + X0)*0.5

                    X_min = min(X0, X1, X2)
                    X_max = max(X0, X1, X2)

                    kw0, kw1, kw2 = 0, 0, 0

                    ReVel =  U_ell/(Pr*Vc)
                    vel = ReVel*(X_max-X_min)

                    aBw = U_ell*Y12*Y0/ (8.0 * Pr) + Vc*X12/2
                    aCw = (U_ell*Y12)/(2.0*Pr)

                    if vel<=1e-8:
                        DW=ReVel*(X0*Y12+X1*Y20+X2*Y01)

                        kw0=(aBw*ReVel*(X2-X1)+aCw*(-Y12+ReVel*((X1-X_max)*Y2-(X2-X_max)*Y1)))/DW
                        kw1=(aBw*ReVel*(X0-X2)+aCw*(-Y20+ReVel*((X2-X_max)*Y0-(X0-X_max)*Y2)))/DW
                        kw2=(aBw*ReVel*(X1-X0)+aCw*(-Y01+ReVel*((X0-X_max)*Y1-(X1-X_max)*Y0)))/DW
                    
                    else:
                        EW0 = exp(U_ell*(X0-X_max)/(Pr*Vc))
                        EW1 = exp(U_ell*(X1-X_max)/(Pr*Vc))
                        EW2 = exp(U_ell*(X2-X_max)/(Pr*Vc))

                        Delta_W   = EW0*Y12 + EW1*Y20 + EW2*Y01
                        one_over_Delta_W = 1 / Delta_W

                        kw0 = (aBw*(EW2-EW1)+ aCw*(EW1*Y2-EW2*Y1))*one_over_Delta_W
                        kw1 = (aBw*(EW0-EW2)+ aCw*(EW2*Y0-EW0*Y2))*one_over_Delta_W
                        kw2 = (aBw*(EW1-EW0)+ aCw*(EW0*Y1-EW1*Y0))*one_over_Delta_W

                    W_BorderIntegral    += W1*(kw1) +  W2*(kw2)
                    W_BorderIntegral_k0 += kw0

                    Delta_TA = T0*y12 + T1*y20 + T2*y01
                    A_T = Delta_TA / Delta
                    W_AreaIntegral += Delta / 6 * (-A_T)

                    kw0_T, kw1_T, kw2_T = 0, 0, 0

                    ReVel_T =  U_ell/(Pr_T*Vc_T)
                    vel_T = ReVel_T*(X_max-X_min)

                    aBw_T = U_ell*Y12*Y0/ (8.0 * Pr_T) + Vc_T*X12/2
                    aCw_T = (U_ell*Y12)/(2.0*Pr_T)

                    if vel_T<=1e-8:
                        DW_T=ReVel_T*(X0*Y12+X1*Y20+X2*Y01)

                        kw0_T=(aBw_T*ReVel_T*(X2-X1)+aCw_T*(-Y12+ReVel_T*((X1-X_max)*Y2-(X2-X_max)*Y1)))/DW_T
                        kw1_T=(aBw_T*ReVel_T*(X0-X2)+aCw_T*(-Y20+ReVel_T*((X2-X_max)*Y0-(X0-X_max)*Y2)))/DW_T
                        kw2_T=(aBw_T*ReVel_T*(X1-X0)+aCw_T*(-Y01+ReVel_T*((X0-X_max)*Y1-(X1-X_max)*Y0)))/DW_T
                    
                    else:
                        EW0_T = exp(U_ell*(X0-X_max)/(Pr_T*Vc_T))
                        EW1_T = exp(U_ell*(X1-X_max)/(Pr_T*Vc_T))
                        EW2_T = exp(U_ell*(X2-X_max)/(Pr_T*Vc_T))

                        Delta_W_T   = EW0_T*Y12 + EW1_T*Y20 + EW2_T*Y01
                        one_over_Delta_W_T = 1 / Delta_W_T

                        kw0_T = (aBw_T*(EW2_T-EW1_T)+ aCw_T*(EW1_T*Y2-EW2_T*Y1))*one_over_Delta_W_T
                        kw1_T = (aBw_T*(EW0_T-EW2_T)+ aCw_T*(EW2_T*Y0-EW0_T*Y2))*one_over_Delta_W_T
                        kw2_T = (aBw_T*(EW1_T-EW0_T)+ aCw_T*(EW0_T*Y1-EW1_T*Y0))*one_over_Delta_W_T

                        


                    T_BorderIntegral    += T1*(kw1_T) +  T2*(kw2_T)
                    T_BorderIntegral_k0 += kw0_T

                Psi_new[n_node] = (-Psi_BorderIntegral_nb + Psi_AreaIntegral)/Psi_BorderIntegral_a0
                W_new[n_node] = (-W_BorderIntegral + W_AreaIntegral)/W_BorderIntegral_k0

                T_new[n_node] = (-T_BorderIntegral + T_AreaIntegral)/T_BorderIntegral_k0

            # SOLID BORDER #
            # ============ #
            elif CONVECTION_BORDER_INDEX in segment_index:            
                # normal #
                # ------ #
                node = nodes[n_node]
                sina, cosa = 0, 0

                normalX, normalY = 0, 0
                if 10 in segment_index:
                    # left  
                    sina, cosa = -1, 0
                if 11 in segment_index:
                    # right
                    sina, cosa = 1, 0
                if 12 in segment_index:
                    # bottom  
                    sina, cosa = 0, 1
                if 13 in segment_index:
                    # upp  
                    sina, cosa = 0, -1

                # LOCAL COORDS #
                # ------------ #

                x_center, y_center = nodes[n_node]
                def ToLocal_Border(x, y):
                    X = (x-x_center)*cosa + (y-y_center)*sina
                    Y = -(x-x_center)*sina + (y-y_center)*cosa
                    return X,Y


                # Integral values #
                # --------------- #
                # Just an area of triangle, nothing more
                W_Source_Integral = 0

                # Integral of W over triangle
                W_Source_Area_Integral = 0

                # Integral a to b
                W_Border_Integral = 0

                T_BorderIntegral = 0
                T_BorderIntegral_k0 = 0
                T_AreaIntegral = 0

                for n_trig_neigbor in trig_neighbors[n_node]:
                    triangle = triangles[n_trig_neigbor]

                    n0, n1, n2 = triangle
                    if n_node == n1:
                        n0, n1, n2 =  n_node, n2, n0
                    elif n_node == n2:
                        n0, n1, n2 = n_node, n0, n1

                    X0, Y0 = ToLocal_Border(*nodes[n0])
                    X1, Y1 = ToLocal_Border(*nodes[n1])
                    X2, Y2 = ToLocal_Border(*nodes[n2])

                    Psi0, Psi1, Psi2 = Psi[n0], Psi[n1], Psi[n2]

                    # Interpolation # 
                    # ------------- #
                    is_moving_border = 13 in segment_index

                    c1 = Psi1-Psi0-Vx*Y1 if is_moving_border else Psi1-Psi0
                    c2 = Psi2-Psi0-Vx*Y2 if is_moving_border else Psi2-Psi0

                    delta = 0.5*Y1*Y2*(X1*Y2 - X2*Y1)
                    delta_a = 0.5*(c1*(Y2**2) - c2*(Y1**2))
                    delta_b = c2*X1*Y1- c1*X2*Y2

                    a = 0
                    b = 0
                    if delta != 0:
                        a = delta_a / delta
                        b = delta_b / delta
                    else:
                        if Y1 == 0 and Y2 == 0:
                            print(f"exception in node {n_node}")
                        if Y1 == 0:
                            b = 2*c2/(Y2**2)
                        elif Y2 == 0:
                            b = 2*c1/(Y1**2)

                    xa, ya = X1*0.5, Y1*0.5
                    xb, yb = X2*0.5, Y2*0.5
                    xc, yc = (X1+X2)/3.0, (Y1+Y2)/3.0 

                    moving_border_add = Vx*(xb-xa) if is_moving_border else 0
                    
                    border_part = 0.5*a*(yb*yb-ya*ya) - 0.5*a*(xb*xb-xa*xa) - \
                                    b*(0.5*(ya+yc)*(xc-xa) + 0.5*(yb+yc)*(xb-xc))  + \
                                    moving_border_add
                    W_Border_Integral += border_part

                    # TODO abs???
                    x0, y0 = nodes[n0]
                    x1, y1 = nodes[n1]
                    x2, y2 = nodes[n2]
                    triangle_delta = abs((x1-x0)*(y2-y0) - (x0-x2)*(y0-y1))

                    W1, W2 = W[n1], W[n2]

                    W_Source_Area_Integral += 11.0*triangle_delta/108.0
                    W_Source_Integral += (7.0*W1+ 7.0*W2)*triangle_delta/216.0


                    x0, y0 = nodes[n0]
                    x1, y1 = nodes[n1]
                    x2, y2 = nodes[n2]
                    x10, y01 = x1-x0, y0-y1
                    x21, y12 = x2-x1, y1-y2
                    x02, y20 = x0-x2, y2-y0
                    Delta = x10*y20-x02*y01

                    Delta_PsiA = Psi0*y12 + Psi1*y20 + Psi2*y01
                    Delta_PsiB = Psi0*x21 + Psi1*x02 + Psi2*x10

                    A_Psi = Delta_PsiA / Delta
                    B_PSi = Delta_PsiB / Delta

                    U_x = B_PSi
                    U_y = -A_Psi
                    U_ell = sqrt(U_x*U_x + U_y*U_y)

                    sina_T = 0
                    cosa_T = 1

                    if U_ell > 0:
                        sina_T = U_y / U_ell
                        cosa_T = U_x / U_ell

                    xc = (x0+x1+x2)/3
                    yc = (y0+y1+y2)/3

                    def ToLocal_T(x, y):
                        X =  (x-xc)*cosa_T + (y-yc)*sina_T
                        Y = -(x-xc)*sina_T + (y-yc)*cosa_T
                        return X,Y

                    X0, Y0 = ToLocal_T(x0,y0)
                    X1, Y1 = ToLocal_T(x1,y1)
                    X2, Y2 = ToLocal_T(x2,y2)

                    kw0_T, kw1_T, kw2_T = 0, 0, 0

                    X_min = min(X0, X1, X2)
                    X_max = max(X0, X1, X2)

                    X12, Y12 = X1 - X2, Y1 - Y2
                    X20, Y20 = X2 - X0, Y2 - Y0
                    X01, Y01 = X0 - X1, Y0 - Y1

                    ReVel_T =  U_ell/(Pr_T*Vc_T)
                    vel_T = ReVel_T*(X_max-X_min)

                    aBw_T = U_ell*Y12*Y0/ (8.0 * Pr_T) + Vc_T*X12/2
                    aCw_T = (U_ell*Y12)/(2.0*Pr_T)

                    if vel_T<=1e-8:
                        DW_T=ReVel_T*(X0*Y12+X1*Y20+X2*Y01)

                        kw0_T=(aBw_T*ReVel_T*(X2-X1)+aCw_T*(-Y12+ReVel_T*((X1-X_max)*Y2-(X2-X_max)*Y1)))/DW_T
                        kw1_T=(aBw_T*ReVel_T*(X0-X2)+aCw_T*(-Y20+ReVel_T*((X2-X_max)*Y0-(X0-X_max)*Y2)))/DW_T
                        kw2_T=(aBw_T*ReVel_T*(X1-X0)+aCw_T*(-Y01+ReVel_T*((X0-X_max)*Y1-(X1-X_max)*Y0)))/DW_T
                    
                    else:
                        EW0_T = exp(U_ell*(X0-X_max)/(Pr_T*Vc_T))
                        EW1_T = exp(U_ell*(X1-X_max)/(Pr_T*Vc_T))
                        EW2_T = exp(U_ell*(X2-X_max)/(Pr_T*Vc_T))

                        Delta_W_T   = EW0_T*Y12 + EW1_T*Y20 + EW2_T*Y01
                        one_over_Delta_W_T = 1 / Delta_W_T

                        kw0_T = (aBw_T*(EW2_T-EW1_T)+ aCw_T*(EW1_T*Y2-EW2_T*Y1))*one_over_Delta_W_T
                        kw1_T = (aBw_T*(EW0_T-EW2_T)+ aCw_T*(EW2_T*Y0-EW0_T*Y2))*one_over_Delta_W_T
                        kw2_T = (aBw_T*(EW1_T-EW0_T)+ aCw_T*(EW0_T*Y1-EW1_T*Y0))*one_over_Delta_W_T

                        
                    T0, T1, T2 = T[n0], T[n1], T[n2]

                    T_BorderIntegral    += T1*(kw1_T) +  T2*(kw2_T)
                    T_BorderIntegral_k0 += kw0_T

                Psi_new[n_node] = 0

                x, y = node
                if x == 0 and y == 0 or \
                x == 0 and y == 1 or \
                x == 1 and y == 1 or \
                x == 1 and y == 0:
                    W_new[n_node] = 0 
                else:
                    W_new[n_node] = (-W_Border_Integral + W_Source_Integral)/W_Source_Area_Integral

                # Add "q" term
                x0, y0 = nodes[n_node]
                if CONSTANT_T_FLOW_BORDER_INDEX in segment_index or True:
                    border_neighbours = []
                    for neighbour in node_neighbours[n_node]:
                        if CONSTANT_T_FLOW_BORDER_INDEX in segment_indices[neighbour]:
                            border_neighbours.append(neighbour)

                    for n_border_neighbour in border_neighbours:
                        x1, y1 = nodes[n_border_neighbour]
                        l = sqrt((x1-x0)**2 + (y1-y0)**2)*0.5

                        q = 0
                        T_BorderIntegral += l * q

                T_new[n_node] = (-T_BorderIntegral + T_AreaIntegral)/T_BorderIntegral_k0   
                if 12 in segment_index and x0 > 0.5:
                    T_new[n_node] = 1             


        # ERRORS # 
        # ------ #
        Delta_Psi_Error_Squared = sum((Psi - Psi_new)**2)/(QPsi*QPsi)/sum(Psi_new**2)
        Delta_Ws_Error_Squared = sum((W - W_new)**2)/(QW*QW)/sum(W_new**2)

        Delta_Ts_Error_Squared = sum((T- T_new)**2)/(QT*QT)/sum(T_new**2)

        Error = max(Delta_Psi_Error_Squared, Delta_Ws_Error_Squared)
        

        if n_cycle % 50 == 0 or True:
            print(f"cycle n = {n_cycle}, dpsi == {sqrt(Delta_Psi_Error_Squared):.2e}, dW = {sqrt(Delta_Ws_Error_Squared):.2e}, dT = {sqrt(Delta_Ts_Error_Squared):.2e}")

        Saver.logger.LogErrors(Psi = sqrt(Delta_Psi_Error_Squared), W = sqrt(Delta_Ws_Error_Squared))

        # NEXT STEP #
        # --------- #
        Psi = Psi*(1-QPsi) + np.copy(Psi_new)*QPsi
        W = W*(1-QW) + np.copy(W_new)*QW

        T = T*(1-QT) + np.copy(T_new)*QT

        n_cycle += 1



    import matplotlib as matplot

    x, y = nodes[:,0], nodes[:,1]
    triangulation = matplot.tri.Triangulation(x,y,triangles)

    # Saver.SaveResults("SavedResults", "ReworkedRe1000", W = W, Psi = Psi)

    from Scripts.Plotter import PlotNodes
    PlotNodes(triangulation, T)
    PlotNodes(triangulation, Psi)

    PlotNodes(triangulation, W)
    # PlotScatter(nodes, W)


if __name__ == "__main__":
    main()