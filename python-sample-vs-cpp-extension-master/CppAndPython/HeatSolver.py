from types import MemberDescriptorType
import matplotlib
import numpy as np
from Scripts.MeshReader import ReadRaw, ReadSaved
from math import exp, sqrt
import matplotlib.tri as tri
from Scripts.ResultAnalysis import PlotMesh, PlotScatter 
import Scripts.ResultFileHandling as files

from time import perf_counter

from pstats import Stats, SortKey

def main():
    Saver = files.ResultSaving("T")

    USUAL_HEAT_MEDIUM_INDEX = 2000
    TEMP_BORDER_INDEX = 1000

    mesh_name = "square"

    nodes, triangles, segment_indices, trig_neighbors, node_neighbours = ReadSaved(f"SavedMeshes/{mesh_name}.dat")

    N_nodes = len(nodes)
    N_trigs = len(triangles)


    # PARAMS #
    # ====== #
    Re_T = 1
    Pr_T = 1
    Vc_T = 1 


    N_CYCLIES_MAX = 100
    MAX_DELTA_ERROR = 1e-5

    QT = 1

    Saver.AddParams(mesh_name = mesh_name, QT = QT)

    # Arrays #
    # ====== #
    T = np.zeros(N_nodes)
    T_new = np.array(T)

    for n_node in range(N_nodes):
        x,y = nodes[n_node]
        if 0.4<x<0.6 and 0.4<y<0.6:
            T_new[n_node] = 1.0

    Error = 2*MAX_DELTA_ERROR

    n_cycle = 0
    while n_cycle < N_CYCLIES_MAX and Error>=MAX_DELTA_ERROR:
        for n_node in range(N_nodes):
            segment_index = segment_indices[n_node]

            # USUAL MEDIUM #
            # ============ #
            if USUAL_HEAT_MEDIUM_INDEX in segment_index or TEMP_BORDER_INDEX in segment_index:
                T_BorderIntegral = 0
                T_BorderIntegral_k0 = 0
                T_AreaIntegral = 0

                for n_trig_neigbor in trig_neighbors[n_node]:
                    n0, n1, n2 = triangles[n_trig_neigbor]
                    if n_node == n1:
                        n0, n1, n2 =  n_node, n2, n0
                    elif n_node == n2:
                        n0, n1, n2 = n_node, n0, n1

                    T0, T1, T2 = T[n0], T[n1], T[n2] 

                    x0, y0 = nodes[n0]
                    x1, y1 = nodes[n1]
                    x2, y2 = nodes[n2]
                    x10, y01 = x1-x0, y0-y1
                    x21, y12 = x2-x1, y1-y2
                    x02, y20 = x0-x2, y2-y0
                    Delta = x10*y20-x02*y01

                    U_x = 0
                    U_y = -0.01
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

                    X_min = min(X0, X1, X2)
                    X_max = max(X0, X1, X2)

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


                # Add "q" term
                x0, y0 = nodes[n_node]
                if TEMP_BORDER_INDEX in segment_index:
                    border_neighbours = []
                    for neighbour in node_neighbours[n_node]:
                        if TEMP_BORDER_INDEX in segment_indices[neighbour]:
                            border_neighbours.append(neighbour)

                    for n_border_neighbour in border_neighbours:
                        x1, y1 = nodes[n_border_neighbour]
                        l = sqrt((x1-x0)**2 + (y1-y0)**2)*0.5

                        q = 1
                        T_BorderIntegral += l * q

                if 0.4<x0<0.6 and 0.4<y0<0.6  and False:
                    T_new[n_node] = 1.0
                else:
                    T_new[n_node] = (-T_BorderIntegral + T_AreaIntegral)/T_BorderIntegral_k0                
                

        # ERRORS # 
        # ------ #
        Delta_Ts_Error_Squared = sum((T- T_new)**2)/(QT*QT)/sum(T_new**2)

        Error = sqrt(Delta_Ts_Error_Squared)        

        if n_cycle % 50 == 0 or True:
            print(f"cycle n = {n_cycle}, dT == {Error:.2e}")

        Saver.logger.LogErrors(T = Error)

        # NEXT STEP #
        # --------- #
        T = T*(1-QT) + np.copy(T_new)*QT

        n_cycle += 1


    import matplotlib as matplot

    x, y = nodes[:,0], nodes[:,1]
    triangulation = matplot.tri.Triangulation(x,y,triangles)

    # Saver.SaveResults("SavedResults", "ReworkedRe1000", W = W, Psi = Psi)

    from Scripts.ResultAnalysis import PlotNodes
    # PlotNodes(triangulation, Psi)
    PlotNodes(triangulation, T)
    # PlotScatter(nodes, W)


if __name__ == "__main__":
    main()