import numpy as np
from pytools import delta
from Scripts.MeshReader import ReadRaw, ReadSaved
from math import atan2, exp, sqrt
import matplotlib.tri as tri
from Scripts.ResultAnalysis import MagneticsAnalysis, PlotElements, PlotMesh, PlotNodes, PlotScatter 
import Scripts.ResultFileHandling as files
import matplotlib as matplot

ONE_THIRD = 1.0 / 3.0
ELEVEN_OVER_108 = 11.0 / 108.0

def solve(*args, **kwargs):
    Saver = files.ResultSaving("W", "Psi", "T")

    # Mesh data #
    # ========= #
    mesh_name = kwargs.get("mesh_name", "N120_n4_R1_dr0.3_extended")
    result_name =  kwargs.get("result_name", f"saved_result_fluid_and_magneticsV0_{mesh_name}")

    nodes, triangles, segment_indices, trig_neighbors, node_neighbours, triangle_indeces = ReadSaved(f"SavedMeshes/{mesh_name}.dat")

    x_coords, y_coords = nodes[:,0], nodes[:,1]
    triangulation = matplot.tri.Triangulation(x_coords, y_coords, triangles)

    N_nodes = len(nodes)
    N_trigs = len(triangles)

    magnetics_result_name = f"magnetic_test_finall_{mesh_name}"
    magnetics_result = MagneticsAnalysis("SavedResults", magnetics_result_name)

    magnetics_result_mesh_name = magnetics_result.GetMeshName()
    if (magnetics_result_mesh_name != mesh_name):
        raise Exception("magnetics results are made for different mesh")


    # Regions #
    # ======= #
    CONDUCTOR_REGION_INDEX = 0
    CONDUCTOR_BORDER_INDEX = 1
    MEDIUM_REGION_INDEX = 2
    MEDIUM_OUTER_BORDER_INDEX = 3
    VOID_REGION_INDEX = 4
    VOID_OUTER_BORDER_INDEX = 5
    # 5 444 3 222 1 000 1 222 3 444 5

    # TODO replace ?
    is_a_fluid_region = lambda node_index : segment_indices[node_index] in [CONDUCTOR_BORDER_INDEX, MEDIUM_REGION_INDEX, MEDIUM_OUTER_BORDER_INDEX]
    is_a_wall = lambda node_index : segment_indices[node_index] in [CONDUCTOR_BORDER_INDEX, MEDIUM_OUTER_BORDER_INDEX]

    is_a_fluid_region_array = [is_a_fluid_region(n_node) for n_node in range(N_nodes)]
    is_a_wall_array = [is_a_wall(n_node) for n_node in range(N_nodes)]
    
    fluid_domain_nodes_indeces_array = list(filter(is_a_fluid_region, range(N_nodes)))

    # PARAMS #
    # ====== #
    # Dynamics
    Pr = kwargs.get("Pr", 10)
    Ra = kwargs.get("Ra", 30000)
    Ram = kwargs.get("Ram", 1)
    chi0 = magnetics_result.GetParam("chi0")

    # Temperature
    Re_T = 1
    Pr_T = 1 # ?
    Vc_T = 1 

    T_inner = 1
    T_outer = 0

    # Cycles
    N_CYCLIES_MAX = kwargs.get("N_CYCLIES_MAX", 2000)
    MAX_DELTA_ERROR = kwargs.get("MAX_DELTA_ERROR", 1e-5)

    PRINT_LOG_EVERY_N_CYCLES = 10

    # Unused, wall velocity
    Vx = 0 

    # Relaxation
    QW = kwargs.get("QW", 0.7)
    QPsi = kwargs.get("QPsi", 1.2)
    QT = kwargs.get("QT", 1.2)

    Saver.AddParams(mesh_name = mesh_name, Ra = Ra, Ram=Ram, magnetics_result_name = magnetics_result_name, chi0=chi0, Pr = Pr, QPsi = QPsi, QW = QW, QT = QT)

    # Arrays #
    # ====== #
    Psi = np.zeros(N_nodes)
    W = np.zeros(N_nodes)
    T = np.zeros(N_nodes)

    H_nodes = magnetics_result.GetH_Nodes()
    H_triangles = magnetics_result.GetH()
    mu_triangles = magnetics_result.GetMu()

    dHdx_triangles = np.zeros(N_trigs)
    dHdy_triangles = np.zeros(N_trigs)

    # Init
    for n_node in range(N_nodes):
        if is_a_fluid_region_array[n_node]:
            tag = segment_indices[n_node]
            x, y = nodes[n_node]
            r = sqrt(x*x+y*y)

            W[n_node] = np.random.rand(1)*0.01+0.01 
            Psi[n_node] = 0.0001*np.sin(3*x*np.pi)*np.sin(3*y*np.pi)
            T[n_node] = T_inner + (T_outer-T_inner)*(r-1.0)

            if tag == CONDUCTOR_BORDER_INDEX:
                Psi[n_node] = 0
                T[n_node] = T_inner
            elif tag == MEDIUM_OUTER_BORDER_INDEX:
                Psi[n_node] = 0
                T[n_node] = T_outer
        else:
            Psi[n_node] = 0.0
            W[n_node] = 0.0
            T[n_node] = 0.0


    for n_trig in range(N_trigs):
        n0, n1, n2 = triangles[n_trig]

        x0, y0 = nodes[n0]
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        x10, y01 = x1-x0, y0-y1
        x21, y12 = x2-x1, y1-y2
        x02, y20 = x0-x2, y2-y0

        H0, H1, H2 = H_nodes[n0], H_nodes[n1], H_nodes[n2]

        AH = H0*y12 + H1*y20 + H2*y01
        BH = H0*x21 + H1*x02 + H2*x10
        
        dHdx_triangles[n_trig] = AH
        dHdy_triangles[n_trig] = BH


    Psi_new = np.array(Psi)
    W_new = np.array(W)
    T_new = np.array(T)

    Psi_errors = np.zeros(N_nodes)
    W_errors = np.zeros(N_nodes)
    T_errors = np.zeros(N_nodes)

    # Cycles #
    # ====== #
    Error = 2*MAX_DELTA_ERROR

    n_cycle = 0 
    while n_cycle < N_CYCLIES_MAX and Error>=MAX_DELTA_ERROR:
        for n_node in fluid_domain_nodes_indeces_array:
            aPsi0 = 0
            aPsinb = 0
            aW0 = 0
            aWnb = 0
            aT0 = 0
            aTnb = 0
            S = 0         
            I = 0
            source_integral = 0

            segment_index = segment_indices[n_node]
            for n_trig_neigbor in trig_neighbors[n_node]:
                # Neighbour triangles #
                # =================== #
                if not all([is_a_fluid_region_array[n_node_local] for n_node_local in triangles[n_trig_neigbor]]):
                    continue

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

                x01, y10 = -x10,-y01
                x12, y21 = -x21, -y12
                x20, y02 = -x02, -y20

                Delta = x10*y20-x02*y01

                Psi0, Psi1, Psi2 = Psi[n0], Psi[n1], Psi[n2]
                W0, W1, W2 = W[n0], W[n1], W[n2]
                T0, T1, T2 = T[n0], T[n1], T[n2]

                if is_a_wall_array[n_node]:
                    # Boundaries # 
                    # ========== #

                    # Local Coords #
                    # ============ #

                    # mormal
                    border_neighbours = []
                    for neighbour in node_neighbours[n_node]:
                        if is_a_wall_array[neighbour]:
                            border_neighbours.append(neighbour)

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

                    if segment_index == MEDIUM_OUTER_BORDER_INDEX:
                        normalX, normalY = -normalX, -normalY

                    sina, cosa = -normalX, normalY                
                    
                    # Local coords
                    X1=(x1-x0)*cosa+(y1-y0)*sina
                    Y1=(y1-y0)*cosa-(x1-x0)*sina

                    X2=(x2-x0)*cosa+(y2-y0)*sina
                    Y2=(y2-y0)*cosa-(x2-x0)*sina

                    a, b = 0, 0
                    if Y1==0: 
                        a=0
                        b=2*(Psi2-Psi0)/Y2**2
                    elif Y2==0:
                        a=0
                        b=2*(Psi1-Psi0)/Y1**2
                    else:
                        DeltaGr = 0.5 * Y1*Y2*(X1*Y2-X2*Y1)
                        DeltaAGr = 0.5 * ( (Psi1-Psi0) * Y2**2 - (Psi2-Psi0) * Y1**2 )
                        DeltaBGr = (Psi2-Psi0)*X1*Y1 - (Psi1-Psi0)*X2*Y2

                        a=DeltaAGr/DeltaGr
                        b=DeltaBGr/DeltaGr
                    
                    Xa = X1*0.5; Ya = Y1*0.5; Xb = X2*0.5; Yb = Y2*0.5
                    Xc = (X1+X2) * ONE_THIRD; Yc = (Y1+Y2) * ONE_THIRD

                    I = I + 0.5*a*((Yb**2-Ya**2)-(Xb**2-Xa**2))
                    I = I - 0.5*b*((Ya+Yc)*(Xc-Xa)+(Yb+Yc)*(Xb-Xc))

                    aWnb = aWnb + (7.0*W1+7.0*W2)/216.0*Delta
                    aW0 = aW0 + 11.0 * Delta / 108.0

                else:
                    # Inner nodes #
                    # =========== #

                    # Stream function #
                    # =============== #
                    aPsi0 = aPsi0 + (x21**2+y12**2)/Delta
                    aPsinb = aPsinb + ( Psi1*(y20*y12+x02*x21) + Psi2*(y01*y12+x10*x21) )/Delta
                    S = S + (22.0*W0 + 7.0*W1 + 7.0*W2) * Delta/216.0

                    # Stream local coords #
                    # =================== #
                    Apsi = Psi0*y12 + Psi1*y20 + Psi2*y01
                    Bpsi = Psi0*x21 + Psi1*x02 + Psi2*x10

                    Vx = Bpsi/Delta
                    Vy = -Apsi/Delta
                    Vel=sqrt(Vx**2 + Vy**2)

                    sina = 0
                    cosa = 1

                    if Vel>1e-8:
                        cosa = Vx/Vel
                        sina = Vy/Vel

                    xc = (x0+x1+x2)*ONE_THIRD
                    yc = (y0+y1+y2)*ONE_THIRD  

                    X0=(x0-xc)*cosa+(y0-yc)*sina
                    Y0=(y0-yc)*cosa-(x0-xc)*sina
                    X1=(x1-xc)*cosa+(y1-yc)*sina
                    Y1=(y1-yc)*cosa-(x1-xc)*sina
                    X2=(x2-xc)*cosa+(y2-yc)*sina
                    Y2=(y2-yc)*cosa-(x2-xc)*sina
                    
                    Xmax=max([X0,X1,X2]); Xmin=min([X0,X1,X2])
                    
                    X12 = X1-X2; Y12 = Y1-Y2
                    Y20 = Y2-Y0; Y01 = Y0-Y1  
                    X10 = X1-X0; X02 = X0-X2


                    # Temperature #
                    # =========== #
                    abT = Vel*Y12*Y0/8.0+X12/2.0
                    acT = Vel*Y12/2.0

                    k0T=0
                    k1T=0
                    k2T=0
                    
                    if abs(Vel)<1e-14:
                        k0T = 0
                        k1T = 0
                        k2T = 0
                    elif abs(Vel*(Xmax-Xmin))<1e-8:
                        DT = Vel*(X0*Y12+X1*Y20+X2*Y01)
                        k0T = (abT*Vel*(X2-X1)+acT*(-Y12+Vel*((X1-Xmax)*Y2-(X2-Xmax)*Y1)))/DT
                        k1T = (abT*Vel*(X0-X2)+acT*(-Y20+Vel*((X2-Xmax)*Y0-(X0-Xmax)*Y2)))/DT
                        k2T = (abT*Vel*(X1-X0)+acT*(-Y01+Vel*((X0-Xmax)*Y1-(X1-Xmax)*Y0)))/DT
                    else:
                        E0T = exp(Vel*(X0-Xmax))
                        E1T = exp(Vel*(X1-Xmax))
                        E2T = exp(Vel*(X2-Xmax))
                        DT = E0T*Y12+E1T*Y20+E2T*Y01
                        k0T = (abT*(E2T-E1T)+acT*(E1T*Y2-E2T*Y1))/DT
                        k1T = (abT*(E0T-E2T)+acT*(E2T*Y0-E0T*Y2))/DT
                        k2T = (abT*(E1T-E0T)+acT*(E0T*Y1-E1T*Y0))/DT

                    AT = T0*y12 + T1*y20 + T2*y01
                    BT = T0*x21 + T1*x02 + T2*x10

                    c_mag = dHdx_triangles[n_trig_neigbor] * BT - dHdy_triangles[n_trig_neigbor] * AT
                    source_integral += Ra*AT/6.0 + Ram*mu_triangles[n_trig_neigbor]*(c_mag)/6.0
                    
                    aT0 = aT0 + k0T
                    aTnb = aTnb + k1T*T1 + k2T*T2


                    # Vorticity #
                    # ========= #
                    VelPr = Vel/Pr
                    abW = VelPr*Y12*Y0/8.0+X12/2.0
                    acW = VelPr*Y12/2.0

                    k0W=0
                    k1W=0
                    k2W=0

                    if abs(Vel)<1e-14:
                        k0W = 0
                        k1W = 0
                        k2W = 0
                    elif abs(VelPr*(Xmax-Xmin))<1e-8:
                            DW = VelPr*(X0*Y12+X1*Y20+X2*Y01)
                            k0W = (abW*VelPr*(X2-X1)+acW*(-Y12+VelPr*((X1-Xmax)*Y2-(X2-Xmax)*Y1)))/DW
                            k1W = (abW*VelPr*(X0-X2)+acW*(-Y20+VelPr*((X2-Xmax)*Y0-(X0-Xmax)*Y2)))/DW
                            k2W = (abW*VelPr*(X1-X0)+acW*(-Y01+VelPr*((X0-Xmax)*Y1-(X1-Xmax)*Y0)))/DW
                    else:
                            E0 = exp(VelPr*(X0-Xmax))
                            E1 = exp(VelPr*(X1-Xmax))
                            E2 = exp(VelPr*(X2-Xmax))
                            DW = (E0*Y12+E1*Y20+E2*Y01)
                            k0W = (abW*(E2-E1)+acW*(E1*Y2-E2*Y1))/DW
                            k1W = (abW*(E0-E2)+acW*(E2*Y0-E0*Y2))/DW
                            k2W = (abW*(E1-E0)+acW*(E0*Y1-E1*Y0))/DW

                    aW0 = aW0 + k0W
                    aWnb = aWnb + k1W*W1+k2W*W2
            
            # Advance #
            # ======= #
            if not is_a_wall_array[n_node]:
                # Medium #
                # ====== #
                Psi_new[n_node] = (-aPsinb+S)/aPsi0
                W_new[n_node] = -(aWnb+source_integral)/aW0
                T_new[n_node] = -aTnb/aT0
            else:
                # Wall #
                # ==== #
                Psi_new[n_node] = 0
                W_new[n_node] = -(I+aWnb)/aW0
                if segment_index == CONDUCTOR_BORDER_INDEX:
                    T_new[n_node] = T_inner
                else:
                    T_new[n_node] = T_outer
            
            Psi_errors[n_node] = (Psi[n_node] - Psi_new[n_node])
            W_errors[n_node] = (W_new[n_node] - W[n_node])
            T_errors[n_node] = (T_new[n_node] - T[n_node])


        Delta_Psi_Error = sqrt(sum(Psi_errors**2)/(QPsi*QPsi)/sum(Psi_new[is_a_fluid_region_array]**2))
        Delta_Ws_Error = sqrt(sum(W_errors**2)/(QW*QW)/sum(W_new[is_a_fluid_region_array]**2))

        Delta_Ts_Error = sqrt(sum(T_errors**2)/(QT*QT)/sum(T_new[is_a_fluid_region_array]**2))

        Error = max(Delta_Psi_Error, Delta_Ws_Error, Delta_Ts_Error)
        
        if n_cycle % PRINT_LOG_EVERY_N_CYCLES == 0:
            print(f"cycle n = {n_cycle}, dpsi == {(Delta_Psi_Error):.5e}, dW = {(Delta_Ws_Error):.5e}, dT = {(Delta_Ts_Error):.5e}")

        Saver.logger.LogErrors(Psi = (Delta_Psi_Error), W = (Delta_Ws_Error), T = (Delta_Ts_Error))

        # NEXT STEP #
        # --------- #
        Psi = Psi*(1-QPsi) + np.copy(Psi_new)*QPsi
        W = W*(1-QW) + np.copy(W_new)*QW

        T = T*(1-QT) + np.copy(T_new)*QT

        n_cycle += 1

    Saver.SaveResults("SavedResults", result_name)
    Saver.SaveResult("SavedResults", result_name, "nodes", W = W, Psi = Psi, T = T)

    # mask = [index != 2 for index in triangle_indeces]
    # triangulation.set_mask(mask)

    # PlotNodes(triangulation, T)
    # PlotNodes(triangulation, Psi)

    # PlotNodes(triangulation, W)

def main():
    ram_range = [1, 100, 10000, 100000, 150000, 200000, 250000]
    for ram in ram_range:    
        solve(Ram = ram, result_name = f"validation_ram_{ram}")

if __name__ == "__main__":
    main()

            
            

