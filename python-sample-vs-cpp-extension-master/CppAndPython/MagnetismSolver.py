import matplotlib
from matplotlib.pyplot import title
import numpy as np
from Scripts.MeshReader import ReadRaw, ReadSaved
import Scripts.ResultAnalysis as ResultAnalysis
from math import exp, sqrt
import matplotlib.tri as tri
import Scripts.ResultFileHandling as files

Saver = files.ResultSaving("Fi")

# Mesh data #
# ========= #
mesh_name = f"N120_n4_R1_dr0.3_extended"
result_name = f"magnetic_test_finall_{mesh_name}"

nodes, triangles, segment_indices, trig_neighbors, node_neighbours, trianlge_indices = ReadSaved(f"SavedMeshes/{mesh_name}.dat")

N_nodes = len(nodes)
N_trigs = len(triangles)


# Regions #
# ======= #
CONDUCTOR_REGION_INDEX = 0
CONDUCTOR_BORDER_INDEX = 1
MEDIUM_REGION_INDEX = 2
MEDIUM_OUTER_BORDER_INDEX = 3
VOID_REGION_INDEX = 4
VOID_OUTER_BORDER_INDEX = 5
# 5 444 3 222 1 000 1 222 3 444 5


# PARAMS #
# ====== #
# Relaxation
QF = 1.0

#Field
chi0 = 3.0
H0 = 1.0
mu0 = 1000

# Cycles
N_CYCLIES_MAX = 2000
MAX_DELTA_ERROR = 1e-5


Saver.AddParams(mesh_name = mesh_name, chi0 = chi0, H0 = H0, mu0 = mu0, QF = QF)

# Arrays #
# ====== #
H = np.zeros(N_trigs)
Mu = np.zeros(N_trigs)
Fi = np.zeros(N_nodes)  

H_nodes = np.zeros(N_nodes)

# Init
for n_triangle in range(N_trigs):
    H[n_triangle] = H0

    segment_index = trianlge_indices[n_triangle]
    if segment_index == 0:
        Mu[n_triangle] = mu0
    elif segment_index == 2:
        Mu[n_triangle] = 1 + chi0*H[n_triangle]/(1+chi0*H[n_triangle])
        # Mu[n_triangle] = 1
    elif segment_index == 4:
        Mu[n_triangle] = 1

H_new = np.copy(H)
Mu_new = np.copy(Mu)

Fi_new = np.array(Fi)
Fi_errors = np.zeros(N_nodes)

# Solve #
# ===== #
Error = 2*MAX_DELTA_ERROR

n_cycle = 0 
while n_cycle < N_CYCLIES_MAX and Error>=MAX_DELTA_ERROR:
    for n_node in range(N_nodes):

        a0F = 0
        anbF = 0

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

        bInfinity = segment_index == VOID_OUTER_BORDER_INDEX
        if not bInfinity:
            if abs(a0F) > 0:
                Fi_new[n_node] = -anbF/a0F
        else:
            Fi_new[n_node] = H0 * nodes[n_node][1]

        Fi_errors[n_node] = Fi_new[n_node]-Fi[n_node]

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
        segment_index = trianlge_indices[n_triangle]
        if segment_index == CONDUCTOR_REGION_INDEX:
            Mu_new[n_triangle] = mu0
        elif segment_index == MEDIUM_REGION_INDEX:
            Mu_new[n_triangle] = 1 + chi0*H_new[n_triangle]/(1+chi0*H_new[n_triangle])
            # Mu_new[n_triangle] = 1 
        elif segment_index == VOID_REGION_INDEX:
            Mu_new[n_triangle] = 1
        else:
            print("aboba")

    Mu = Mu_new
    H = H_new
    Fi = Fi*(1-QF) + np.copy(Fi_new)*QF

    Delta_Fi_Error = sqrt(sum(Fi_errors**2)/(QF*QF)/sum(Fi_new**2))
    Error = Delta_Fi_Error

    if n_cycle % 50 == 0 or True:
        print(f"cycle n = {n_cycle}, dFi = {(Delta_Fi_Error):.2e}")

    Saver.logger.LogErrors(Fi = (Delta_Fi_Error))

    n_cycle += 1

for n_node in range(N_nodes):
    numerator_sum = 0.0
    denominator_sum = 0.0

    for n_trig_neigbor in trig_neighbors[n_node]:
        n0, n1, n2 = triangles[n_trig_neigbor]
        if n_node == n1:
            n0, n1, n2 =  n_node, n2, n0
        elif n_node == n2:
            n0, n1, n2 = n_node, n0, n1

        x0, y0 = nodes[n0]
        x1, y1 = nodes[n1]
        x2, y2 = nodes[n2]

        Mx = (x1+x2)*0.5
        My = (y1+y2)*0.5

        dx = Mx-x0
        dy = My-y0 
        mr = sqrt(dx*dx+dy*dy)
        r = 2.0*mr/3.0

        numerator_sum += H[n_trig_neigbor]/r
        denominator_sum += 1.0/r
        
    H_nodes[n_node] = numerator_sum / denominator_sum


Saver.SaveResults("SavedMagnetics", result_name)
Saver.SaveResult("SavedMagnetics", result_name, "triangles",  H = H, Mu = Mu)
Saver.SaveResult("SavedMagnetics", result_name, "nodes", Fi = Fi, H_nodes = H_nodes)

x, y = nodes[:, 0], nodes[:, 1]
triangulation = tri.Triangulation(x, y, triangles) 

ResultAnalysis.PlotNodes(triangulation, Fi)
ResultAnalysis.PlotElements(triangulation, H)
ResultAnalysis.PlotElements(triangulation, Mu)