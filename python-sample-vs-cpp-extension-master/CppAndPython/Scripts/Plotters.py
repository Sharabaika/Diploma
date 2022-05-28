from matplotlib import pyplot as plt
import matplotlib
from Scripts.ResultAnalysis import DynamycsAnalysis, MagneticsAnalysis
import numpy as np

def PlotMesh(points, triangles, segment_idex, index_nodes = False, scatted_nodes = False, index_regions = False):
    x, y = points[:, 0], points[:, 1]

    w, h = matplotlib.rcParams["figure.figsize"] 
    fig, ax = plt.subplots(figsize=(w*2.7, h*2.7))
    

    if scatted_nodes:
        m = ax.scatter(x, y, s=100, c=segment_idex, cmap='Dark2')   
        plt.colorbar(m)
    else:        
        ax.triplot(x, y, triangles, color='green')
        ax.set_aspect('equal')
        
    if index_nodes:
        for point_index in range(len(x)):
            ax.text(x=x[point_index], y=y[point_index], s = point_index, color='red', fontsize=10)

    if index_regions:
        for point_index in range(len(x)):
            ax.text(x=x[point_index], y=y[point_index], s = segment_idex[point_index], color='red', fontsize=8)

    plt.xlim((-2.1,2.1))
    plt.ylim((-2.1,2.1))
    plt.show()

def CompairMeshes(*meshes):
    w, h = matplotlib.rcParams["figure.figsize"] 
    fig, axs = plt.subplots(1, 2, sharey = True, figsize=(w*2.5, h*1.25))

    for n, mesh in enumerate(meshes):
        nodes, triangles, segment_indices, trig_neighbors, node_neighbours, triangle_indeces = mesh

        x, y = nodes[:, 0], nodes[:, 1]

        axs[n].triplot(x, y, triangles, color='green')
        axs[n].set_aspect('equal')
        axs[n].set_xlim((-2.1,2.1))
        axs[n].set_ylim((-2.1,2.1))
    
    axs[0].set_title("250 точек на границе, 6900 в области 2", fontsize=30)
    axs[1].set_title("375 точек на границе, 14900 в области 2", fontsize=30)


    plt.show()


def SavePlot(path):
    plt.savefig(path)

def ShowPlot():
    plt.show()

def Clear():
    plt.clf()

def PlotNodes(triangulation, Fi, **kwargs):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    cf = ax.tricontourf(triangulation, Fi)
    
    ax.set(
        xlim = kwargs.get("xlim", (-16, 16)),
        ylim = kwargs.get("ylim", (-16, 16))
    )

    fig.colorbar(cf, ax=ax)

def PlotScatter(points, z, **kwargs):
    x, y = points[:, 0], points[:, 1]

    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    sc = ax.scatter(x, y, s=100, c=z) 

    ax.set(
        xlim = kwargs.get("xlim", (-16, 16)),
        ylim = kwargs.get("ylim", (-16, 16))
    )

    plt.colorbar(sc)


def PlotElements(triang, z, **kwargs):
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    tpc = ax1.tripcolor(triang, z, shading='flat')
    
    ax1.set(
        xlim = kwargs.get("xlim", (-16, 16)),
        ylim = kwargs.get("ylim", (-16, 16))
    )    
    
    fig1.colorbar(tpc)


def PlotLevel(x, y, F, **kwargs):
    import numpy as np
    import matplotlib.pyplot as plot
    import pylab

    title = kwargs.get("title", "")
    xrange =  kwargs.get("xrange", (-10,10))
    yrange =  kwargs.get("yrange", (-10,10))

    manual = kwargs.get("manual", False)

    levels = []
    if "nlevels" in kwargs:
        nlevels = kwargs["nlevels"]
        indecies = list(filter(lambda i: xrange[0]<=x[i]<=xrange[1] and yrange[0]<=y[i]<=yrange[1], range(len(F))))
        maxF = max(F[indecies])
        minF = min(F[indecies])

        step = (maxF-minF)/(nlevels+2)
        levels = np.arange(minF + step/2, maxF - step/2, step)
        

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    plot.title(title)

    pylab.xlim(xrange)
    pylab.ylim(yrange)

    plot.xlabel('X')
    plot.ylabel('Y')


    if len(levels):
        contours = plot.tricontour(x, y, F, levels = levels)
    else:
        contours = plot.tricontour(x, y, F)

    plot.clabel(contours, inline=1, fontsize=10, manual = manual)

    plot.show()

def PlotLevelNodes(nodes, F, **kwargs):
    x, y = nodes[:,0], nodes[:,1]
    PlotLevel(x,y,F, **kwargs)

def PlotLevelTriangles(nodes, triangles, F, **kwargs):
    XPoints, YPoints = [], []
    for triangle in triangles:
        points = [nodes[node] for node in triangle]
        x, y = sum(points)
        XPoints.append(x/3)
        YPoints.append(y/3)
    
    PlotLevel(XPoints, YPoints, F, **kwargs)

class MagneticsPlot:
    def __init__(self, analysis : MagneticsAnalysis):
        self.analysis = analysis
        
        nodes, triangles, tags, trig_neighbours, node_neighbours, trianlge_indices = analysis.GetMesh()

        self.x = nodes[:,0]
        self.y = nodes[:,1]
        triangulation = matplotlib.tri.Triangulation(self.x,self.y,triangles)
        self.triangulation = triangulation

        mask = [index != 2 for index in trianlge_indices]
        self.inner_triangulation = matplotlib.tri.Triangulation(self.x,self.y,triangles)
        self.inner_triangulation.set_mask(mask)


    def PlotFi(self, b_inner_only = True, **kwargs):
        if b_inner_only:
            PlotNodes(self.inner_triangulation, self.analysis.GetFi(), **kwargs, xlim = (-2,2), ylim=(-2,2))
        else:
            PlotNodes(self.triangulation, self.analysis.GetFi(), **kwargs)

    def FillAx(ax, fun, tri, color):
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        cmap = plt.cm.get_cmap(color)
        level = ax.tricontourf(tri, fun, cmap = cmap)
        contours = ax.tricontour(tri, fun, linewidths = 1, colors = "black")

        driver = make_axes_locatable(ax)
        cax = driver.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(level, cax=cax)

        return contours

    def PlotH(self, **kwargs):            
        w, h = matplotlib.rcParams["figure.figsize"] 
        fig, ax = plt.subplots(sharey = True, figsize=(w*1.5, h*1.5))

        plt.xlabel('X')
        plt.ylabel('Y')

        inner = MagneticsPlot.FillAx(ax,self.analysis.GetH_Nodes(), self.inner_triangulation, "viridis")

        ax.set_xlim((-2,2))
        ax.set_ylim((-2,2))
        ax.set_aspect('equal')

        plt.clabel(inner, inline=1, fontsize=10, manual = True)

        plt.show()

    def PlotH_Nodes(self, b_inner_only = True, **kwargs):
        if b_inner_only:
            PlotNodes(self.inner_triangulation, self.analysis.GetH_Nodes(), **kwargs, xlim = (-2,2), ylim=(-2,2))
        else:
            PlotNodes(self.triangulation, self.analysis.GetH_Nodes(), **kwargs)

    def PlotMu(self, **kwargs):
        PlotElements(self.inner_triangulation, self.analysis.GetMu(), **kwargs, xlim = (-2,2), ylim=(-2,2))


class DynamycsPlot:
    def __init__(self, analysis : DynamycsAnalysis):
        self.analysis = analysis
        
        nodes, triangles, tags, trig_neighbours, node_neighbours, trianlge_indices = analysis.GetMesh()

        self.x = nodes[:,0]
        self.y = nodes[:,1]

        mask = [index != 2 for index in trianlge_indices]

        triangulation = matplotlib.tri.Triangulation(self.x,self.y,triangles)
        triangulation.set_mask(mask)
        self.triangulation = triangulation

    def FillAx(ax, fun, color, tri):
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        cmap = plt.cm.get_cmap(color)
        level = ax.tricontourf(tri, fun, cmap = cmap)
        contours = ax.tricontour(tri.x, tri.y, fun, linewidths = 1, colors = "black")

        driver = make_axes_locatable(ax)
        cax = driver.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(level, cax=cax)

        return contours

    def PlotPsiT(self):

        fig, axs = plt.subplots(1,2)

        plt.xlabel('X')
        plt.ylabel('Y')

        # Psi
        psi_contours = DynamycsPlot.FillAx(axs[0], self.analysis.GetPsi(), "viridis", self.triangulation)
        axs[0].set_title("Psi")

        T_contours = DynamycsPlot.FillAx(axs[1], self.analysis.GetT(), "turbo", self.triangulation)
        axs[1].set_title("T")

        for ax in axs:
            ax.set_xlim((-2,2))
            ax.set_ylim((-2,2))
            ax.set_aspect('equal')

        plt.clabel(psi_contours, inline=1, fontsize=10, manual = False)
        plt.clabel(T_contours, inline=1, fontsize=10, manual = True)

        plt.show()

    def PlotPsiLevel(self):        
        PlotLevel(self.x, self.y, self.analysis.GetPsi(), nlevels = 10, xrange = (-2,2), yrange = (-2,2), manual = True)


    def PlotPsi(self, **kwargs):
        PlotNodes(self.triangulation, self.analysis.GetPsi(), **kwargs, xlim=(-2, 2), ylim=(-2,2))


    def PlotW(self, **kwargs):
        PlotNodes(self.triangulation, self.analysis.GetW(), **kwargs, xlim=(-2, 2), ylim=(-2,2))


    def PlotT(self, **kwargs):
        PlotNodes(self.triangulation, self.analysis.GetT(), **kwargs, xlim=(-2, 2), ylim=(-2,2))

def PlotPsiT_Table(*results):
    n = len(results)

    fig, axs = plt.subplots(n, 2)
    plt.xlabel('X')
    plt.ylabel('Y')

    for i, result in enumerate(results):
        plotter = DynamycsPlot(result)

        psi_contours = DynamycsPlot.FillAx(axs[i, 0], result.GetPsi(), "viridis", plotter.triangulation)
        axs[i, 0].set_ylabel(f"Ram = { result.GetParam('Ram') }")

        T_contours = DynamycsPlot.FillAx(axs[i, 1], result.GetT(), "turbo", plotter.triangulation)

    for ax in axs.flatten():
        ax.set_xlim((-2,2))
        ax.set_ylim((-2,2))
        ax.set_aspect('equal')

    axs[0, 0].set_title("Psi")
    axs[0, 1].set_title("T")        

    plt.show()

def PlotPsiTable(*results):
    n = len(results)

    fig, axs = plt.subplots(n, len(results[0]))
    plt.xlabel('X')
    plt.ylabel('Y')

    for i, row_result in enumerate(results):
        angle = row_result[0].path.split("/")[-2].split("_")[-1]
        axs[i, 0].set_ylabel(f"alpha = {angle}")

        for j, result in enumerate(row_result):
            plotter = DynamycsPlot(result)

            psi_contours = DynamycsPlot.FillAx(axs[i, j], result.GetPsi(), "viridis", plotter.triangulation)


            if i == 0:
                axs[0, j].set_title(f"Ram = { result.GetParam('Ram') }")
            

    for ax in axs.flatten():
        ax.set_xlim((-2,2))
        ax.set_ylim((-2,2))
        ax.set_aspect('equal')

    plt.show()

def ComapairH(*results):
    n = len(results)

    n_col = 2

    import math
    rows = math.ceil(n / n_col)

    fig, axs = plt.subplots(rows,n_col)
    plt.xlabel('X')
    plt.ylabel('Y')

    for i, result in enumerate(results):   
        print(i)
        plotter = MagneticsPlot(result)

        ax = axs.flatten()[i]

        angle = result.path.split("/")[-2].split("_")[-1]

        mag = MagneticsPlot.FillAx(ax,result.GetH_Nodes(), plotter.inner_triangulation, "viridis")
        # ax.set_title(f"alpha = {angle}")

    for ax in axs.flatten():
        ax.set_xlim((-2,2))
        ax.set_ylim((-2,2))
        ax.set_aspect('equal')

    axs.flatten()[0].set_ylabel("H")    

    plt.show()

def PlotMaxH(results):
    alphas = []
    max_H = []

    fig, axs = plt.subplots()
    plt.xlabel('alpha')
    plt.ylabel('H max')

    for i, result in enumerate(results):   
        print(i)
        plotter = MagneticsPlot(result)

        angle = result.path.split("/")[-2].split("_")[-1]
        alphas.append(angle)
        H = result.GetH_Nodes()
        max_H.append(max(H))

    axs.plot(alphas, max_H)
    axs.scatter(alphas, max_H)


    plt.show()

def CompairLocalNuselts(results):
    n = len(results)

    fig, axs = plt.subplots(len(results))
    plt.xlabel('X')
    plt.ylabel('Y')

    for i, result_row in enumerate(results):
        ram = result_row[0].GetParam("Ram")
        axs[i].set_ylabel(f"Ram = {ram}")

        for j, result in enumerate(result_row):
            fis, nus = result.CalculateLocalNulselt()


            axs[i].plot(fis, nus, linestyle = 'dashed' if j == 0 else 'solid')

    plt.show()