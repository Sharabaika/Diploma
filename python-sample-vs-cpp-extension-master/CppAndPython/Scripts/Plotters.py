from matplotlib import pyplot as plt
import matplotlib
from Scripts.ResultAnalysis import DynamycsAnalysis, MagneticsAnalysis


def PlotMesh(points, triangles, segment_idex, index_nodes = False, scatted_nodes = False, index_regions = False):
    x, y = points[:, 0], points[:, 1]

    fig, ax = plt.subplots()
    

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

        step = (maxF-minF)/nlevels
        levels = np.arange(minF, maxF + step/2, step)
        

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
        tri = self.triangulation
        tri.set_mask(mask)
        self.inner_triangulation = tri

    def PlotFi(self, b_inner_only = True, **kwargs):
        if b_inner_only:
            PlotNodes(self.inner_triangulation, self.analysis.GetFi(), **kwargs, xlim = (-2,2), ylim=(-2,2))
        else:
            PlotNodes(self.triangulation, self.analysis.GetFi(), **kwargs)

    def PlotH(self, b_inner_only = True, **kwargs):
        if b_inner_only:
            PlotElements(self.inner_triangulation, self.analysis.GetH(), **kwargs, xlim = (-2,2), ylim=(-2,2))
        else:
            PlotElements(self.triangulation, self.analysis.GetH(), **kwargs)

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


    def PlotPsiLevel(self):        
        PlotLevel(self.x, self.y, self.analysis.GetPsi(), nlevels = 30, xrange = (-2,2), yrange = (-2,2), manual = True)


    def PlotPsi(self, **kwargs):
        PlotNodes(self.triangulation, self.analysis.GetPsi(), **kwargs, xlim=(-2, 2), ylim=(-2,2))


    def PlotW(self, **kwargs):
        PlotNodes(self.triangulation, self.analysis.GetW(), **kwargs, xlim=(-2, 2), ylim=(-2,2))


    def PlotT(self, **kwargs):
        PlotNodes(self.triangulation, self.analysis.GetT(), **kwargs, xlim=(-2, 2), ylim=(-2,2))
