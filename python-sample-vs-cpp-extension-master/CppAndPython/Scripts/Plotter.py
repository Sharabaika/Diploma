import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

from Scripts.MeshReader import ReadSaved

def PlotMesh(points, triangles, segment_idex, index_nodes = False, scatted_nodes = False, index_regions = False):
    x, y = points[:, 0], points[:, 1]

    fig, ax = plt.subplots()
    
    ax.triplot(x, y, triangles, color='green')
    ax.set_aspect('equal')
    
    if scatted_nodes:
        ax.scatter(x, y, s=100, c=[arr[0] for arr in segment_idex])     

    if index_nodes:
        for point_index in range(len(x)):
            ax.text(x=x[point_index], y=y[point_index], s = point_index, color='red', fontsize=10)

    if index_regions:
        for point_index in range(len(x)):
            ax.text(x=x[point_index], y=y[point_index], s = segment_idex[point_index], color='red', fontsize=8)

    x, y = points[0]
    ax.arrow(0, 0, x, y, head_width=0.005, head_length=0.01, fc='k', ec='k')

    plt.show()

def PlotNodes(triangulation, Fi):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    cf = ax.tricontourf(triangulation, Fi)
    fig.colorbar(cf, ax=ax)
    plt.show()

def PlotScatter(points, z):
    x, y = points[:, 0], points[:, 1]

    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    sc = ax.scatter(x, y, s=100, c=z) 
    plt.colorbar(sc)
    plt.show()

def PlotElements(triang, z):
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    tpc = ax1.tripcolor(triang, z, shading='flat')
    fig1.colorbar(tpc)
    plt.show()


class CoolPlots:
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
        CoolPlots.PlotLevel(x,y,F, **kwargs)

    def PlotLevelTriangles(nodes, triangles, F, **kwargs):
        XPoints, YPoints = [], []
        for triangle in triangles:
            points = [nodes[node] for node in triangle]
            x, y = sum(points)
            XPoints.append(x/3)
            YPoints.append(y/3)
        
        CoolPlots.PlotLevel(XPoints, YPoints, F, **kwargs)


class ResultAnalysis:
    def __init__(self, folder, result_name):
        self.path = os.path.join(folder, result_name)
        
        self.logs = None
        self.saved = None
        self.params = None

        self.loaded_mesh = None

    def GetLogs(self):
        if self.logs is None:
            logs_path = os.path.join(self.path, "logs.csv")
            self.logs = pd.read_csv(logs_path, index_col=0)
        return self.logs

    def GetSavedResults(self):
        if self.saved is None:
            saved_path = os.path.join(self.path, "saved.csv")
            self.saved = pd.read_csv(saved_path)
        return self.saved

    def GetParams(self):
        if self.params is None:
            params_path = os.path.join(self.path, "params.csv")
            self.params = pd.read_csv(params_path)
        return self.params

    def GetMesh(self):
        if self.loaded_mesh is None:
            mesh_name = self.GetParams().iloc[0]['mesh_name']
            self.loaded_mesh = ReadSaved(f"SavedMeshes/{mesh_name}.dat")
        return self.loaded_mesh

    def PlotErrors(self, *args, **kwargs):       
        xmin, xmax = kwargs.get("xrange", (0,-1))
        xmax = len(self.GetLogs().index) if xmax == -1 else xmax 

        traces = kwargs.get("traces", self.GetLogs().columns)

        ax = self.GetLogs()[traces][xmin:xmax].plot()
        ax.set_yscale("log")
        ax.set_ylim((0,1))

        plt.show()

class DynamycsAnalysis(ResultAnalysis):
    def __init__(self, folder, result_name):
        super().__init__(folder, result_name)

    def GetPsi(self):
        return self.GetSavedResults()["Psi"]

    def GetW(self):
        return self.GetSavedResults()["W"]
    
    def PlotPsi(self):
        nodes, triangles, segment_indices, trig_neighbors, node_neighbours = self.GetMesh()

        x, y = nodes[:,0], nodes[:,1]
        triangulation = matplotlib.tri.Triangulation(x,y,triangles)
        
        CoolPlots.PlotLevel(x, y, self.GetPsi(), nlevels = 30, xrange = (0,1), yrange = (0,1), manual = True)

    def PlotW(self):
        nodes, triangles, segment_indices, trig_neighbors, node_neighbours = self.GetMesh()
        x, y = nodes[:,0], nodes[:,1]
        triangulation = matplotlib.tri.Triangulation(x,y,triangles)

        PlotNodes(triangulation, self.GetW())

        

def main():
    pass

if __name__ == "__main__":
    main()

