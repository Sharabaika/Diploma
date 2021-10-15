
def PlotMesh(points, triangles, segment_idex, index_nodes = False, scatted_nodes = False):
    import matplotlib.pyplot as plt

    x, y = points[:, 0], points[:, 1]
    plt.triplot(x, y, triangles, color='green')

    if scatted_nodes:
        plt.scatter(x, y, s=100, c=segment_idex)    

    if index_nodes:
        for point_index in range(len(x)):
            plt.text(x=x[point_index], y=y[point_index], s = point_index, color='red', fontsize=14)

    plt.show()

def PlotT(points, T):
    import matplotlib.pyplot as plt
    import numpy as np

    x, y = points[:, 0], points[:, 1]

    v = np.linspace(0, 1, 11, endpoint=True)
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    tcf = ax1.tricontourf(x, y, T, v)
    fig1.colorbar(tcf)

    plt.show()

def PlotFi(points, Fi):
    import matplotlib.pyplot as plt
    import numpy as np

    x, y = points[:, 0], points[:, 1]

    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    cf = ax.tricontourf(x,y,Fi)
    fig.colorbar(cf, ax=ax)
    plt.show()

def PlotH(triangles, H):
    pass