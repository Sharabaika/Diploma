import numpy as np

def GenerateCurvyMesh(N, n, R, dr):
    res = np.zeros((N,3))

    for i,fi in enumerate(np.linspace(0, 2*np.pi, N)):
        r = R+dr*np.sin(fi*n)
        res[i,0] = r*np.cos(fi)
        res[i,1] = r*np.sin(fi)

    res[-1] = res[0]
    return res


def main():
    import matplotlib.pyplot as plt
    
    points = GenerateCurvyMesh(1000, 5, 0.3, 0.1)

    x,y = points[:,0], points[:,1]

    plt.scatter(x, y)
    plt.show()


if __name__ == "__main__":
    main()