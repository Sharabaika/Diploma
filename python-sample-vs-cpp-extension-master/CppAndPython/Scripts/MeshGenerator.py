import numpy as np

def GenerateCurvyMesh(N, n, R, dr):
    res = np.zeros((N,3))

    for i,fi in enumerate(np.linspace(0, 2*np.pi, N)):
        r = R+dr*np.cos(fi*n)
        res[i,0] = r*np.sin(fi)
        res[i,1] = r*np.cos(fi)

    res[-1] = res[0]
    return res

def RotateMeshByAngle(mesh, angle):
    nodes, triangles, segment_indices, trig_neighbors, node_neighbours, triangle_indeces = mesh

    theta = np.deg2rad(angle)
    rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    for (n_node, node) in enumerate(nodes):
        x,y = node
        
        new_node = np.dot(rot, node)
        nodes[n_node] = new_node

        if n_node % 5000 == 0:
            print(f"rotated {n_node} nodes")

    return nodes, triangles, segment_indices, trig_neighbors, node_neighbours, triangle_indeces


def main():
    import matplotlib.pyplot as plt
    
    points = GenerateCurvyMesh(600, 5, 1, 0.3)

    x,y = points[:,0], points[:,1]

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    
    plt.scatter(x, y)
    plt.show()


if __name__ == "__main__":
    main()