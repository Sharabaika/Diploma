class MeshConfig:
    MESH_FILE_HEADER = "mesh file\n"

    MESH_FILE_NODES_HEADER = "nodes\n"
    MESH_FILE_NODES_FORMAT = "{0} {1} \n"

    MESH_FILE_TRIANGLES_HEADER = "triangles\n"
    MESH_FILE_TRIANGLES_FORMAT = "{0} {1} {2}\n"

    MESH_FILE_SEGMENT_TAGS_HEADER = "segment tags\n"
    MESH_FILE_SEGMENT_TAGS_FORMAT = "{0}\n"

    MESH_FILE_TRIANGLE_NEIGBOUR_HEADER = "trig neighbours\n"
    MESH_FILE_TRIANGLE_NEIGBOUR_FORMAT = "{0}\n"

    MESH_FILE_NODE_NEIGBOUR_HEADER = "node neighbours\n"
    MESH_FILE_NODE_NEIGBOUR_FORMAT = "{0}\n"

    MESH_FILE_TRIANGLE_INDICES_HEADER = "trianlge indices\n"
    MESH_FILE_TRIANGLE_INDICES_FORMAT = "{0}\n"

def WriteMeshToDAT(filename, mesh):
    new_file=open("{filename}.dat".format(filename=filename),mode="w",encoding="utf-8")

    
    new_file.write("{0}\n".format(len(mesh)))
    for v in mesh:
        line = "{0} {1} {2}\n".format(*v)
        new_file.write(line)

    new_file.close()

def SaveMesh(path, filename, nodes, triangles, segment_tags, trig_neighbors, node_neighbours, trianlge_indices):
    import numpy as np
    import os
    
    if not os.path.exists(path):
        os.makedirs(path)

    name = "{filename}.dat".format(filename=filename)
    with open(os.path.join(path, name), 'w') as new_file:
        new_file.write(MeshConfig.MESH_FILE_HEADER)

        new_file.write(MeshConfig.MESH_FILE_NODES_HEADER)
        new_file.write(str(len(nodes)) + '\n')
        new_file.writelines(map(lambda node:MeshConfig.MESH_FILE_NODES_FORMAT.format(*node), nodes))

        new_file.write(MeshConfig.MESH_FILE_TRIANGLES_HEADER)
        new_file.write(str(len(triangles)) + '\n')
        new_file.writelines(map(lambda trig:MeshConfig.MESH_FILE_TRIANGLES_FORMAT.format(*trig), triangles))

        new_file.write(MeshConfig.MESH_FILE_SEGMENT_TAGS_HEADER)
        new_file.write(str(len(segment_tags)) + '\n')
        new_file.writelines(map(lambda tags:MeshConfig.MESH_FILE_SEGMENT_TAGS_FORMAT.format(tags), segment_tags))

        new_file.write(MeshConfig.MESH_FILE_TRIANGLE_NEIGBOUR_HEADER)
        new_file.write(str(len(trig_neighbors)) + '\n')
        new_file.writelines(map(lambda trigs:MeshConfig.MESH_FILE_TRIANGLE_NEIGBOUR_FORMAT.format(trigs), trig_neighbors))

        new_file.write(MeshConfig.MESH_FILE_NODE_NEIGBOUR_HEADER)
        new_file.write(str(len(node_neighbours)) + '\n')
        new_file.writelines(map(lambda nodes:MeshConfig.MESH_FILE_NODE_NEIGBOUR_FORMAT.format(nodes), node_neighbours))

        new_file.write(MeshConfig.MESH_FILE_TRIANGLE_INDICES_HEADER)
        new_file.write(str(len(trianlge_indices)) + '\n')
        new_file.writelines(map(lambda index:MeshConfig.MESH_FILE_TRIANGLE_INDICES_FORMAT.format(index), trianlge_indices))


def main():
    from MeshGenerator import GenerateCurvyMesh

    N = 250
    n = 3
    R = 1
    dr = 0.3

    points = GenerateCurvyMesh(N, n, R, dr)
    WriteMeshToDAT(f"SavedCurves/N{N}_n{n}_R{R}_dr{dr}", points)
    # for n in [3,4]:
    #     points = GenerateCurvyMesh(N, n, R, dr)
    #     WriteMeshToDAT(f"SavedCurves/N{N}_n{n}_R{R}_dr{dr}", points)

if __name__ == "__main__":
    main()