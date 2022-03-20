def ReadRaw(mesh, default_tags, *segment_files):
    import numpy as np

    points_len = int(mesh.readline())
    points = np.empty((points_len, 2))
    region_tags = [default_tags for _ in range(points_len)]

    for i in range(points_len):
        x, y, z = np.asfarray(mesh.readline().split(), dtype = float)
        points[i] = x,y

    for  bordet_segment, segment_tags in segment_files:
        is_complete = False
        border_len = int(bordet_segment.readline())

        while not is_complete:
            for _ in range(border_len):
                line = np.asfarray(bordet_segment.readline().split(), dtype = float)
                if len(line) == 3:
                    x, y, z = line
                    for point_index in [i for i in range(points_len) if points[i][0]==x and points[i][1]==y]:
                        region_tags[point_index] = segment_tags
            
            last_line = bordet_segment.readline()
            if last_line:
                border_len = int(last_line)
            else:
                is_complete = True
            

    triangles = []
    triangles_len = int(mesh.readline())
    for i in range(triangles_len):
        a, b, c, _ = mesh.readline().split()
        triangles.append([int(a)-1, int(b)-1, int(c)-1])

    trig_neighbours = [[] for _ in range(points_len)]
    node_neighbours = [[] for _ in range(points_len)]
    for point_index in range(points_len):
        for n_triangle, triangle in enumerate(triangles):
            if point_index in triangle:
                trig_neighbours[point_index].append(n_triangle)

        for n_trig_neighbour in trig_neighbours[point_index]:
            for p in triangles[n_trig_neighbour]:
                node_neighbours[point_index].append(p)
        node_neighbours[point_index] = (list(dict.fromkeys(node_neighbours[point_index])))
        if point_index in node_neighbours[point_index]:
            node_neighbours[point_index].remove(point_index)


    def most_common(lst):
        return max(set(lst), key=lst.count)
    
    N_triangles = len(triangles)
    triangle_indecies = [-1 for _ in range(N_triangles)]
    for n_triangle in range(N_triangles):
        indices = [region_tags[i] for i in triangles[n_triangle]]

        bConductor = any([index == 0 for index in indices])
        if bConductor:
            triangle_indecies[n_triangle] = 0
            continue

        bMedium = any([index == 2 for index in indices])
        if bMedium:
            triangle_indecies[n_triangle] = 2 
            continue

        bVoid = any([index == 4 for index in indices])
        if bVoid:
            triangle_indecies[n_triangle] = 4
            continue          

    return points, triangles, region_tags, trig_neighbours, node_neighbours, triangle_indecies


def ReadSaved(filename):
    import numpy as np
    import os
    file = open(filename, "r")

    def ReadArrays(dim):
        array_header = file.readline()
        array_len =  int(file.readline())
        array = np.empty((array_len, dim))
        for i in range(array_len):
            array[i] = np.asfarray(file.readline().split(), dtype = float)
        return array

    def ReadLists(data_type = float):
        list_header = file.readline()
        list_len =  int(file.readline())
        res_list = []
        for i in range(list_len):
            data = list(data_type(string) for string in file.readline().replace(',','').strip("[]\n").split(" "))
            res_list.append(data)
        return res_list
    
    def ReadInts():
        header = file.readline()
        list_len = int(file.readline())
        res_list = [-1 for _ in range(list_len)]
        
        for i in range(list_len):
            res_list[i] = int(file.readline())

        return res_list


    file_header = file.readline()

    nodes = ReadArrays(2)
    triangles = ReadLists(int)
    tags = ReadInts()
    trig_neighbours = ReadLists(int)
    node_neighbours = ReadLists(int)
    trianlge_indices = ReadInts()

    return nodes, triangles, tags, trig_neighbours, node_neighbours, trianlge_indices
