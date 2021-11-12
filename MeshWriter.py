import numpy as np

def WriteMesh(filename, mesh):
    new_file=open("{filename}.dat".format(filename=filename),mode="w",encoding="utf-8")

    
    new_file.write("{0}\n".format(len(mesh)))
    for v in mesh:
        line = "{0} {1} {2}\n".format(*v)
        new_file.write(line)

    new_file.close()


def main():
    import MeshGenerator as gen
    
    points = gen.GenerateCurvyMesh(100, 5, 0.3, 0.1)

    WriteMesh("testfile", points)

if __name__ == "__main__":
    main()