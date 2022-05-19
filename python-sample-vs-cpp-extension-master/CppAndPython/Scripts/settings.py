
class ParamsSettings:
    ram_range = [1000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000]
    ram_range_short = [1000, 20000, 40000, 60000, 80000, 100000]

class MeshNames:
    n0 = "Computational/n0_N100-500-500-100"
    n0_N1000 = "Computational/n0_N1000"

    n0_250 = "Computational/n0_N50-250-250-50"
    n0_375 = "Computational/n0_N75-375-375-75"
    n0_500 = "Computational/n0_N100-500-500-100"
    n0_600 = "Computational/n0_N150-600-600-150"
    n0_700 = "Computational/n0_N100-700-700-100"

    n2_600_dr_03 = "Computational/n2_N100-600-600-100_dr_03"
    n2_600_dr_03_rot = "Computational/n2_N100-600-600-100_dr_03_rot"

    n3_600_dr_03 = "Computational/n3_N100-600-600-100_dr_03"
#old
    n_2_dr_03 = "Computational/n_2_dr_0.3"
    n_2_dr_03_r = "Computational/n_2_dr_0.3_r"

    n_3_dr_03 = "Computational/n_3_dr_0.3"
    n_3_dr_03_N_500 = "Computational/n_3_dr_0.3_N_500"
    n_3_dr_03_r = "Computational/n_3_dr_0.3_r"

    mesh_list = [n0_600, n2_600_dr_03, n2_600_dr_03_rot]
    mesh_list_n0 = [n0_250, n0_375, n0_500, n0_600, n0_700]

    def GetShortName(mesh):
        prefix = "Computational/"
        if mesh.startswith(prefix):
            return mesh[len(prefix):]
        return mesh

class ResultName:
    def MakeName(mesh_name, ram):
        return f"{mesh_name}/Ram_{ram}"

    def ExtendName(name):
        return f"SavedResults/{name}"

    def MakeExtendedName(mesh_name, ram):
        return ResultName.ExtendName(ResultName.MakeName(mesh_name, ram))


class MagneticsResultName:
    def MakeName(mesh):
        return f"{mesh}/magnetics_H_5_chi0_2_mu_1000"