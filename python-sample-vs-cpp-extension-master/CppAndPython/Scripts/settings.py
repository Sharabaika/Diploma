
class ParamsSettings:
    ram_range = [1000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000]

class MeshNames:
    n0 = "Computational/n0_N100-500-500-100"
    n_2_dr_03 = "Computational/n_2_dr_0.3"
    n_2_dr_03_r = "Computational/n_2_dr_0.3_r"
    n_3_dr_03 = "Computational/n_3_dr_0.3"

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