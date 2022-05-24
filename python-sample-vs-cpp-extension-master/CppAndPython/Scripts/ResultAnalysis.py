from cmath import nan
from math import atan2
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from os import path
from Scripts.MeshReader import ReadSaved
from Scripts.settings import MeshNames


class ResultAnalysis:
    def __init__(self, folder, result_name):
        self.path = os.path.join(folder, result_name)
        
        self.logs = None
        self.saved = {}
        self.params = None

        self.loaded_mesh = None

    def GetLogs(self):
        if self.logs is None:
            logs_path = os.path.join(self.path, "logs.csv")
            self.logs = pd.read_csv(logs_path, index_col=0)
        return self.logs

    def GetSavedResults(self, category="nodes"):
        if not category in self.saved:
            saved_path = os.path.join(self.path, f"{category}.csv")
            self.saved[category] = pd.read_csv(saved_path)
        return self.saved[category]

    def GetParams(self):
        if self.params is None:
            params_path = os.path.join(self.path, "params.csv")
            self.params = pd.read_csv(params_path)
        return self.params

    def GetParam(self, name):
        return self.GetParams().iloc[0][name]

    def GetMeshName(self):
        return self.GetParams().iloc[0]['mesh_name']

    def GetMesh(self):
        if self.loaded_mesh is None:
            mesh_name = self.GetMeshName()
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

class MagneticsAnalysis(ResultAnalysis):
    def __init__(self, folder, result_name):
        super().__init__(folder, result_name)

    def MakeExplicit(H_triangles, H_nodes, Fi, Mu, nodes, triangles, segment_indices, trig_neighbors, node_neighbours, triangle_indeces):
        res = MagneticsAnalysis("", "")

        res.saved["nodes"] = {
            "H_nodes" : H_nodes,
            "Fi" : Fi
        }

        res.saved["triangles"] = {
            "H" : H_triangles,
            "Mu" : Mu
        }

        res.loaded_mesh = (nodes, triangles, segment_indices, trig_neighbors, node_neighbours, triangle_indeces)
        
        return res

    def GetFi(self):
        return np.array(self.GetSavedResults("nodes")["Fi"])

    def GetH(self):
        return np.array(self.GetSavedResults("triangles")["H"])
    
    def GetH_Nodes(self):
        return np.array(self.GetSavedResults("nodes")["H_nodes"])

    def GetMu(self):
        return np.array(self.GetSavedResults("triangles")["Mu"])


class DynamycsAnalysis(ResultAnalysis):
    def __init__(self, folder, result_name):
        super().__init__(folder, result_name)

    def MakeExplicit(Psi, W, T, nodes, triangles, segment_indices, trig_neighbors, node_neighbours, triangle_indeces):
        res = DynamycsAnalysis("", "")
        res.saved["nodes"] = {
            "Psi" : Psi,
            "W" : W,
            "T" : T
        }

        res.loaded_mesh = (nodes, triangles, segment_indices, trig_neighbors, node_neighbours, triangle_indeces)
        
        return res

    def GetPsi(self):
        return self.GetSavedResults()["Psi"]

    def GetW(self):
        return self.GetSavedResults()["W"]

    def GetT(self):
        return self.GetSavedResults()["T"]

    def GetError(self, variable):
        return self.GetLogs()[variable]

    def CalculateLocalNulselt(self):
        nodes, triangles, segment_indices, trig_neighbors, node_neighbours, trianlge_indices = self.GetMesh()
        T = self.GetT()
        
        OUTER_BORDER_INDEX = 11   # 2
        def is_wall(n):
            return OUTER_BORDER_INDEX in segment_indices[n]

        fis = []
        nuls = []

        n_triangles = len(triangles)
        for n_triangle in range(n_triangles):
            triangle = triangles[n_triangle]
            
            border_nodes = []
            for n_node in triangle:
                if is_wall(n_node):
                    border_nodes.append(n_node)
            
            if len(border_nodes) == 2:
                a,b,c = triangle
                if is_wall(a) and is_wall(b):
                    pass
                elif is_wall(b) and is_wall(c):
                    a,b,c = b,c,a
                elif is_wall(c) and is_wall(a):
                    a,b,c = c,a,b
                
                ax,ay = nodes[a]
                bx,by = nodes[b]
                cx,cy = nodes[c]

                ab = np.sqrt((ax-bx)**2+(ay-by)**2)
                bc = np.sqrt((bx-cx)**2+(by-cy)**2)
                ca = np.sqrt((cx-ax)**2+(cy-ay)**2)

                p = (ab+bc+ca)*0.5

                h = 2.0*np.sqrt(p*(p-ab)*(p-bc)*(p-ca))/ab

                dtdy = T[c]/h
                nuls.append(dtdy)
                
                midx, midy = (ax+bx)*0.5, (ay+by)*0.5
                fi = np.rad2deg(atan2(midy, midx))
                fis.append(fi)

        order = np.argsort(fis)
        return np.array(fis)[order], np.array(nuls)[order]

        

    def CalculateNulselt(self):
        nodes, triangles, segment_indices, trig_neighbors, node_neighbours, trianlge_indices = self.GetMesh()
        T = self.GetT()

        MEDIUM_REGION_INDEX = 2        
        OUTER_BORDER_INDEX = 3
        def is_wall(n):
            return OUTER_BORDER_INDEX == segment_indices[n]


        integral = 0.0
        border_length = 0.0
        n_triangles = len(triangles)
        for n_triangle in range(n_triangles):
            triangle = triangles[n_triangle]
            
            border_nodes = []
            for n_node in triangle:
                if is_wall(n_node):
                    border_nodes.append(n_node)
            
            if len(border_nodes) == 2:
                a,b,c = triangle
                if is_wall(a) and is_wall(b):
                    pass
                elif is_wall(b) and is_wall(c):
                    a,b,c = b,c,a
                elif is_wall(c) and is_wall(a):
                    a,b,c = c,a,b

                if segment_indices[c] != MEDIUM_REGION_INDEX:
                    continue

                ax,ay = nodes[a]
                bx,by = nodes[b]
                cx,cy = nodes[c]

                ab = np.sqrt((ax-bx)**2+(ay-by)**2)
                bc = np.sqrt((bx-cx)**2+(by-cy)**2)
                ca = np.sqrt((cx-ax)**2+(cy-ay)**2)

                p = (ab+bc+ca)*0.5

                h = 2.0*np.sqrt(p*(p-ab)*(p-bc)*(p-ca))/ab

                dtdy = T[c]/h
                integral += dtdy*ab
                border_length += ab

        return integral
                

class NuseltTable:
    path_to_table = "Nus.csv"

    result_name_literal = "result_name"
    mesh_name_literal = "mesh_name"
    mesh_name_short_literal = "mesh_name_short"
    ra_literal = "Ra"
    ram_literal = "Ram"
    h_literal = "H"
    nuselt_result_literal = "nu"

    columns=[
            result_name_literal,
            mesh_name_literal,
            mesh_name_short_literal,
            ra_literal,
            ram_literal, 
            h_literal,
            nuselt_result_literal
        ]
    # index | result_name | mesh_name |nu
    # 0     | "aboa"      | "n0"      |40

    # index | result_name | mesh_name                            | mesh_name_short | Ra | Ram | H | nu
    # 0     | "aboa"      | "Computational/n0_N100-500-500-100"  | "n0"            | 0  | 100 | 5 | 40


    def __init__(self, table, path = path_to_table) -> None:
        self.table = table
        self.path = path
    
    def LoadFromCSV(table_path = path_to_table):
        table = pd.read_csv(table_path)
        return NuseltTable(table, table_path)

    def SaveToCSV(self):
        self.table.to_csv(self.path, index=False)

    def MakeRow(result_name, nu):            
            analysis = DynamycsAnalysis("SavedResults", result_name) 
            magnetics = MagneticsAnalysis("SavedMagnetics", analysis.GetParam("magnetics_result_name"))

            mesh_name = analysis.GetMeshName()
            short_name = MeshNames.GetShortName(mesh_name)
            ra = analysis.GetParam("Ra")
            ram = analysis.GetParam("Ram")
            H = magnetics.GetParam("H0")
            
            return [result_name, mesh_name, short_name, ra, ram, H, nu]

    def GetNuselt(self, result_name, b_calculate_if_missing = True, b_resafe = True):
        if NuseltTable.result_name_literal in self.table.columns.values.tolist():
            results = self.table.loc[self.table[NuseltTable.result_name_literal] == result_name]
            lres = len(results)
            if lres > 0:
                if lres > 1:
                    print(f"more than 1 result for {result_name}")
                return results.iloc[0][NuseltTable.nuselt_result_literal]
        if b_calculate_if_missing and path.exists(f"SavedResults/{result_name}/nodes.csv"):
            analysis = DynamycsAnalysis("SavedResults", result_name) 
            nuselt = analysis.CalculateNulselt()
            df = pd.DataFrame([NuseltTable.MakeRow(result_name, nuselt)], columns=NuseltTable.columns)
            print(f"CALCULATING NUSELT = {nuselt} FOR {result_name}")

            self.table = pd.concat([self.table, df])

            if b_resafe:
                self.SaveToCSV()

            return nuselt

        return np.NaN

    def RedoTable(self):
        result_lst = []
        for index, row in self.table.iterrows():
            result_name = row[NuseltTable.result_name_literal]
            nu = row[NuseltTable.nuselt_result_literal]

            result_lst.append(NuseltTable.MakeRow(result_name, nu))
        
        result_df = pd.DataFrame(result_lst, columns=NuseltTable.columns)
        print(result_df)

        self.table = result_df
        self.SaveToCSV()


class NuseltFitTable:
    path_to_table = "Nus_fit.csv"
    mesh_name_literal = "mesh_name"
    mult_literal = "mult"
    pow_literal = "pow"

    def __init__(self, table, path) -> None:
        self.table = table
        self.path = path
    
    def LoadFromCSV(table_path = path_to_table):
        table = pd.read_csv(table_path)
        return NuseltTable(table, table_path)

    def SaveToCSV(self):
        self.table.to_csv(self.path, index=False)

    

        
            


