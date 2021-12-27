#pragma once
#include <vector>

using std::vector;

typedef vector<vector<int>> JaggedArr;
typedef vector<double> Arr;

#define MeshData Arr x, Arr y, JaggedArr triangles, JaggedArr segment_indices, JaggedArr trig_neighbors, JaggedArr node_neighbours
#define FluidDynamicsParams double Re, double Vx
#define FluidDynamicsNumericParams double QPsi, double QW, double max_error, int max_cycles
#define FluidDynamicsOUTResults Arr& Psi, Arr& W
#define FluidDynamicsOutErrors Arr& Delta_Psi, Arr& Delta_W

void SolveFluid_Implementation(MeshData, FluidDynamicsParams, FluidDynamicsNumericParams, FluidDynamicsOUTResults, FluidDynamicsOutErrors);