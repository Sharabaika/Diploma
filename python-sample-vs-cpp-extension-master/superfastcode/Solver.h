#pragma once
#include <vector>

using std::vector;

typedef vector<double> Arr;
typedef vector<int> ArrInt;
typedef vector<ArrInt> JaggedArr;

const static int SOLID_BORDER_INDEX = 1000;
const static int MEDIUM_INDEX = 2000;

const double PI = 3.141592653589793;

void SolveFluid_Implementation(const Arr x, const Arr y, const JaggedArr triangles, const JaggedArr segment_indices, const JaggedArr trig_neighbors, const JaggedArr node_neighbours,
	const double Re, const double Vx,
	const double QPsi, const double QW, const double max_error, const int max_cycles,
	Arr& Psi, Arr& W,
	Arr& Delta_Psi, Arr& Delta_W);