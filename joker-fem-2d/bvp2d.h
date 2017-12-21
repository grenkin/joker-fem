#ifndef BVP2D_H_INCLUDED
#define BVP2D_H_INCLUDED

#include "splines2d.h"

struct ProblemData {
    int N;  // number of equations
    Mesh& mesh;
    std::vector<double> a;  // diffusion coefficients
    std::vector<BoundaryFunctionP0> b;  // coefficient in boundary conditions
    std::vector<BoundaryFunctionP1> w;  // function in boundary conditions
    std::vector<std::vector<FunctionP1> > q;  // reaction terms functions
    std::vector<FunctionP1> g;  // right-hand sides

    ProblemData (int _N, Mesh& _mesh)
        : N(_N), mesh(_mesh)
    {
        a.resize(N);
        b.resize(N, mesh);
        w.resize(N, mesh);
        q.resize(N);
        for (int i = 0; i < N; ++i)
            q[i].resize(N, mesh);
        g.resize(N, mesh);
    }
};

#endif // BVP2D_H_INCLUDED
