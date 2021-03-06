#ifndef BVP3D_H_INCLUDED
#define BVP3D_H_INCLUDED

#include "splines3d.h"

struct ProblemData {
    int N;  // number of equations
    Mesh& mesh;
    std::vector<double> a;  // diffusion coefficients
    std::vector<BoundaryFunctionP0> b;  // coefficients in boundary conditions
    std::vector<BoundaryFunctionP1> w;  // functions in boundary conditions
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

struct Parameters {
    double linear_sys_tol;
    int max_linear_sys_iterations;

    Parameters ()
    {
        linear_sys_tol = 1e-7;  // relative error in the linear system solver
        // maximal number of iterations in the linear system solver
        max_linear_sys_iterations = 100000;
    }
};

// sol contains the initial guess when calling the function
// and contains the solution after function call
void SolveBVP (const ProblemData&, const Parameters&, std::vector<FunctionP1>& sol);

#endif // BVP3D_H_INCLUDED
