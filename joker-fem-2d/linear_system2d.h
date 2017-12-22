#ifndef LINEAR_SYSTEM2D_H_INCLUDED
#define LINEAR_SYSTEM2D_H_INCLUDED

#include "bvp2d.h"
#include <boost/numeric/mtl/mtl.hpp>

class LinearSystem {
    typedef mtl::mat::inserter<mtl::compressed2D<double>,
        mtl::update_plus<double> > inserter_type;

    int N, nodes_num;
    int unknowns;
    mtl::compressed2D<double> A;
    mtl::dense_vector<double> b;
    inserter_type* ins;
    const Parameters& param;

    // calculate the index of an unknown in the linear system
    int get_index (int i, int j)
    {
        return i * nodes_num + j;
    }
public:
    LinearSystem (int _N, int _nodes_num, const Parameters& _param)
        : N(_N), nodes_num(_nodes_num), unknowns(N * nodes_num),
            A(unknowns, unknowns), b(unknowns), param(_param)
    {
        int avg_vertex_degree = 6;
        int elements_per_row = N * (avg_vertex_degree + 1);
        ins = new inserter_type(A, elements_per_row);

        for (int i = 0; i < unknowns; ++i)
            b[i] = 0.0;
    }

    // add a coefficient to the linear system
    // i - differential equation, j - node
    void AddCoeff (int eq_i, int eq_j, int var_i, int var_j, double val)
    {
        int eq = get_index(eq_i, eq_j);
        int var = get_index(var_i, var_j);
        (*ins)[eq][var] << val;
    }

    // add a value to the right-hand side of the linear system
    // i - differential equation, j - node
    void AddRhs (int eq_i, int eq_j, double val)
    {
        int eq = get_index(eq_i, eq_j);
        b[eq] += val;
    }

    // sol contains the initial guess and
    // will contain the solution after function call
    void Solve (std::vector<FunctionP1>& sol);
};


#endif // LINEAR_SYSTEM2D_H_INCLUDED
