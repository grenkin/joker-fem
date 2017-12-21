#ifndef LINEAR_SYSTEM2D_H_INCLUDED
#define LINEAR_SYSTEM2D_H_INCLUDED

class LinearSystem {
    int N, nodes_num;
public:
    LinearSystem (int _N, int _nodes_num)
        : N(_N), nodes_num(_nodes_num)
    {

    }

    // add a coefficient to the linear system
    // i - differential equation, j - node
    void AddCoeff (int eq_i, int eq_j, int var_i, int var_j, double val)
    {

    }

    // add a value to the right-hand side of the linear system
    // i - differential equation, j - node
    void AddRhs (int eq_i, int eq_j, double val)
    {

    }
};


#endif // LINEAR_SYSTEM2D_H_INCLUDED
