#include "linear_system2d.h"
#include <boost/numeric/itl/itl.hpp>

void LinearSystem::Solve (std::vector<FunctionP1>& sol)
{
    if (!ins)
        throw "The system was already solved.";
    else {
        delete ins;
        ins = 0;

        // Set the initial guess
        mtl::dense_vector<double> x(unknowns);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < nodes_num; ++j) {
                int eq = get_index(i, j);
                x[eq] = sol[i].values[j];
            }
        }

        // Solve the linear system
        itl::pc::ilu_0<mtl::compressed2D<double> > P(A);
        itl::basic_iteration<double> iter(b, param.max_linear_sys_iterations,
            param.linear_sys_tol);
        bicgstab(A, x, b, P, iter);

        // Set the solution
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < nodes_num; ++j) {
                int eq = get_index(i, j);
                sol[i].values[j] = x[eq];
            }
        }
    }
}
