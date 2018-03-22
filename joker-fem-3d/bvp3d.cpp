#include "bvp3d.h"
#include "integrands3d.h"
#include "linear_system3d.h"

void SolveBVP (const ProblemData& data, const Parameters& param,
    std::vector<FunctionP1>& sol)
{
    int nodes_num = data.mesh.nodes_num;
    std::vector<BasisFunction> phi(nodes_num, data.mesh);
    for (int i = 0; i < nodes_num; ++i)
        phi[i] = BasisFunction(data.mesh, i);

    LinearSystem sys(data.N, nodes_num, param);
    for (int i = 0; i < data.N; ++i) {
        for (int j = 0; j < data.mesh.nodes_num; ++j) {
            for (auto s : data.mesh.adjacent_nodes[j]) {
                sys.AddCoeff(i, j, i, s,
                    data.a[i] * Integrate(grad(phi[s]) * grad(phi[j]))
                        + BoundaryIntegrate(data.b[i] * phi[s] * phi[j])
                );
            }
            for (int k = 0; k < data.N; ++k) {
                for (auto s : data.mesh.adjacent_nodes[j]) {
                    sys.AddCoeff(i, j, k, s,
                        Integrate(data.q[i][k] * phi[s] * phi[j])
                    );
                }
            }
            sys.AddRhs(i, j,
                Integrate(data.g[i] * phi[j])
                    + BoundaryIntegrate(data.b[i] * data.w[i] * phi[j])
            );
        }
    }

    sys.Solve(sol);
}
