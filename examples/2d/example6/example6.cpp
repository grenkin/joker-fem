/*
-\Delta u = xy
Cirle: {(x,y): x^2 + y^2 < 4}
u_n + u = - x * y / 3
Exact solution: u(x, y) = 1/12 * (4 - r^2) * x * y
http://sci.alnam.ru/book_eqf.php?id=8
*/

#include <joker-fem-2d/bvp2d.h>
#include <iostream>
#include <cmath>

using namespace std;

double gfun(double x, double y)
{
    return x * y;
}

double exact(double x, double y)
{
    double r = sqrt(x * x + y * y);
    return 1./12. * (4 - r * r) * x * y;
}

int main() {
    int tests_num = 4;
    int n_values[tests_num] = {20, 40, 80, 160};
    string files_names[tests_num] =
        {"tr20.mesh", "tr40.mesh", "tr80.mesh", "tr160.mesh"};

    int N = 1;

    double rms_old;
    for (int test = 0; test < tests_num; ++test) {
        Mesh mesh(files_names[test]);
        ProblemData data(N, mesh);
        data.a[0] = 1.0;
        for (int i = 0; i < mesh.boundary_edges_num; ++i)
            data.b[0].values[i] = 1.0;
        for (int i = 0; i < mesh.boundary_nodes_num; ++i) {
            double x = mesh.nodes[mesh.boundary_nodes[i]].x;
            double y = mesh.nodes[mesh.boundary_nodes[i]].y;
            data.w[0].values[i] = - x * y / 3.;
        }
        for (int i = 0; i < mesh.nodes_num; ++i) {
            data.g[0].values[i] = gfun(mesh.nodes[i].x, mesh.nodes[i].y);
            data.q[0][0].values[i] = 0.0;
        }

        vector<FunctionP1> sol(N, mesh);
        for (int i = 0; i < mesh.nodes_num; ++i)
            sol[0].values[i] = 0.0;
        SolveBVP(data, Parameters(), sol);

        double rms = 0.0;
        for (int i = 0; i < mesh.nodes_num; ++i) {
            rms += pow(sol[0].values[i] - exact(mesh.nodes[i].x, mesh.nodes[i].y), 2);
        }
        rms = sqrt(rms / mesh.nodes_num);
        cout << "n = " << n_values[test] << "   rms = " << rms << endl;
        if (test > 0)
            cout << "rms_old / rms = " << rms_old / rms << endl;
        cout << endl;
        rms_old = rms;
    }

    return 0;
}
