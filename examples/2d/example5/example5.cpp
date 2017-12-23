/*
A nonlinear problem in 2D
-\Delta u + u^2 = -4 + (x^2 + y^2)^2
(0, 1) \times (0, 1)
-u_x + u|_{x=0} = y^2
u_x + u|_{x=1} = 3 + y^2
-u_y + u|_{y=0} = x^2
u_y + u|_{y=1} = 3 + x^2
Exact solution: u(x, y) = x^2 + y^2
*/

#include <joker-fem-2d/bvp2d.h>
#include <iostream>
#include <cmath>

using namespace std;

double gfun(double x, double y)
{
    return -4 + pow(x * x + y * y, 2);
}

double exact(double x, double y)
{
    return x * x + y * y;
}

bool eq (double a, double b)
{
    return fabs(a - b) < 1e-5;
}

int main() {
    int tests_num = 5;
    int n_values[tests_num] = {20, 40, 80, 160, 320};
    string files_names[tests_num] =
        {"tr20.mesh", "tr40.mesh", "tr80.mesh", "tr160.mesh", "tr320.mesh"};

    double L = 1.0;
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
            if (eq(x, 0) && eq(y, 0))
                data.w[0].values[i] = 0.5 * y * y + 0.5 * x * x;
            else if (eq(x, 0) && eq(y, L))
                data.w[0].values[i] = 0.5 * y * y + 0.5 * (3 + x * x);
            else if (eq(x, L) && eq(y, 0))
                data.w[0].values[i] = 0.5 * (3 + y * y) + 0.5 * x * x;
            else if (eq(x, L) && eq(y, L))
                data.w[0].values[i] = 0.5 * (3 + y * y) + 0.5 * (3 + x * x);
            else if (eq(x, 0))
                data.w[0].values[i] = y * y;
            else if (eq(x, L))
                data.w[0].values[i] = 3 + y * y;
            else if (eq(y, 0))
                data.w[0].values[i] = x * x;
            else if (eq(y, L))
                data.w[0].values[i] = 3 + x * x;
            else
                throw;
        }
        vector<FunctionP1> sol(N, mesh);
        for (int i = 0; i < mesh.nodes_num; ++i)
            sol[0].values[i] = 0.0;
        // set the initial guess
        FunctionP1 sol_old(mesh);
        for (int i = 0; i < mesh.nodes_num; ++i)
            sol_old.values[i] = 0.0;

        int iterations = 0;
        while (1) {
            cout << ++iterations << " ";
            // linearized equation:
            // -\Delta u + 2 * u0 * u = -4 + (x^2 + y^2)^2 + u0^2
            for (int i = 0; i < mesh.nodes_num; ++i) {
                data.g[0].values[i] = gfun(mesh.nodes[i].x, mesh.nodes[i].y)
                    + pow(sol_old.values[i], 2);
                data.q[0][0].values[i] = 2 * sol_old.values[i];
            }

            SolveBVP(data, Parameters(), sol);

            double max_diff = 0.0;
            for (int i = 0; i < mesh.nodes_num; ++i)
                max_diff = fmax(max_diff, fabs(sol[0].values[i] - sol_old.values[i]));
            if (max_diff < 1e-5)
                break;

            for (int i = 0; i < mesh.nodes_num; ++i)
                sol_old.values[i] = sol[0].values[i];
        }

        double rms = 0.0;
        for (int i = 0; i < mesh.nodes_num; ++i) {
            rms += pow(sol[0].values[i] - exact(mesh.nodes[i].x, mesh.nodes[i].y), 2);
        }
        rms = sqrt(rms / mesh.nodes_num);
        cout << "\nn = " << n_values[test] << "   rms = " << rms << endl;
        if (test > 0)
            cout << "rms_old / rms = " << rms_old / rms << endl;
        cout << endl;
        rms_old = rms;
    }

    return 0;
}
