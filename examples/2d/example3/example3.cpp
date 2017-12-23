/*
-\Delta u = 2 * sin(x) * sin(y)
(0, pi) \times (0, pi)
-u_x + u|_{x=0} = - sin(y)
u_x + u|_{x=pi} = - sin(y)
-u_y + u|_{y=0} = - sin(x)
u_y + u|_{y=pi} = - sin(x)
Exact solution: u(x, y) = sin(x) * sin(y)
*/

#include <joker-fem-2d/bvp2d.h>
#include <iostream>
#include <cmath>

using namespace std;

double gfun(double x, double y)
{
    return 2 * sin(x) * sin(y);
}

double exact(double x, double y)
{
    return sin(x) * sin(y);
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

    double L = 3.14159265358979323846;
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
            if (eq(x, 0) || eq(x, L))
                data.w[0].values[i] = - sin(y);
            else if (eq(y, 0) || eq(y, L))
                data.w[0].values[i] = - sin(x);
            else
                throw;
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
