/*
-\Delta u + u^2 = -6 + (x^2 + y^2 + z^2)^2
-u_x + u|_{x=-1} = 3 + y^2 + z^2
u_x + u|_{x=1} = 3 + y^2 + z^2
-u_y + u|_{y=-1} = 3 + x^2 + z^2
u_y + u|_{y=1} = 3 + x^2 + z^2
-u_z + u|_{z=-1} = 3 + x^2 + y^2
u_z + u|_{z=1} = 3 + x^2 + y^2
Exact solution: u(x, y, z) = x^2 + y^2 + z^2
*/

#include <joker-fem-3d/bvp3d.h>
#include <iostream>
#include <cmath>

using namespace std;

double gfun(double x, double y, double z)
{
    return -6 + pow(x * x + y * y + z * z, 2);
}

double exact(double x, double y, double z)
{
    return x * x + y * y + z * z;
}

bool eq (double a, double b)
{
    return fabs(a - b) < 1e-5;
}

int main()
{
    int tests_num = 3;
    int n_values[tests_num] = {10, 20, 40};
    string files_names[tests_num] =
        {"tr10_cube.mesh", "tr20_cube.mesh", "tr40_cube.mesh"};

    int N = 1;
    double rms_old;
    for (int test = 0; test < tests_num; ++test) {
        Mesh mesh(files_names[test]);
        ProblemData data(N, mesh);
        data.a[0] = 1.0;
        for (int i = 0; i < mesh.boundary_triangles_num; ++i)
            data.b[0].values[i] = 1.0;
        for (int i = 0; i < mesh.boundary_nodes_num; ++i) {
            double x = mesh.nodes[mesh.boundary_nodes[i]].x;
            double y = mesh.nodes[mesh.boundary_nodes[i]].y;
            double z = mesh.nodes[mesh.boundary_nodes[i]].z;
            if (eq(x, -1) || eq(x, 1))
                data.w[0].values[i] = 3 + y * y + z * z;
            else if (eq(y, -1) || eq(y, 1))
                data.w[0].values[i] = 3 + x * x + z * z;
            else if (eq(z, -1) || eq(z, 1))
                data.w[0].values[i] = 3 + x * x + y * y;
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
            // -\Delta u + 2 * u0 * u = -6 + (x^2 + y^2 + z^2)^2 + u0^2
            for (int i = 0; i < mesh.nodes_num; ++i) {
                data.g[0].values[i] =
                    gfun(mesh.nodes[i].x, mesh.nodes[i].y, mesh.nodes[i].z)
                    + pow(sol_old.values[i], 2);
                data.q[0][0].values[i] = 2 * sol_old.values[i];
            }

            SolveBVP(data, Parameters(), sol);

            double max_diff = 0.0;
            for (int i = 0; i < mesh.nodes_num; ++i)
                max_diff = fmax(max_diff,
                    fabs(sol[0].values[i] - sol_old.values[i]));
            if (max_diff < 1e-5)
                break;

            for (int i = 0; i < mesh.nodes_num; ++i)
                sol_old.values[i] = sol[0].values[i];
        }

        double rms = 0.0;
        for (int i = 0; i < mesh.nodes_num; ++i) {
            rms += pow(sol[0].values[i]
                - exact(mesh.nodes[i].x, mesh.nodes[i].y, mesh.nodes[i].z), 2);
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
