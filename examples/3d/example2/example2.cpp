/*
-\Delta u + u - 2v = 0
-0.5\Delta v = 3 * sin(x) * sin(y) * sin(z)
-u_x + 3u|_{x=0} = 3 * (-1/3 * sin(y) * sin(z))
u_x + 3u|_{x=pi} = 3 * (-1/3 * sin(y) * sin(z))
-u_y + 3u|_{y=0} = 3 * (-1/3 * sin(x) * sin(z))
u_y + 3u|_{y=pi} = 3 * (-1/3 * sin(x) * sin(z))
-u_z + 3u|_{z=0} = 3 * (-1/3 * sin(x) * sin(y))
u_z + 3u|_{z=pi} = 3 * (-1/3 * sin(x) * sin(y))
-0.5v_x + v|_{x=0} = - sin(y) * sin(z)
0.5v_x + v|_{x=pi} = - sin(y) * sin(z)
-0.5v_y + v|_{y=0} = - sin(x) * sin(z)
0.5v_y + v|_{y=pi} = - sin(x) * sin(z)
-0.5v_z + v|_{z=0} = - sin(x) * sin(y)
0.5v_z + v|_{z=pi} = - sin(x) * sin(y)
Exact solution: u(x, y, z) = sin(x) * sin(y) * sin(z),
                v(x, y, z) = 2 * sin(x) * sin(y) * sin(z)
*/

#include <joker-fem-3d/bvp3d.h>
#include <iostream>
#include <cmath>

using namespace std;

double gfun2(double x, double y, double z)
{
    return 3 * sin(x) * sin(y) * sin(z);
}

double exact1(double x, double y, double z)
{
    return sin(x) * sin(y) * sin(z);
}

double exact2(double x, double y, double z)
{
    return 2 * sin(x) * sin(y) * sin(z);
}

bool eq (double a, double b)
{
    return fabs(a - b) < 1e-5;
}

int main() {
    int tests_num = 3;
    int n_values[tests_num] = {10, 20, 40};
    string files_names[tests_num] =
        {"tr10_cube.mesh", "tr20_cube.mesh", "tr40_cube.mesh"};

    double L = 3.14159265358979323846;
    int N = 2;

    double rms_old;
    for (int test = 0; test < tests_num; ++test) {
        Mesh mesh(files_names[test]);
        ProblemData data(N, mesh);
        data.a[0] = 1.0;  data.a[1] = 0.5;
        for (int i = 0; i < mesh.boundary_triangles_num; ++i) {
            data.b[0].values[i] = 3.0;
            data.b[1].values[i] = 1.0;
        }
        for (int i = 0; i < mesh.boundary_nodes_num; ++i) {
            double x = mesh.nodes[mesh.boundary_nodes[i]].x;
            double y = mesh.nodes[mesh.boundary_nodes[i]].y;
            double z = mesh.nodes[mesh.boundary_nodes[i]].z;
            if (eq(x, 0) || eq(x, L)) {
                data.w[0].values[i] = - 1./3. * sin(y) * sin(z);
                data.w[1].values[i] = - sin(y) * sin(z);
            }
            else if (eq(y, 0) || eq(y, L)) {
                data.w[0].values[i] = - 1./3. * sin(x) * sin(z);
                data.w[1].values[i] = - sin(x) * sin(z);
            }
            else if (eq(z, 0) || eq(z, L)) {
                data.w[0].values[i] = - 1./3. * sin(x) * sin(y);
                data.w[1].values[i] = - sin(x) * sin(y);
            }
            else
                throw;
        }
        for (int i = 0; i < mesh.nodes_num; ++i) {
            data.g[0].values[i] = 0.0;
            data.g[1].values[i] =
                gfun2(mesh.nodes[i].x, mesh.nodes[i].y, mesh.nodes[i].z);
            data.q[0][0].values[i] = 1.0;
            data.q[0][1].values[i] = -2.0;
            data.q[1][0].values[i] = data.q[1][1].values[i] = 0.0;
        }

        vector<FunctionP1> sol(N, mesh);
        for (int i = 0; i < mesh.nodes_num; ++i) {
            sol[0].values[i] = 0.0;
            sol[1].values[i] = 0.0;
        }
        SolveBVP(data, Parameters(), sol);

        double rms = 0.0;
        for (int i = 0; i < mesh.nodes_num; ++i) {
            rms += pow(sol[0].values[i]
                    - exact1(mesh.nodes[i].x, mesh.nodes[i].y, mesh.nodes[i].z), 2)
                + pow(sol[1].values[i]
                    - exact2(mesh.nodes[i].x, mesh.nodes[i].y, mesh.nodes[i].z), 2);
        }
        rms = sqrt(rms / (2 * mesh.nodes_num));
        cout << "n = " << n_values[test] << "   rms = " << rms << endl;
        if (test > 0)
            cout << "rms_old / rms = " << rms_old / rms << endl;
        cout << endl;
        rms_old = rms;
    }

    return 0;
}
