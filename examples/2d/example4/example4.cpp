/*
-2\Delta u - 2v = 2sin(y)
-\Delta v = sin(x)
-2u_x + u|_{x=0} = -2 + sin(y)
2u_x + 3u|_{x=pi} = 3(-2/3 + sin(y))
-2u_y + 5u|_{y=0} = 5(-2/5 + sin(x))
2u_y + u|_{y=pi} = -2 + sin(x)
-v_x + 3v|_{x=0} = 3 * (-1/3)
v_x + 3v|_{x=pi} = 3 * (-1/3)
-v_y|{y=0} = 0
v_y{y=pi} = 0
Exact solution: u(x, y) = sin(x) + sin(y), v(x, y) = sin(x)
*/

#include <joker-fem-2d/bvp2d.h>
#include <iostream>
#include <cmath>

using namespace std;

double gfun1(double x, double y)
{
    return 2 * sin(y);
}

double gfun2(double x, double y)
{
    return sin(x);
}

double exact1(double x, double y)
{
    return sin(x) + sin(y);
}

double exact2(double x, double y)
{
    return sin(x);
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
    int N = 2;

    double rms_old;
    for (int test = 0; test < tests_num; ++test) {
        Mesh mesh(files_names[test]);
        ProblemData data(N, mesh);
        data.a[0] = 2.0;
        data.a[1] = 1.0;
        for (int i = 0; i < mesh.boundary_edges_num; ++i) {
            int* boundary_nodes = mesh.boundary_edges[i].boundary_nodes;
            double x1 = mesh.nodes[mesh.boundary_nodes[boundary_nodes[0]]].x;
            double y1 = mesh.nodes[mesh.boundary_nodes[boundary_nodes[0]]].y;
            double x2 = mesh.nodes[mesh.boundary_nodes[boundary_nodes[1]]].x;
            double y2 = mesh.nodes[mesh.boundary_nodes[boundary_nodes[1]]].y;
            if (eq(x1, 0) && eq(x2, 0)) {
                data.b[0].values[i] = 1.0;
                data.b[1].values[i] = 3.0;
            }
            else if (eq(x1, L) && eq(x2, L)) {
                data.b[0].values[i] = 3.0;
                data.b[1].values[i] = 3.0;
            }
            else if (eq(y1, 0) && eq(y2, 0)) {
                data.b[0].values[i] = 5.0;
                data.b[1].values[i] = 0.0;
            }
            else if (eq(y1, L) && eq(y2, L)) {
                data.b[0].values[i] = 1.0;
                data.b[1].values[i] = 0.0;
            }
            else
                throw;
        }
        for (int i = 0; i < mesh.boundary_nodes_num; ++i) {
            double x = mesh.nodes[mesh.boundary_nodes[i]].x;
            double y = mesh.nodes[mesh.boundary_nodes[i]].y;
            if (eq(x, 0)) {
                data.w[0].values[i] = -2 + sin(y);
                data.w[1].values[i] = -1./3.;
            }
            else if (eq(x, L)) {
                data.w[0].values[i] = -2./3. + sin(y);
                data.w[1].values[i] = -1./3.;
            }
            else if (eq(y, 0)) {
                data.w[0].values[i] = -2./5. + sin(x);
                data.w[1].values[i] = 1.0;
            }
            else if (eq(y, L)) {
                data.w[0].values[i] = -2 + sin(x);
                data.w[1].values[i] = 1.0;
            }
            else
                throw;
        }
        for (int i = 0; i < mesh.nodes_num; ++i) {
            data.g[0].values[i] = gfun1(mesh.nodes[i].x, mesh.nodes[i].y);
            data.g[1].values[i] = gfun2(mesh.nodes[i].x, mesh.nodes[i].y);
            data.q[0][0].values[i] = 0.0;
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
            rms += pow(sol[0].values[i] - exact1(mesh.nodes[i].x, mesh.nodes[i].y), 2)
                + pow(sol[1].values[i] - exact2(mesh.nodes[i].x, mesh.nodes[i].y), 2);
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
