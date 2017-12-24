/*
Steady-state radiative heat transfer model

-a\Delta theta + b * kappaa * (theta^4 - phi) = 0
-alpha\Delta phi + kappaa * (phi - theta^4) = 0
a * dtheta/dn + beta * (theta - theta_b) = 0
alpha * dphi / dn + gamma * (phi - theta_b^4) = 0
0 < x < L,  0 < y < L
theta_b(0, y) = theta_b(L, y) = 0.5 + 0.5 * y / L
theta_b(x, 0) = 0.5,  theta_b(x, L) = 1
*/

#include <joker-fem-2d/bvp2d.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double Tmax = 773;
double StBol = 5.67e-8;
double a = 0.0515;
double beta = 10;
double kappa = 1;
double kappafull = 10;
double alpha = 1. / (3 * kappafull);
double n = 1;
double b = 4 * StBol * pow(Tmax, 3) * beta * pow(n, 2);
double gamma = 0.3;
double L = 1;

bool eq (double a, double b)
{
    return fabs(a - b) < 1e-5;
}

int main() {
    int N = 2;
    Mesh mesh("tr.mesh");
    ProblemData data(N, mesh);
    data.a[0] = a;
    data.a[1] = alpha;
    for (int i = 0; i < mesh.boundary_edges_num; ++i) {
        data.b[0].values[i] = beta;
        data.b[1].values[i] = gamma;
    }
    for (int i = 0; i < mesh.boundary_nodes_num; ++i) {
        double x = mesh.nodes[mesh.boundary_nodes[i]].x;
        double y = mesh.nodes[mesh.boundary_nodes[i]].y;
        double thetab;
        if (eq(x, 0) || eq(x, L))
            thetab = 0.5 + 0.5 * y / L;
        else if (eq(y, 0))
            thetab = 0.5;
        else if (eq(y, L))
            thetab = 1.;
        else
            throw;
        data.w[0].values[i] = thetab;
        data.w[1].values[i] = pow(thetab, 4);
    }
    vector<FunctionP1> sol(N, mesh), sol_old(N, mesh);
    // set the initial guess
    for (int i = 0; i < mesh.nodes_num; ++i) {
        sol[0].values[i] = 1.0;
        sol[1].values[i] = 1.0;
        sol_old[0].values[i] = 1.0;
        sol_old[1].values[i] = 1.0;
    }
    for (int iter = 1; iter <= 8; ++iter) {
        cout << iter;
        /* Linearized equations:
        -a\Delta theta + b * kappaa * (4*theta0^3 * theta - phi) = b * kappaa * 3*theta0^4
        -alpha\Delta phi + kappaa * (phi - 4*theta0^3 * theta) = - kappaa * 3*theta0^4
        */
        for (int i = 0; i < mesh.nodes_num; ++i) {
            data.g[0].values[i] = b * kappa * 3 * pow(sol_old[0].values[i], 4);
            data.g[1].values[i] = - kappa * 3 * pow(sol_old[0].values[i], 4);
            data.q[0][0].values[i] = b * kappa * 4 * pow(sol_old[0].values[i], 3);
            data.q[0][1].values[i] = - b * kappa;
            data.q[1][0].values[i] = - kappa * 4 * pow(sol_old[0].values[i], 3);
            data.q[1][1].values[i] = kappa;
        }

        SolveBVP(data, Parameters(), sol);

        double max_diff = 0.0;
        for (int i = 0; i < mesh.nodes_num; ++i)
            max_diff = fmax(max_diff, fabs(sol[0].values[i] - sol_old[0].values[i]));
        cout << "   diff = " << max_diff << endl;

        for (int i = 0; i < mesh.nodes_num; ++i) {
            sol_old[0].values[i] = sol[0].values[i];
            sol_old[1].values[i] = sol[1].values[i];
        }
    }

    ofstream fout("output.txt");
    for (int i = 0; i <= 100; ++i) {
        for (int j = 0; j <= 100; ++j) {
            double x = i * L / 100;
            double y = j * L / 100;
            fout << x << "  " << y << "  " << sol_old[1].ValueXY(x, y)
                << "  " << sol_old[0].ValueXY(x, y) << endl;
        }
    }

    return 0;
}
