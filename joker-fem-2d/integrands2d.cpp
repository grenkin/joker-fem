#include "integrands2d.h"
#include <cmath>

// See Zienkiewicz et al., section 5.11
double Integrand::Integrate (int triangle_index)
{
    const int POINTS_NUM = 4;
    const double POINTS_L0[POINTS_NUM] = {1./3., 0.6, 0.2, 0.2};
    const double POINTS_L1[POINTS_NUM] = {1./3., 0.2, 0.6, 0.2};
    const double WEIGHTS[POINTS_NUM] = {-27./48., 25./48., 25./48., 25./48.};

    double ans = 0;
    for (int i = 0; i < POINTS_NUM; ++i)
        ans += WEIGHTS[i] * Value(triangle_index, POINTS_L0[i], POINTS_L1[i]);
    return ans;
}

// Gauss quadrature, see Zienkiewicz et al., section 5.9.2
double BoundaryIntegrand::Integrate (int boundary_edge_index)
{
    const int POINTS_NUM = 3;
    const double POINTS[POINTS_NUM] = {- sqrt(0.6), 0, sqrt(0.6)};
    const double WEIGHTS[POINTS_NUM] = {5./9., 8./9., 5./9.};

    double ans = 0;
    for (int i = 0; i < POINTS_NUM; ++i)
        ans += WEIGHTS[i] * Value(boundary_edge_index, POINTS[i]);
    return ans;
}

double Integrate (Integrand& integrand)
{
    double ans = 0.0;
    for (std::list<int>::iterator i = integrand.support.begin();
        i != integrand.support.end(); ++i)
    {
        int triangle_index = *i;
        ans += integrand.Integrate(triangle_index);
    }
    return ans;
}

double BoundaryIntegrate (BoundaryIntegrand& integrand)
{
    double ans = 0.0;
    for (std::list<int>::iterator i = integrand.boundary_support.begin();
        i != integrand.boundary_support.end(); ++i)
    {
        int boundary_edge_index = *i;
        ans += integrand.Integrate(boundary_edge_index);
    }
    return ans;
}

std::list<int> intersect_supports (std::list<int>& s1, std::list<int>& s2)
{
    std::list<int> ans;
    for (std::list<int>::iterator i = s1.begin(); i != s1.end(); ++i) {
        if (std::find(s2.begin(), s2.end(), *i) != s2.end())
            ans.push_back(*i);
    }
    return ans;
}

Mult_P1_Basis operator* (FunctionP1& p1, BasisFunction& basis)
{
    return Mult_P1_Basis(p1.mesh, p1, basis);
}

BoundaryMult_P0_Basis_Basis operator* (AuxBoundaryMult_P0_Basis& aux,
    BasisFunction& basis)
{
    return BoundaryMult_P0_Basis_Basis(aux.mesh, aux.p0, aux.basis, basis);
}
