#include "integrands2d.h"
#include <cmath>

// See Zienkiewicz et al., section 5.11
double Integrand::Integrate (int triangle_index)
{
    const int POINTS_NUM = 4;
    const double POINTS_L0[POINTS_NUM] = {1./3., 0.6, 0.2, 0.2};
    const double POINTS_L1[POINTS_NUM] = {1./3., 0.2, 0.6, 0.2};
    const double WEIGHTS[POINTS_NUM] = {-27./48., 25./48., 25./48., 25./48.};

    double ans = 0.0;
    for (int i = 0; i < POINTS_NUM; ++i)
        ans += WEIGHTS[i] * Value(triangle_index, POINTS_L0[i], POINTS_L1[i]);
    return 2 * mesh.TriangleArea(triangle_index) * ans;
}

// Gauss quadrature, see Zienkiewicz et al., section 5.9.2
double BoundaryIntegrand::Integrate (int boundary_edge_index)
{
    const int POINTS_NUM = 3;
    const double POINTS[POINTS_NUM] = {- sqrt(0.6), 0, sqrt(0.6)};
    const double WEIGHTS[POINTS_NUM] = {5./9., 8./9., 5./9.};

    double ans = 0.0;
    for (int i = 0; i < POINTS_NUM; ++i)
        ans += WEIGHTS[i] * Value(boundary_edge_index, POINTS[i]);
    return 0.5 * mesh.BoundaryEdgeLength(boundary_edge_index) * ans;
}

double Integrate (Integrand& integrand)
{
    double ans = 0.0;
    for (auto triangle_index : integrand.support)
        ans += integrand.Integrate(triangle_index);
    return ans;
}

double BoundaryIntegrate (BoundaryIntegrand& integrand)
{
    double ans = 0.0;
    for (auto boundary_edge_index : integrand.boundary_support)
        ans += integrand.Integrate(boundary_edge_index);
    return ans;
}

std::list<int> intersect_supports (std::list<int>& s1, std::list<int>& s2)
{
    std::list<int> ans;
    for (auto i = s1.begin(); i != s1.end(); ++i) {
        if (std::find(s2.begin(), s2.end(), *i) != s2.end())
            ans.push_back(*i);
    }
    return ans;
}

Mult_P1_Basis operator* (FunctionP1& p1, BasisFunction& basis)
{
    return Mult_P1_Basis(p1.mesh, p1, basis);
}

Mult_P1_Basis_Basis operator* (Mult_P1_Basis& p1_basis, BasisFunction& basis)
{
    return Mult_P1_Basis_Basis(p1_basis.mesh, p1_basis.p1, p1_basis.basis, basis);
}

GradBasis grad (BasisFunction& basis)
{
    return GradBasis(basis);
}

Mult_GradBasis_GradBasis operator* (GradBasis& grad_basis1,
    GradBasis& grad_basis2)
{
    return Mult_GradBasis_GradBasis(grad_basis1.basis.mesh, grad_basis1.basis,
        grad_basis2.basis);
}

AuxBoundaryMult_P0_P1 operator* (BoundaryFunctionP0& p0, BoundaryFunctionP1& p1)
{
    return AuxBoundaryMult_P0_P1(p0, p1);
}

BoundaryMult_P0_P1_Basis operator* (AuxBoundaryMult_P0_P1& aux, BasisFunction& basis)
{
    return BoundaryMult_P0_P1_Basis(aux.p0.mesh, aux.p0, aux.p1, basis);
}

AuxBoundaryMult_P0_Basis operator* (BoundaryFunctionP0& p0,
    BasisFunction& basis)
{
    return AuxBoundaryMult_P0_Basis(p0, basis);
}

BoundaryMult_P0_Basis_Basis operator* (AuxBoundaryMult_P0_Basis& aux,
    BasisFunction& basis)
{
    return BoundaryMult_P0_Basis_Basis(aux.p0.mesh, aux.p0, aux.basis, basis);
}
