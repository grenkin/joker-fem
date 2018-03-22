#include "integrands3d.h"
#include <cmath>
#include <algorithm>

double Integrand::Integrate (int tetrahedron_index) const
{
    // See Zienkiewicz et al., section 5.11
    const int POINTS_NUM = 5;
    const double POINTS_L0[POINTS_NUM] = {1./4., 1./2., 1./6., 1./6., 1./6.};
    const double POINTS_L1[POINTS_NUM] = {1./4., 1./6., 1./2., 1./6., 1./6.};
    const double POINTS_L2[POINTS_NUM] = {1./4., 1./6., 1./6., 1./2., 1./6.};
    const double WEIGHTS[POINTS_NUM] = {-4./5., 9./20., 9./20., 9./20., 9./20.};

    double ans = 0.0;
    for (int i = 0; i < POINTS_NUM; ++i)
        ans += WEIGHTS[i] * Value(
            tetrahedron_index, POINTS_L0[i], POINTS_L1[i], POINTS_L2[i]);
    return mesh->TetrahedronVolume(tetrahedron_index) * ans;
}

double BoundaryIntegrand::Integrate (int boundary_triangle_index) const
{
    // See Zienkiewicz et al., section 5.11
    const int POINTS_NUM = 4;
    const double POINTS_L0[POINTS_NUM] = {1./3., 0.6, 0.2, 0.2};
    const double POINTS_L1[POINTS_NUM] = {1./3., 0.2, 0.6, 0.2};
    const double WEIGHTS[POINTS_NUM] = {-27./48., 25./48., 25./48., 25./48.};

    double ans = 0.0;
    for (int i = 0; i < POINTS_NUM; ++i)
        ans += WEIGHTS[i] * Value(boundary_triangle_index,
            POINTS_L0[i], POINTS_L1[i]);
    return mesh->BoundaryTriangleArea(boundary_triangle_index) * ans;
}

double Integrate (const Integrand& integrand)
{
    double ans = 0.0;
    for (auto tetrahedron_index : integrand.support)
        ans += integrand.Integrate(tetrahedron_index);
    return ans;
}

double BoundaryIntegrate (const BoundaryIntegrand& integrand)
{
    double ans = 0.0;
    for (auto boundary_triangle_index : integrand.boundary_support)
        ans += integrand.Integrate(boundary_triangle_index);
    return ans;
}

std::list<int> intersect_supports (const std::list<int>& s1,
    const std::list<int>& s2)
{
    std::list<int> ans;
    for (auto i = s1.begin(); i != s1.end(); ++i) {
        if (std::find(s2.begin(), s2.end(), *i) != s2.end())
            ans.push_back(*i);
    }
    return ans;
}

Mult_P1_Basis operator* (const FunctionP1& p1, const BasisFunction& basis)
{
    return Mult_P1_Basis(p1.mesh, p1, basis);
}

Mult_P1_Basis_Basis operator* (const Mult_P1_Basis& p1_basis,
    const BasisFunction& basis)
{
    return Mult_P1_Basis_Basis(p1_basis.mesh, p1_basis.p1,
        p1_basis.basis, basis);
}

AuxGradBasis grad (const BasisFunction& basis)
{
    return AuxGradBasis(basis);
}

Mult_GradBasis_GradBasis operator* (const AuxGradBasis& grad_basis1,
    const AuxGradBasis& grad_basis2)
{
    return Mult_GradBasis_GradBasis(grad_basis1.basis.mesh, grad_basis1.basis,
        grad_basis2.basis);
}

AuxBoundaryMult_P0_P1 operator* (const BoundaryFunctionP0& p0,
    const BoundaryFunctionP1& p1)
{
    return AuxBoundaryMult_P0_P1(p0, p1);
}

BoundaryMult_P0_P1_Basis operator* (const AuxBoundaryMult_P0_P1& aux,
    const BasisFunction& basis)
{
    return BoundaryMult_P0_P1_Basis(aux.p0.mesh, aux.p0, aux.p1, basis);
}

AuxBoundaryMult_P0_Basis operator* (const BoundaryFunctionP0& p0,
    const BasisFunction& basis)
{
    return AuxBoundaryMult_P0_Basis(p0, basis);
}

BoundaryMult_P0_Basis_Basis operator* (const AuxBoundaryMult_P0_Basis& aux,
    const BasisFunction& basis)
{
    return BoundaryMult_P0_Basis_Basis(aux.p0.mesh, aux.p0, aux.basis, basis);
}
