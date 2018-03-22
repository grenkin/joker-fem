#ifndef INTEGRANDS3D_H_INCLUDED
#define INTEGRANDS3D_H_INCLUDED

#include "splines3d.h"

struct Integrand {
    Mesh* mesh;
    // list of indices of tetrahedrons where the function is not zero
    std::list<int> support;
    // function value at the point with local coordinates L0, L1, L2
    virtual double Value (int tetrahedron_index,
        double L0, double L1, double L2) const = 0;

    Integrand (Mesh* _mesh)
        : mesh(_mesh)
    {}

    // integrate the function over a given tetrahedron
    double Integrate (int tetrahedron_index) const;
};

struct BoundaryIntegrand {
    Mesh* mesh;
    // list of indices of boundary triangles where the function is not zero
    std::list<int> boundary_support;
    // function value at the point with local coordinates L0, L1
    virtual double Value (int boundary_triangle_index,
        double L0, double L1) const = 0;

    BoundaryIntegrand (Mesh* _mesh)
        : mesh(_mesh)
    {}

    // integrate the function over a given boundary triangle
    double Integrate (int boundary_triangle_index) const;
};

// integrate the function over the domain
double Integrate (const Integrand& integrand);
// integrate the function over the boundary
double BoundaryIntegrate (const BoundaryIntegrand& integrand);

struct Vector3 {
    double x, y, z;

    Vector3 (double _x, double _y, double _z)
        : x(_x), y(_y), z(_z)
    {}

    double operator* (Vector3 vec)
    {
        return x * vec.x + y * vec.y + z * vec.z;
    }
};

struct BasisFunction {
    Mesh* mesh;
    int node_index;
    // list of indices of tetrahedrons where the function is not zero
    std::list<int> support;
    // list of indices of boundary triangles where the function is not zero
    std::list<int> boundary_support;

    BasisFunction (Mesh& _mesh)
        : mesh(&_mesh), node_index(-1)
    {}

    BasisFunction (Mesh& _mesh, int _node_index)
        : mesh(&_mesh), node_index(_node_index)
    {
        support = mesh->tetrahedrons_for_nodes[node_index];

        int boundary_node_index = mesh->boundary_indices[node_index];
        if (boundary_node_index != -1) {
            boundary_support =
                mesh->triangles_for_boundary_nodes[boundary_node_index];
        }
    }

    // function value at the point with local coordinates L0, L1, L2
    double Value (int tetrahedron_index, double L0, double L1, double L2) const
    {
        double values[4];  // function values at the triangle vertices
        for (int i = 0; i < 4; ++i) {
            values[i] =
                (mesh->tetrahedrons[tetrahedron_index].nodes[i] == node_index)
                    ? 1.0 : 0.0;
        }
        return L0 * values[0] + L1 * values[1] + L2 * values[2]
            + (1 - L0 - L1 - L2) * values[3];
    }

    // function value at the point with local coordinates L0, L1
    double BoundaryValue (int boundary_triangle_index,
        double L0, double L1) const
    {
        double values[3];  // function values at the triangle vertices
        for (int i = 0; i < 3; ++i) {
            values[i] =
                (mesh->boundary_triangles[boundary_triangle_index].boundary_nodes[i]
                    == mesh->boundary_indices[node_index]) ? 1.0 : 0.0;
        }
        return L0 * values[0] + L1 * values[1] + (1 - L0 - L1) * values[2];
    }

    // function gradient value at a given triangle
    Vector3 GradValue (int tetrahedron_index) const
    {
        double values[4];  // function values at the tetrahedron vertices
        for (int i = 0; i < 4; ++i) {
            values[i] =
                (mesh->tetrahedrons[tetrahedron_index].nodes[i] == node_index)
                    ? 1.0 : 0.0;
        }
        double a[4], b[4], c[4], d[4];
        double A = mesh->SignedTetrahedronVolume(tetrahedron_index);
        mesh->LocalCoefficients(tetrahedron_index, a, b, c, d);
        return Vector3(
            1 / (6 * A) * (values[0] * b[0] + values[1] * b[1] + values[2] * b[2]
                + values[3] * b[3]),
            1 / (6 * A) * (values[0] * c[0] + values[1] * c[1] + values[2] * c[2]
                + values[3] * c[3]),
            1 / (6 * A) * (values[0] * d[0] + values[1] * d[1] + values[2] * d[2]
                + values[3] * d[3]));
    }
};

std::list<int> intersect_supports (const std::list<int>&,
    const std::list<int>&);

struct Mult_P1_Basis : public Integrand {
    const FunctionP1& p1;
    const BasisFunction& basis;

    Mult_P1_Basis (Mesh* _mesh, const FunctionP1& _p1,
        const BasisFunction& _basis)
        : Integrand(_mesh), p1(_p1), basis(_basis)
    {
        support = basis.support;
    }

    double Value (int tetrahedron_index, double L0, double L1, double L2) const
    {
        return p1.Value(tetrahedron_index, L0, L1, L2)
            * basis.Value(tetrahedron_index, L0, L1, L2);
    }
};

Mult_P1_Basis operator* (const FunctionP1&, const BasisFunction&);

struct Mult_P1_Basis_Basis : public Integrand {
    const FunctionP1& p1;
    const BasisFunction& basis1;
    const BasisFunction& basis2;

    Mult_P1_Basis_Basis (Mesh* _mesh, const FunctionP1& _p1,
        const BasisFunction& _basis1, const BasisFunction& _basis2)
        : Integrand(_mesh), p1(_p1), basis1(_basis1), basis2(_basis2)
    {
        support = intersect_supports(basis1.support, basis2.support);
    }

    double Value (int tetrahedron_index, double L0, double L1, double L2) const
    {
        return p1.Value(tetrahedron_index, L0, L1, L2)
            * basis1.Value(tetrahedron_index, L0, L1, L2)
            * basis2.Value(tetrahedron_index, L0, L1, L2);
    }
};

Mult_P1_Basis_Basis operator* (const Mult_P1_Basis&, const BasisFunction&);

struct AuxGradBasis {
    const BasisFunction& basis;

    AuxGradBasis (const BasisFunction& _basis)
        : basis(_basis)
    {}
};

AuxGradBasis grad (const BasisFunction&);

struct Mult_GradBasis_GradBasis : public Integrand {
    const BasisFunction& basis1;
    const BasisFunction& basis2;

    Mult_GradBasis_GradBasis (Mesh* _mesh, const BasisFunction& _basis1,
        const BasisFunction& _basis2)
        : Integrand(_mesh), basis1(_basis1), basis2(_basis2)
    {
        support = intersect_supports(basis1.support, basis2.support);
    }

    double Value (int tetrahedron_index, double L0, double L1, double L2) const
    {
        return basis1.GradValue(tetrahedron_index)
            * basis2.GradValue(tetrahedron_index);
    }
};

Mult_GradBasis_GradBasis operator* (const AuxGradBasis&, const AuxGradBasis&);

struct AuxBoundaryMult_P0_P1 {
    const BoundaryFunctionP0& p0;
    const BoundaryFunctionP1& p1;

    AuxBoundaryMult_P0_P1 (const BoundaryFunctionP0& _p0,
        const BoundaryFunctionP1& _p1)
        : p0(_p0), p1(_p1)
    {}
};

AuxBoundaryMult_P0_P1 operator* (const BoundaryFunctionP0&,
    const BoundaryFunctionP1&);

struct BoundaryMult_P0_P1_Basis : public BoundaryIntegrand {
    const BoundaryFunctionP0& p0;
    const BoundaryFunctionP1& p1;
    const BasisFunction& basis;

    BoundaryMult_P0_P1_Basis (Mesh* _mesh, const BoundaryFunctionP0& _p0,
        const BoundaryFunctionP1& _p1, const BasisFunction& _basis)
        : BoundaryIntegrand(_mesh), p0(_p0), p1(_p1), basis(_basis)
    {
        boundary_support = basis.boundary_support;
    }

    double Value (int boundary_triangle_index, double L0, double L1) const
    {
        return p0.Value(boundary_triangle_index)
            * p1.Value(boundary_triangle_index, L0, L1)
            * basis.BoundaryValue(boundary_triangle_index, L0, L1);
    }
};

BoundaryMult_P0_P1_Basis operator* (const AuxBoundaryMult_P0_P1&,
    const BasisFunction&);

struct AuxBoundaryMult_P0_Basis {
    const BoundaryFunctionP0& p0;
    const BasisFunction& basis;

    AuxBoundaryMult_P0_Basis (const BoundaryFunctionP0& _p0,
        const BasisFunction& _basis)
        : p0(_p0), basis(_basis)
    {}
};

AuxBoundaryMult_P0_Basis operator* (const BoundaryFunctionP0&,
    const BasisFunction&);

struct BoundaryMult_P0_Basis_Basis : public BoundaryIntegrand {
    const BoundaryFunctionP0& p0;
    const BasisFunction& basis1;
    const BasisFunction& basis2;

    BoundaryMult_P0_Basis_Basis (Mesh* _mesh, const BoundaryFunctionP0& _p0,
        const BasisFunction& _basis1, const BasisFunction& _basis2)
        : BoundaryIntegrand(_mesh), p0(_p0), basis1(_basis1), basis2(_basis2)
    {
        boundary_support = intersect_supports(basis1.boundary_support,
            basis2.boundary_support);
    }

    double Value (int boundary_triangle_index, double L0, double L1) const
    {
        return p0.Value(boundary_triangle_index)
            * basis1.BoundaryValue(boundary_triangle_index, L0, L1)
            * basis2.BoundaryValue(boundary_triangle_index, L0, L1);
    }
};

BoundaryMult_P0_Basis_Basis operator* (const AuxBoundaryMult_P0_Basis&,
    const BasisFunction&);

#endif // INTEGRANDS3D_H_INCLUDED
