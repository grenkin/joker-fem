#ifndef SPLINES3D_H_INCLUDED
#define SPLINES3D_H_INCLUDED

#include "mesh3d.h"

// piecewise linear function in the domain
struct FunctionP1 {
    Mesh* mesh;
    std::vector<double> values;  // values of the function at mesh nodes

    FunctionP1 (Mesh& _mesh)
        : mesh(&_mesh)
    {
        values.resize(mesh->nodes_num);
    }

    // function value at the point with local coordinates L0, L1, L2
    double Value (int tetrahedron_index, double L0, double L1, double L2) const
    {
        int* nodes = mesh->tetrahedrons[tetrahedron_index].nodes;
        return L0 * values[nodes[0]] + L1 * values[nodes[1]]
            + L2 * values[nodes[2]] + (1 - L0 - L1 - L2) * values[nodes[3]];
    }

    // function value at point (x, y, z)
    double ValueXYZ (double x, double y, double z) const
    {
        int tetrahedron_index = mesh->TetrahedronForPoint(x, y, z);
        if (tetrahedron_index == -1)
            throw "The point does not belong to the domain.";
        double L0, L1, L2;
        mesh->LocalCoordinates(tetrahedron_index, x, y, z, L0, L1, L2);
        return Value(tetrahedron_index, L0, L1, L2);
    }
};

// piecewise linear function on the boundary
struct BoundaryFunctionP1 {
    Mesh* mesh;
    // values of the function at boundary nodes of the mesh
    std::vector<double> values;

    BoundaryFunctionP1 (Mesh& _mesh)
        : mesh(&_mesh)
    {
        values.resize(mesh->boundary_triangles_num);
    }

    // function value at the point with local coordinates L0, L1
    double Value (int boundary_triangle_index, double L0, double L1) const
    {
        int* nodes =
            mesh->boundary_triangles[boundary_triangle_index].boundary_nodes;
        return L0 * values[nodes[0]] + L1 * values[nodes[1]]
            + (1 - L0 - L1) * values[nodes[2]];
    }
};

// piecewise constant function on the boundary
struct BoundaryFunctionP0 {
    Mesh* mesh;
    // values of the function at boundary triangles of the mesh
    std::vector<double> values;

    BoundaryFunctionP0 (Mesh& _mesh)
        : mesh(&_mesh)
    {
        values.resize(mesh->boundary_triangles_num);
    }

    double Value (int boundary_triangle_index) const
    {
        return values[boundary_triangle_index];
    }
};

#endif // SPLINES3D_H_INCLUDED
