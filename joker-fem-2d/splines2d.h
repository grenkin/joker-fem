#ifndef SPLINES2D_H_INCLUDED
#define SPLINES2D_H_INCLUDED

#include "mesh2d.h"

// piecewise linear function in the domain
struct FunctionP1 {
    Mesh* mesh;
    std::vector<double> values;  // values of the function at mesh nodes

    FunctionP1 (Mesh& _mesh)
        : mesh(&_mesh)
    {
        values.resize(mesh->nodes_num);
    }

    // function value at the point with local coordinates L0, L1
    double Value (int triangle_index, double L0, double L1) const
    {
        int* nodes = mesh->triangles[triangle_index].nodes;
        return L0 * values[nodes[0]] + L1 * values[nodes[1]]
            + (1 - L0 - L1) * values[nodes[2]];
    }

    // function value at point (x, y)
    double ValueXY (double x, double y) const
    {
        int triangle_index = mesh->TriangleForPoint(x, y);
        if (triangle_index == -1)
            throw "The point does not belong to the domain.";
        double L0, L1;
        mesh->LocalCoordinates(triangle_index, x, y, L0, L1);
        return Value(triangle_index, L0, L1);
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
        values.resize(mesh->boundary_nodes_num);
    }

    // function value at the point with parameter t (from -1 to 1)
    double Value (int boundary_edge_index, double t) const
    {
        int* nodes = mesh->boundary_edges[boundary_edge_index].boundary_nodes;
        return 0.5 * (1 - t) * values[nodes[0]]
            + 0.5 * (1 + t) * values[nodes[1]];
    }
};

// piecewise constant function on the boundary
struct BoundaryFunctionP0 {
    Mesh* mesh;
    // values of the function at boundary edges of the mesh
    std::vector<double> values;

    BoundaryFunctionP0 (Mesh& _mesh)
        : mesh(&_mesh)
    {
        values.resize(mesh->boundary_edges_num);
    }

    double Value (int boundary_edge_index) const
    {
        return values[boundary_edge_index];
    }
};

#endif // SPLINES2D_H_INCLUDED
