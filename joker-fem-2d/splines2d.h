#ifndef SPLINES2D_H_INCLUDED
#define SPLINES2D_H_INCLUDED

#include "mesh2d.h"

// piecewise linear function in the domain
class FunctionP1 {
    Mesh& mesh;
public:
    std::vector<double> values;  // values of the function at mesh nodes

    FunctionP1 (Mesh& _mesh)
        : mesh(_mesh)
    {
        values.resize(mesh.nodes_num);
    }

    // function value at the point with local coordinates L0, L1
    double Value (int triangle_index, double L0, double L1)
    {
        int* nodes = mesh.triangles[triangle_index].nodes;
        return L0 * values[nodes[0]] + L1 * values[nodes[1]]
            + (1 - L0 - L1) * values[nodes[2]];
    }
};

// piecewise linear function on the boundary
class BoundaryFunctionP1 {
    Mesh& mesh;
public:
    // values of the function at boundary nodes of the mesh
    std::vector<double> values;

    BoundaryFunctionP1 (Mesh& _mesh)
        : mesh(_mesh)
    {
        values.resize(mesh.boundary_nodes_num);
    }

    // function value at the point with parameter t (from -1 to 1)
    double Value (int boundary_edge_index, double t)
    {
        int* nodes = mesh.boundary_edges[boundary_edge_index].boundary_nodes;
        return 0.5 * (1 - t) * values[nodes[0]]
            + 0.5 * (1 + t) * values[nodes[1]];
    }
};

// piecewise constant function on the boundary
class BoundaryFunctionP0 {
    Mesh& mesh;
public:
    // values of the function at boundary edges of the mesh
    std::vector<double> values;

    BoundaryFunctionP0 (Mesh& _mesh)
        : mesh(_mesh)
    {
        values.resize(mesh.boundary_edges_num);
    }

    double Value (int boundary_edge_index)
    {
        return values[boundary_edge_index];
    }
};

#endif // SPLINES2D_H_INCLUDED