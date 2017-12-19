#ifndef INTEGRANDS2D_H_INCLUDED
#define INTEGRANDS2D_H_INCLUDED

#include "splines2d.h"

class Integrand {
    Mesh& mesh;
public:
    // list of indices of triangles where the function is not zero
    std::list<int> support;
    // function value at the point with local coordinates L0, L1
    virtual double Value (int triangle_index, double L0, double L1) = 0;

    Integrand (Mesh& _mesh)
        : mesh(_mesh)
    {}

    // integrate the function by a given triangle
    double Integrate (int triangle_index);
};

class BoundaryIntegrand {
    Mesh& mesh;
public:
    // list of indices of boundary edges where the function is not zero
    std::list<int> boundary_support;
    // function value at the point with parameter t (from -1 to 1)
    virtual double Value (int boundary_edge_index, double t) = 0;

    BoundaryIntegrand (Mesh& _mesh)
        : mesh(_mesh)
    {}

    // integrate the function by a given boundary edge
    double Integrate (int boundary_edge_index);
};

struct BasisFunction {
    Mesh& mesh;
public:
    int node_index;
    // list of indices of triangles where the function is not zero
    std::list<int> support;
    // list of indices of boundary edges where the function is not zero
    std::list<int> boundary_support;
    std::vector<double> values;  // values of the function at mesh nodes
    // values of the function at boundary nodes of the mesh
    std::vector<double> boundary_values;

    BasisFunction (Mesh& _mesh, int _node_index)
        : mesh(_mesh), node_index(_node_index)
    {
        values.resize(mesh.nodes_num, 0.0);
        values[node_index] = 1.0;
        support = mesh.triangles_for_nodes[node_index];

        boundary_values.resize(mesh.boundary_nodes_num, 0.0);
        int boundary_node_index = mesh.boundary_indices[node_index];
        if (boundary_node_index != -1) {
            boundary_values[boundary_node_index] = 1.0;
            boundary_support =
                mesh.edges_for_boundary_nodes[boundary_node_index];
        }
    }

    // function value at the point with local coordinates L0, L1
    double Value (int triangle_index, double L0, double L1)
    {
        int* nodes = mesh.triangles[triangle_index].nodes;
        return L0 * values[nodes[0]] + L1 * values[nodes[1]]
            + (1 - L0 - L1) * values[nodes[2]];
    }

    // function value at the point with parameter t (from -1 to 1)
    double BoundaryValue (int boundary_edge_index, double t)
    {
        int* nodes = mesh.boundary_edges[boundary_edge_index].boundary_nodes;
        return 0.5 * (1 - t) * boundary_values[nodes[0]]
            + 0.5 * (1 + t) * boundary_values[nodes[1]];
    }

};

#endif // INTEGRANDS2D_H_INCLUDED
