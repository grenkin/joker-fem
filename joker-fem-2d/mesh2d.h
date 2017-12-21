#ifndef MESH2D_H_INCLUDED
#define MESH2D_H_INCLUDED

#include <vector>
#include <list>
#include <string>

struct Node {
    double x, y;
};

struct Triangle {
    // indices of triangle vertices in the vector of nodes
    int nodes[3];
};

struct BoundaryEdge {
    // indices of boundary edge vertices in the vector of boundary nodes
    int boundary_nodes[2];
};

struct Mesh {
    // numbers of mesh elements
    int nodes_num, triangles_num, boundary_edges_num;
    std::vector<Node> nodes;
    std::vector<Triangle> triangles;
    std::vector<BoundaryEdge> boundary_edges;

    int boundary_nodes_num;
    // indices of boundary nodes in the vector of nodes
    std::vector<int> boundary_nodes;
    // indices of nodes in the vector boundary_nodes (-1 for non-boundary nodes)
    std::vector<int> boundary_indices;
    // for each node - a list of indices of triangles containing this node
    std::vector<std::list<int> > triangles_for_nodes;
    // for each boundary node - a list of boundary edges containing this node
    std::vector<std::list<int> > edges_for_boundary_nodes;
    // for each node - a list of nodes connected with this node by an edge
    // and the node itself
    std::vector<std::list<int> > adjacent_nodes;

    Mesh (std::string file_name);
    double SignedTriangleArea (int triangle_index);
    double TriangleArea (int triangle_index);
    double BoundaryEdgeLength (int boundary_edge_index);
    void LocalCoefficients (int triangle_index,
        double (&a)[3], double (&b)[3], double (&c)[3]);
    // local coordinates of a given point
    void LocalCoordinates (int triangle_index, double x, double y,
        double& L0, double& L1);
    bool PointInTriangle (int triangle_index, double x, double y);
    // index of the triangle containing a given point
    int TriangleForPoint (double x, double y);
};

#endif // MESH2D_H_INCLUDED
