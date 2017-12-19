#ifndef MESH2D_H_INCLUDED
#define MESH2D_H_INCLUDED

#include <vector>
#include <list>

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
    // for each node - a list of indices of triangles containing this node
    std::vector<std::list<int> > triangles_for_nodes;

    Mesh (std::string file_name);
};

#endif // MESH2D_H_INCLUDED
