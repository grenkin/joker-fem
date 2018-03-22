#ifndef MESH3D_H_INCLUDED
#define MESH3D_H_INCLUDED

#include <vector>
#include <list>
#include <string>

struct Node {
    double x, y, z;
};

struct Tetrahedron {
    // indices of tetrahedron vertices in the vector of nodes
    int nodes[4];
};

struct BoundaryTriangle {
    // indices of boundary triangle vertices in the vector of boundary nodes
    int boundary_nodes[3];
};

class Mesh {
    void add_boundary_node (int node_index);
    void add_adjacent_node (int node, int node_to_add);
public:
    // numbers of mesh elements
    int nodes_num, tetrahedrons_num, boundary_triangles_num;
    std::vector<Node> nodes;
    std::vector<Tetrahedron> tetrahedrons;
    std::vector<BoundaryTriangle> boundary_triangles;

    int boundary_nodes_num;
    // indices of boundary nodes in the vector of nodes
    std::vector<int> boundary_nodes;
    // indices of nodes in the vector boundary_nodes (-1 for non-boundary nodes)
    std::vector<int> boundary_indices;
    // for each node - a list of indices of tetrahedrons containing this node
    std::vector<std::list<int> > tetrahedrons_for_nodes;
    // for each boundary node - a list of boundary triangles containing this node
    std::vector<std::list<int> > triangles_for_boundary_nodes;
    // for each node - a list of nodes connected with this node by an edge
    // and the node itself
    std::vector<std::list<int> > adjacent_nodes;

    Mesh (std::string file_name);
    double SignedTetrahedronVolume (int tetrahedron_index);
    double TetrahedronVolume (int tetrahedron_index);
    double BoundaryTriangleArea (int boundary_triangle_index);
    void LocalCoefficients (int tetrahedron_index,
        double (&a)[4], double (&b)[4], double (&c)[4], double (&d)[4]);
    // local coordinates of a given point
    void LocalCoordinates (int tetrahedron_index, double x, double y, double z,
        double& L0, double& L1, double& L2);
    bool PointInTetrahedron (int tetrahedron_index,
        double x, double y, double z);
    // index of the tetrahedron containing a given point
    int TetrahedronForPoint (double x, double y, double z);
};

#endif // MESH3D_H_INCLUDED
