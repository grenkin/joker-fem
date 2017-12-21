#include "mesh2d.h"
#include <cmath>

Mesh::Mesh (std::string file_name)
{
    // TODO
}

double det_3_3 (double a11, double a12, double a13,
    double a21, double a22, double a23, double a31, double a32, double a33)
{
    return a11 * a22 * a33  +  a12 * a23 * a31  +  a13 * a21 * a32
        -  a13 * a22 * a31  -  a11 * a23 * a32  -  a12 * a21 * a33;
}

double Mesh::SignedTriangleArea (int triangle_index)
{
    int* nodes_ind = triangles[triangle_index].nodes;
    return 0.5 * det_3_3(
        1, nodes[nodes_ind[0]].x, nodes[nodes_ind[0]].y,
        1, nodes[nodes_ind[1]].x, nodes[nodes_ind[1]].y,
        1, nodes[nodes_ind[2]].x, nodes[nodes_ind[2]].y
    );
}

double Mesh::TriangleArea (int triangle_index)
{
    return fabs(SignedTriangleArea(triangle_index));
}

double Mesh::BoundaryEdgeLength (int boundary_edge_index)
{
    int* nodes_ind = boundary_edges[boundary_edge_index].boundary_nodes;
    double x1 = nodes[boundary_nodes[nodes_ind[0]]].x;
    double y1 = nodes[boundary_nodes[nodes_ind[0]]].y;
    double x2 = nodes[boundary_nodes[nodes_ind[1]]].x;
    double y2 = nodes[boundary_nodes[nodes_ind[1]]].y;
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

void Mesh::LocalCoefficients (int triangle_index, double (&a)[3], double (&b)[3],
    double (&c)[3])
{
    int* nodes_ind = triangles[triangle_index].nodes;
    for (int i = 0; i < 3; ++i) {
        a[i] = nodes[nodes_ind[(i + 1) % 3]].x * nodes[nodes_ind[(i + 2) % 3]].y
            - nodes[nodes_ind[(i + 2) % 3]].x * nodes[nodes_ind[(i + 1) % 3]].y;
        b[i] = nodes[nodes_ind[(i + 1) % 3]].y - nodes[nodes_ind[(i + 2) % 3]].y;
        c[i] = nodes[nodes_ind[(i + 2) % 3]].x - nodes[nodes_ind[(i + 1) % 3]].x;
    }
}

bool point_in_triangle (double x, double y, int triangle_index)
{
    // TODO
}

int Mesh::TriangleForPoint (double x, double y)
{
    // TODO: hash
    for (int i = 0; i < triangles_num; ++i) {
        if (point_in_triangle(x, y, i))
            return i;
    }
    return -1;
}
