#include "mesh3d.h"
#include <array>
#include <cmath>
#include <fstream>
#include <algorithm>

void Mesh::add_boundary_node (int node_index)
{
    if (boundary_indices[node_index] == -1) {
        boundary_indices[node_index] = boundary_nodes_num;
        boundary_nodes.resize(boundary_nodes_num + 1);
        boundary_nodes[boundary_nodes_num] = node_index;
        ++ boundary_nodes_num;
    }
}

void Mesh::add_adjacent_node (int node, int node_to_add)
{
    std::list<int>& l = adjacent_nodes[node];
    if (std::find(l.begin(), l.end(), node_to_add) == l.end())
        l.push_back(node_to_add);
}

Mesh::Mesh (std::string file_name)
{
    std::ifstream fin(file_name);
    if (!fin.is_open())
        throw "Error opening the file " + file_name;

    std::string dummy_str;
    int dummy_int;
    for (int i = 0; i < 3; ++i)
        std::getline(fin, dummy_str);

    // Vertices
    fin >> dummy_str;
    fin >> nodes_num;
    nodes.resize(nodes_num);
    for (int i = 0; i < nodes_num; ++i)
        fin >> nodes[i].x >> nodes[i].y >> nodes[i].z >> dummy_int;

    // Tetrahedra
    fin >> dummy_str;
    fin >> tetrahedrons_num;
    tetrahedrons.resize(tetrahedrons_num);
    for (int i = 0; i < tetrahedrons_num; ++i) {
        for (int j = 0; j < 4; ++j) {
            fin >> tetrahedrons[i].nodes[j];
            -- tetrahedrons[i].nodes[j];
        }
        fin >> dummy_int;
    }
    tetrahedrons_for_nodes.resize(nodes_num);
    for (int i = 0; i < tetrahedrons_num; ++i) {
        for (int j = 0; j < 4; ++j)
            tetrahedrons_for_nodes[tetrahedrons[i].nodes[j]].push_back(i);
    }
    adjacent_nodes.resize(nodes_num);
    for (int i = 0; i < nodes_num; ++i)
        adjacent_nodes[i].push_back(i);
    for (int i = 0; i < tetrahedrons_num; ++i) {
        int* nodes = tetrahedrons[i].nodes;
        for (int j = 0; j < 4; ++j) {
            for (int k = j + 1; k < 4; ++k) {
                add_adjacent_node(nodes[j], nodes[k]);
                add_adjacent_node(nodes[k], nodes[j]);
            }
        }
    }

    // Triangles
    fin >> dummy_str;
    fin >> boundary_triangles_num;
    std::vector<std::array<int, 3> > triangles;
    triangles.resize(boundary_triangles_num);
    for (int i = 0; i < boundary_triangles_num; ++i) {
        for (int j = 0; j < 3; ++j) {
            fin >> triangles[i][j];
            -- triangles[i][j];
        }
        fin >> dummy_int;
    }
    boundary_nodes_num = 0;
    boundary_indices.resize(nodes_num);
    for (int i = 0; i < nodes_num; ++i)
        boundary_indices[i] = -1;
    for (int i = 0; i < boundary_triangles_num; ++i) {
        for (int j = 0; j < 3; ++j)
            add_boundary_node(triangles[i][j]);
    }
    boundary_triangles.resize(boundary_triangles_num);
    triangles_for_boundary_nodes.resize(boundary_nodes_num);
    for (int i = 0; i < boundary_triangles_num; ++i) {
        for (int j = 0; j < 3; ++j) {
            int boundary_node = boundary_indices[triangles[i][j]];
            boundary_triangles[i].boundary_nodes[j] = boundary_node;
            triangles_for_boundary_nodes[boundary_node].push_back(i);
        }
    }
}

double det_3_3 (double a11, double a12, double a13,
    double a21, double a22, double a23, double a31, double a32, double a33)
{
    return a11 * a22 * a33  +  a12 * a23 * a31  +  a13 * a21 * a32
        -  a13 * a22 * a31  -  a11 * a23 * a32  -  a12 * a21 * a33;
}

double det_4_4 (double a11, double a12, double a13, double a14,
    double a21, double a22, double a23, double a24,
    double a31, double a32, double a33, double a34,
    double a41, double a42, double a43, double a44)
{
    return a11 * det_3_3(a22, a23, a24,  a32, a33, a34,  a42, a43, a44)
        - a21 * det_3_3(a12, a13, a14,  a32, a33, a34,  a42, a43, a44)
        + a31 * det_3_3(a12, a13, a14,  a22, a23, a24,  a42, a43, a44)
        - a41 * det_3_3(a12, a13, a14,  a22, a23, a24,  a32, a33, a34);
}

double Mesh::SignedTetrahedronVolume (int tetrahedron_index)
{
    int* nodes_ind = tetrahedrons[tetrahedron_index].nodes;
    return (1. / 6.) * det_4_4(
        1, nodes[nodes_ind[0]].x, nodes[nodes_ind[0]].y, nodes[nodes_ind[0]].z,
        1, nodes[nodes_ind[1]].x, nodes[nodes_ind[1]].y, nodes[nodes_ind[1]].z,
        1, nodes[nodes_ind[2]].x, nodes[nodes_ind[2]].y, nodes[nodes_ind[2]].z,
        1, nodes[nodes_ind[3]].x, nodes[nodes_ind[3]].y, nodes[nodes_ind[3]].z);
}

double Mesh::TetrahedronVolume (int tetrahedron_index)
{
    return fabs(SignedTetrahedronVolume(tetrahedron_index));
}

double Mesh::BoundaryTriangleArea (int boundary_triangle_index)
{
    int* nodes_ind =
        boundary_triangles[boundary_triangle_index].boundary_nodes;
    double x1 = nodes[boundary_nodes[nodes_ind[0]]].x;
    double y1 = nodes[boundary_nodes[nodes_ind[0]]].y;
    double z1 = nodes[boundary_nodes[nodes_ind[0]]].z;
    double x2 = nodes[boundary_nodes[nodes_ind[1]]].x;
    double y2 = nodes[boundary_nodes[nodes_ind[1]]].y;
    double z2 = nodes[boundary_nodes[nodes_ind[1]]].z;
    double x3 = nodes[boundary_nodes[nodes_ind[2]]].x;
    double y3 = nodes[boundary_nodes[nodes_ind[2]]].y;
    double z3 = nodes[boundary_nodes[nodes_ind[2]]].z;
    return 0.5 * sqrt(pow(det_3_3(1, y1, z1, 1, y2, z2, 1, y3, z3), 2)
        + pow(det_3_3(1, x1, z1, 1, x2, z2, 1, x3, z3), 2)
        + pow(det_3_3(1, x1, y1, 1, x2, y2, 1, x3, y3), 2));
}

void Mesh::LocalCoefficients (int tetrahedron_index,
        double (&a)[4], double (&b)[4], double (&c)[4], double (&d)[4])
{
    // See Zienkiewicz et al., section 4.11
    int* nodes_ind = tetrahedrons[tetrahedron_index].nodes;
    for (int i = 0; i < 4; ++i) {
        a[i] = pow(-1, i) * det_3_3(nodes[nodes_ind[(i + 1) % 4]].x,
            nodes[nodes_ind[(i + 1) % 4]].y, nodes[nodes_ind[(i + 1) % 4]].z,
            nodes[nodes_ind[(i + 2) % 4]].x,
            nodes[nodes_ind[(i + 2) % 4]].y, nodes[nodes_ind[(i + 2) % 4]].z,
            nodes[nodes_ind[(i + 3) % 4]].x,
            nodes[nodes_ind[(i + 3) % 4]].y, nodes[nodes_ind[(i + 3) % 4]].z
        );
        b[i] = - pow(-1, i) * det_3_3(1, nodes[nodes_ind[(i + 1) % 4]].y,
            nodes[nodes_ind[(i + 1) % 4]].z,
            1, nodes[nodes_ind[(i + 2) % 4]].y, nodes[nodes_ind[(i + 2) % 4]].z,
            1, nodes[nodes_ind[(i + 3) % 4]].y, nodes[nodes_ind[(i + 3) % 4]].z
        );
        c[i] = - pow(-1, i) * det_3_3(nodes[nodes_ind[(i + 1) % 4]].x, 1,
            nodes[nodes_ind[(i + 1) % 4]].z,
            nodes[nodes_ind[(i + 2) % 4]].x, 1, nodes[nodes_ind[(i + 2) % 4]].z,
            nodes[nodes_ind[(i + 3) % 4]].x, 1, nodes[nodes_ind[(i + 3) % 4]].z
        );
        d[i] = - pow(-1, i) * det_3_3(nodes[nodes_ind[(i + 1) % 4]].x,
            nodes[nodes_ind[(i + 1) % 4]].y, 1,
            nodes[nodes_ind[(i + 2) % 4]].x, nodes[nodes_ind[(i + 2) % 4]].y, 1,
            nodes[nodes_ind[(i + 3) % 4]].x, nodes[nodes_ind[(i + 3) % 4]].y, 1
        );
    }
}

void Mesh::LocalCoordinates (int tetrahedron_index, double x, double y, double z,
    double& L0, double& L1, double& L2)
{
    double A = SignedTetrahedronVolume(tetrahedron_index);
    double a[4], b[4], c[4], d[4];
    LocalCoefficients(tetrahedron_index, a, b, c, d);
    L0 = 1 / (6 * A) * (a[0] + b[0] * x + c[0] * y + d[0] * z);
    L1 = 1 / (6 * A) * (a[1] + b[1] * x + c[1] * y + d[1] * z);
    L2 = 1 / (6 * A) * (a[2] + b[2] * x + c[2] * y + d[2] * z);
}

bool Mesh::PointInTetrahedron (int tetrahedron_index,
    double x, double y, double z)
{
    const double eps = 1e-5;
    double L[4];
    LocalCoordinates(tetrahedron_index, x, y, z, L[0], L[1], L[2]);
    L[3] = 1 - L[0] - L[1] - L[2];
    bool ans = true;
    for (int i = 0; i < 4; ++i)
        ans = ans && (0 - eps < L[i]) && (L[i] < 1 + eps);
    return ans;
}

int Mesh::TetrahedronForPoint (double x, double y, double z)
{
    // TODO: hash
    for (int i = 0; i < tetrahedrons_num; ++i) {
        if (PointInTetrahedron(i, x, y, z))
            return i;
    }
    return -1;
}
