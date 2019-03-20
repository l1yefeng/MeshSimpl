//
// Created by nickl on 03/20/19.
//

#ifndef LIB_MESH_SIMPL_EDGE_H
#define LIB_MESH_SIMPL_EDGE_H

#include "types.hpp"
#include <array>
#include <cassert>

namespace MeshSimpl {

namespace Internal {

// Edge defines the struct of an edge
struct Edge {
    static const order INVALID = -1;
    enum BOUNDARY_V { V0 = 0, V1 = 1, BOTH, NONE };

    // sum of quadrics of two endpoints
    Quadric q;
    // where this edge collapse into
    vec3d center;

    // index of two endpoints,
    // when constructing a set of edges, this is used as key by forcing v0 < v1,
    // the increasing order is unnecessary if not in a set, but mind "boundary_v"
    vec2i vertices;

    // index of two incident faces; boundary edges have faces[0]
    vec2i faces;
    // order of this edge in two incident faces, match the order to member "faces"
    std::array<order, 2> idx_in_face;

    // quadric error value
    double error;

    // vertices on boundary
    BOUNDARY_V boundary_v;

    // public methods for convenient retrieval of information
    bool on_boundary() const { return idx_in_face[1] == INVALID; }

    order v_order(idx v) const {
        assert(vertices[0] == v || vertices[1] == v);
        return vertices[0] == v ? 0 : 1;
    }

    order f_order(idx f) const {
        assert(faces[0] == f || (faces[1] == f && idx_in_face[1] != INVALID));
        return faces[0] == f ? 0 : 1;
    }

    order v_del_order() const { return boundary_v == V0 ? 1 : 0; }

    // comparison
    bool operator<(const Edge& rhs) const {
        if (vertices[0] < rhs.vertices[0])
            return true;
        if (vertices[0] > rhs.vertices[0])
            return false;
        return vertices[1] < rhs.vertices[1];
    }
};

typedef std::vector<Edge> E;
typedef std::vector<Quadric> Q;
typedef std::vector<vec3i> F2E;

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_EDGE_H
