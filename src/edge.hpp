//
// Created by nickl on 03/20/19.
//

#ifndef LIB_MESH_SIMPL_EDGE_HPP
#define LIB_MESH_SIMPL_EDGE_HPP

#include "types.hpp"
#include <array>
#include <cassert>

namespace MeshSimpl {
namespace Internal {

// Edge defines the struct of an edge
//
// Only some convenient methods are defined for the purpose of less repetition.
// Methods to initialize (adding first and second face), and update are not provided.
struct Edge {
    static const order INVALID = -1;
    enum BOUNDARY_V { V0 = 0, V1 = 1, BOTH, NEITHER };

    // sum of quadrics of two endpoints
    Quadric q;
    // where this edge collapse into
    vec3d center;

    // index of two endpoints
    vec2i vertices;

    // index of two incident faces;
    // boundary edges have faces[0] and ord_in_face[1] is INVALID;
    // all edges have at least face[0] after edges are constructed in `construct_edges()`
    vec2i faces;
    // order of this edge in two incident faces, match the order to member "faces"
    std::array<order, 2> ord_in_faces;

    // quadric error value
    double error;

    // vertices on boundary
    BOUNDARY_V boundary_v;

    //
    // public methods for convenient retrieval of information
    //

    bool on_boundary() const { return ord_in_faces[1] == INVALID; }

    order v_order(idx v) const {
        assert(vertices[0] == v || vertices[1] == v);
        return vertices[0] == v ? 0 : 1;
    }

    order f_order(idx f) const {
        assert(faces[0] == f || (faces[1] == f && ord_in_faces[1] != INVALID));
        return faces[0] == f ? 0 : 1;
    }

    order v_del_order() const { return boundary_v == V0 ? 1 : 0; }

    //
    // public methods for convenient update
    //

    void swap_faces() {
        std::swap(faces[0], faces[1]);
        std::swap(ord_in_faces[0], ord_in_faces[1]);
    }
};

} // namespace Internal
} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_EDGE_HPP
