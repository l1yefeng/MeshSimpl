//
// Created by nickl on 3/21/19.
//

#ifndef LIB_MESH_SIMPL_CONNECTIVITY_HPP
#define LIB_MESH_SIMPL_CONNECTIVITY_HPP

#include "types.hpp"
#include <cassert>

namespace MeshSimpl {
namespace Internal {

struct Connectivity {
    F& indices;
    E edges;
    F2E face2edge;

    Edge& edge_of_face(idx f, order k) { return edges[face2edge[f][k]]; }
    Edge const& edge_of_face(idx f, order k) const { return edges[face2edge[f][k]]; }

    order v_ord_in_face(idx f, idx v) const {
        assert(indices[f][0] == v || indices[f][1] == v || indices[f][2] == v);
        return indices[f][0] == v ? 0 : (indices[f][1] == v ? 1 : 2);
    }

    idx edge_idx_across_from_v(idx f, idx v) const {
        return face2edge[f][v_ord_in_face(f, v)];
    }
};

} // namespace Internal
} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_CONNECTIVITY_HPP
