//
// Created by nickl on 3/13/19.
//

#include "neighbor.hpp"
#include "util.hpp"

namespace MeshSimpl {

namespace Internal {

void Neighbor::rotate(const E& edges, const F2E& face2edge) {
    const auto& curr_edge = edges[face2edge[face][vi]];
    const idx f_idx_to_edge = f_in_edge(curr_edge, face);
    const idx f = curr_edge.faces[1 - f_idx_to_edge];
    assert(face != f);
    face = f;
    vj = curr_edge.idx_in_face[1 - f_idx_to_edge];
    vi = get_i_from_j(vj);
}

} // namespace Internal

} // namespace MeshSimpl
