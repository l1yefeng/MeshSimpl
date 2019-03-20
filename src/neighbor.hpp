//
// Created by nickl on 3/13/19.
//

#ifndef LIB_MESH_SIMPL_NEIGHBOR_HPP
#define LIB_MESH_SIMPL_NEIGHBOR_HPP

#include "edge.hpp"
#include "types.hpp"

namespace MeshSimpl {

namespace Internal {

// This is a compact structure to represent an incident face around a center vertex.
// It is created for the traverse of the incident faces of a vertex (center)
// during edge-collapse operation.
//
// It contains a face index and two vertex (the third is the center) order: i and j.
// i, j, and center are local orders of three vertices.
// In typical (clockwise) cases, the orientation is center -> i -> j -> center
//
// Mind that in a face, any edge and the vertex across from it have the same order
class Neighbor {
private:
    idx face; // face index in indices
    bool ccw; // set to false when constructed
    order vi; // i = (j + 1) % 3 if clockwise
    order vj; // j = (center + 1) % 3 if clockwise

    order get_i_from_j(order j) {
        if (!ccw)
            return (j + 1) % 3;
        else
            return (j + 2) % 3;
    }

public:
    Neighbor(idx face, order vj, bool ccw)
        : face(face), ccw(ccw), vi(get_i_from_j(vj)), vj(vj) {}
    idx f() const { return face; }
    order i() const { return vi; }
    order j() const { return vj; }
    order center() const { return 3 - vi - vj; }

    void rotate(const E& edges, const F2E& face2edge) {
        const auto& curr_edge = edges[face2edge[face][vi]];
        assert(!curr_edge.on_boundary());

        const order ford = curr_edge.f_order(face);
        const idx next_face = curr_edge.faces[1 - ford];
        assert(face != next_face);
        face = next_face;
        vj = curr_edge.idx_in_face[1 - ford];
        vi = get_i_from_j(vj);
    }
};

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_NEIGHBOR_HPP
