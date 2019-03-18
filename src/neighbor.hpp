//
// Created by nickl on 3/13/19.
//

#ifndef LIB_MESH_SIMPL_NEIGHBOR_H
#define LIB_MESH_SIMPL_NEIGHBOR_H

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
    Neighbor(idx face, order vj, bool ccw = false) : face(face), ccw(ccw), vi(get_i_from_j(vj)), vj(vj) {}
    idx f() const { return face; }
    order i() const { return vi; }
    order j() const { return vj; }
    order center() const { return 3 - vi - vj; }

    // Start traversing in counter clockwise, used when encounter special case (boundary)
    void counter_clockwise() { ccw = true; }

    void rotate(const E& edges, const F2E& face2edge);
};

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_NEIGHBOR_H