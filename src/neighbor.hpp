//
// Created by nickl on 3/13/19.
//

#ifndef MESH_SIMPL_NEIGHBOR_HPP
#define MESH_SIMPL_NEIGHBOR_HPP

#include "edge.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

// This is a compact structure to represent an incident face around a center
// vertex. It is created for the traverse of the incident faces of a vertex
// (center) during edge-collapse operation.
//
// It contains a face index and two vertex (the third is the center) order: i
// and j. i, j, and center are local orders of three vertices. In typical
// (clockwise) cases, the orientation is center -> i -> j -> center
//
// Mind that in a face, any edge and the vertex across from it have the same
// order
class Neighbor {
 private:
  idx _f;  // face index in indices
  bool _ccw;  // determines the direction of rotation
  order _i;  // i = (j + 1) % 3 if clockwise
  order _j;  // j = (center + 1) % 3 if clockwise

  order iFromJ(order j) { return _ccw ? prev(j) : next(j); }

 public:
  Neighbor(idx face, order j, bool ccw)
      : _f(face), _ccw(ccw), _i(iFromJ(j)), _j(j) {}

  idx f() const { return _f; }

  order i() const { return _i; }

  order j() const { return _j; }

  order center() const { return 3 - _i - _j; }

  void rotate(const Faces& faces) {
    const Edge* currEdge = secondEdge(faces);
    assert(!currEdge->onBoundary());

    const idx prevCenter = faces[_f][center()];

    const order ford = currEdge->wingOrder(_f);
    const idx nextFace = currEdge->face(1 - ford);
    assert(_f != nextFace);
    _f = nextFace;
    _j = currEdge->ordInF(1 - ford);
    _i = iFromJ(_j);

    assert(faces[_f][center()] == prevCenter);
  }

  Edge* firstEdge(const Faces& faces) const { return faces.side(f(), _j); }

  Edge* secondEdge(const Faces& faces) const { return faces.side(f(), _i); }

  idx firstV(const Faces& faces) const { return faces[f()][_i]; }

  idx secondV(const Faces& faces) const { return faces[f()][_j]; }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_NEIGHBOR_HPP
