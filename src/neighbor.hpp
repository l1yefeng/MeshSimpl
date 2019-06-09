//
// Created by nickl on 3/13/19.
//

#ifndef MESH_SIMPL_NEIGHBOR_HPP
#define MESH_SIMPL_NEIGHBOR_HPP

#include "edge.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

// Neighbor represents an incident face around a center
// vertex. It is created for the traverse of the incident faces of a vertex
// (center) during edge-collapse operation.
//
// Mind that in a face, any edge and the vertex across from it have the same
// order. Therefore there is no firstEdgeOrd(), use secondVOrd().
class Neighbor {
 private:
  const Faces& faces;
  idx _f;         // face index in indices
  bool _ccw;      // determines the direction of rotation
  order _second;  // = (first + 1) % 3 if counter-clockwise
  order _first;   // = (second - 1) % 3 = (center + 1) % 3 if counter-clockwise

  order getFirst(order second) { return _ccw ? prev(second) : next(second); }

 public:
  Neighbor(const Edge* firstEdge, order toWingOrder, idx vCenter,
           const Faces& faces)
      : faces(faces),
        _f(firstEdge->face(toWingOrder)),
        _ccw(firstEdge->ordInF(toWingOrder) !=
             next(faces.orderOf(firstEdge->face(toWingOrder), vCenter))),
        _second(firstEdge->ordInF(toWingOrder)),
        _first(getFirst(_second)) {
    assert(firstEdge->exists());
    assert(firstEdge->ordInF(toWingOrder) != INVALID);
  }

  void replace(const Edge* firstEdge, order toWingOrder, idx vCenter) {
    _f = firstEdge->face(toWingOrder);
    _ccw = firstEdge->ordInF(toWingOrder) !=
           next(faces.orderOf(firstEdge->face(toWingOrder), vCenter));
    _second = firstEdge->ordInF(toWingOrder);
    _first = getFirst(_second);

    assert(firstEdge->exists());
    assert(firstEdge->ordInF(toWingOrder) != INVALID);
  }

  idx f() const { return _f; }

  order firstVOrd() const { return _first; }

  order secondVOrd() const { return _second; }

  order center() const { return 3 - _first - _second; }

  void rotate() {
    const Edge* currEdge = secondEdge();
    assert(!currEdge->onBoundary());

    const idx prevCenter = faces[_f][center()];

    const order ford = currEdge->wingOrder(_f);
    const idx nextFace = currEdge->face(1 - ford);
    assert(_f != nextFace);
    _f = nextFace;
    _second = currEdge->ordInF(1 - ford);
    _first = getFirst(_second);

    assert(faces[_f][center()] == prevCenter);
  }

  Edge* firstEdge() const { return faces.side(f(), _second); }

  Edge* secondEdge() const { return faces.side(f(), _first); }

  idx firstV() const { return faces[f()][_first]; }

  idx secondV() const { return faces[f()][_second]; }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_NEIGHBOR_HPP
