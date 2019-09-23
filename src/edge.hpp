//
// Created by nickl on 03/20/19.
//

#ifndef MESH_SIMPL_EDGE_HPP
#define MESH_SIMPL_EDGE_HPP

#include <array>
#include <cassert>
#include <limits>
#include <utility>

#include "erasable.hpp"
#include "quadric.hpp"
#include "types.hpp"
#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

class Edge : public Erasable {
 private:
  Vertices &vertices;

  // sum of quadrics of two endpoints
  Quadric _q = {};

  // where this edge collapse into
  vec3d _center = {};

  // quadric error value
  double _error = 0;

  // index of two endpoints
  vec2i _vv;

  // wings
  vec2i _ff = {};
  std::array<order, 2> _ordInFF = {{INVALID, INVALID}};

  void swapWings() {
    std::swap(_ff[0], _ff[1]);
    std::swap(_ordInFF[0], _ordInFF[1]);
  }

 public:
  // Construct an edge with two endpoints indexes.
  // After constructed this edge, call setWing() once or twice
  // Then this edge is correctly initialized
  Edge(Vertices &vertices, idx v0, idx v1)
      : Erasable(), vertices(vertices), _vv({v0, v1}) {}

  // Attach a face to edge.
  // No other modifier member functions should be called before this call
  void setWing(order wingOrd, idx f, order ordInF) {
    _ff[wingOrd] = f;
    _ordInFF[wingOrd] = ordInF;
  }

  //
  // public methods for retrieval of information
  //

  // Returns the collapse center of next collapse
  const vec3d &center() const { return _center; }

  // Returns the error value made of next collapse
  double error() const { return _error; }

  // Returns the sum of quadrics of two endpoints
  const Quadric &q() const { return _q; }

  bool bothEndsOnBoundary() const {
    return vertices.isBoundary(_vv[0]) && vertices.isBoundary(_vv[1]);
  }

  bool neitherEndOnBoundary() const {
    return !vertices.isBoundary(_vv[0]) && !vertices.isBoundary(_vv[1]);
  }

  idx face(order ord) const { return _ff[ord]; }

  // Returns the order of this edge in face(ord)
  order ordInF(order ord) const { return _ordInFF[ord]; }

  // Returns true if edge is on border
  // NOTE: different from both endpoints on border
  bool onBoundary() const { return ordInF(1) == INVALID; }

  // Returns the order of face f on this edge
  order wingOrder(idx f) const {
    assert(face(0) == f || (face(1) == f && ordInF(1) != INVALID));
    return face(0) == f ? 0 : 1;
  }

  // Returns the order of endpoint v on this edge
  order endpointOrder(idx v) const {
    assert(_vv[0] == v || _vv[1] == v);
    return _vv[0] == v ? 0 : 1;
  }

  const vec2i &endpoints() const { return _vv; }

  idx endpoint(order ord) const { return endpoints()[ord]; }

  //
  // public methods for update
  //

  // Plane next collapse.
  // Will set
  //  - what is current sum of endpoints quadrics
  //  - which position to collapse into (center)
  //  - what will be the error
  bool planCollapse(bool fixBoundary);

  void setErrorInfty() { _error = std::numeric_limits<double>::max(); }

  void replaceEndpoint(idx prevV, idx newV);

  // Drop one wing: if it has two wings before then now it has one; if it
  // was a boundary edge (has one wing) it will have no valid wing thus an
  // invalid edge.
  // Returns true if still has one wing; false if invalid
  // User should probably erase this edge after receiving a false return value
  bool dropWing(idx face);
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_EDGE_HPP
