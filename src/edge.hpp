//
// Created by nickl on 03/20/19.
//

#ifndef MESH_SIMPL_EDGE_HPP
#define MESH_SIMPL_EDGE_HPP

#include <array>      // for array
#include <cassert>    // for assert
#include <limits>     // for numeric_limits
#include <utility>    // for swap
#include "types.hpp"  // for order, vec2i, idx, Quadric, vec3d
#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

// Edge defines the struct of an edge
//
// Only some convenient methods are defined for the purpose of less repetition.
// Methods to initialize (adding first and second face), and update are not
// provided.
class Edge {
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

 public:
  // Construct an edge with two endpoints indexes.
  // After constructed this edge, call setWing() once or twice
  // Then this edge is correctly initialized
  Edge(Vertices &vertices, idx v0, idx v1)
      : vertices(vertices), _vv({v0, v1}) {}

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

  bool oneEndOnBoundary() const {
    return !bothEndsOnBoundary() && !neitherEndOnBoundary();
  }

  idx face(order ord) const { return _ff[ord]; }

  // Returns the order of this edge in face(ord)
  order ordInF(order ord) const { return _ordInFF[ord]; }

  // Returns true if egde is on border
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

  // Returns the order of the endpoint next collapse should delete.
  // By keeping one endpoint and deleting the other, followed by update on the
  // kept endpoint (position, quadric, etc), the edge is collapsed into a new
  // vertex (the center)
  order delEndpointOrder() const { return vertices.isBoundary(_vv[0]) ? 1 : 0; }

  const vec2i &endpoints() const { return _vv; }

  idx operator[](order ord) const { return _vv[ord]; }

  //
  // public methods for update
  //

  // Plane next collapse.
  // Will set
  //  - what is current sum of endpoints quadrics
  //  - which position to collapse into (center)
  //  - what will be the error
  void planCollapse(bool fixBoundary);

  void setErrorInfty() { _error = std::numeric_limits<double>::max(); }

  void swapWings() {
    std::swap(_ff[0], _ff[1]);
    std::swap(_ordInFF[0], _ordInFF[1]);
  }

  void replaceEndpoint(idx prevV, idx newV);

  void dropWing(idx face);

  void replaceWing(idx prevF, idx newF, order newOrdInFace);
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_EDGE_HPP
