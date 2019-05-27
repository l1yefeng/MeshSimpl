//
// Created by nickl on 5/11/19.
//

#include "edge.hpp"
#include "util.hpp"
#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

Vertices* Edge::vertices = nullptr;

void Edge::planCollapse(bool fixBoundary) {
  _q = vertices->q(_vv[0]) + vertices->q(_vv[1]);

  if (fixBoundary && bothEndsOnBoundary()) {
    // the plan is: no plan is needed because it will never by modified
    return;
  } else if (fixBoundary && oneEndOnBoundary()) {
    // the plan is: new position is the position of the vertex on border
    _center = vertices->position(vertices->isBoundary(_vv[0]) ? _vv[0] : _vv[1]);
    _error = qError(_q, _center);
  } else {
    // the plan is: new position leads to the lowest error
    const vec3d b{_q[6], _q[7], _q[8]};
    const double c = _q[9];

    // computes the inverse of matrix A in quadric
    const double aDet = _q[0] * (_q[3] * _q[5] - _q[4] * _q[4]) -
                         _q[1] * (_q[1] * _q[5] - _q[4] * _q[2]) +
                         _q[2] * (_q[1] * _q[4] - _q[3] * _q[2]);

    if (aDet != 0) {
      // invertible, find position yielding minimal error
      const double aDetInv = 1.0 / aDet;
      const std::array<double, 6> aInv{
          (_q[3] * _q[5] - _q[4] * _q[4]) * aDetInv,
          (_q[2] * _q[4] - _q[1] * _q[5]) * aDetInv,
          (_q[1] * _q[4] - _q[2] * _q[3]) * aDetInv,
          (_q[0] * _q[5] - _q[2] * _q[2]) * aDetInv,
          (_q[1] * _q[2] - _q[0] * _q[4]) * aDetInv,
          (_q[0] * _q[3] - _q[1] * _q[1]) * aDetInv,
      };
      _center = {-dot({aInv[0], aInv[1], aInv[2]}, b),
                -dot({aInv[1], aInv[3], aInv[4]}, b),
                -dot({aInv[2], aInv[4], aInv[5]}, b)};
      _error = dot(b, _center) + c;
    } else {
      // not invertible, choose from endpoints and midpoint
      _center = midpoint(vertices->position(_vv[0]), vertices->position(_vv[1]));
      _error = qError(_q, _center);
      for (const idx v : _vv) {
        const double err = qError(_q, vertices->position(v));
        if (err < _error) {
          _center = vertices->position(v);
          _error = err;
        }
      }
    }
  }
}

void Edge::replaceEndpoint(idx prevV, idx newV) {
  order ord = endpointOrder(prevV);
  _vv[ord] = newV;

  assert(_vv[0] != _vv[1]);
  if (_vv[0] > _vv[1]) {
    std::swap(_vv[0], _vv[1]);
  }
}

void Edge::dropWing(idx face) {
  assert(!onBoundary());
  order ord = wingOrder(face);
  setWing(ord, 0, INVALID);
  if (ord == 0) swapWings();

  assert(ordInF(0) != INVALID);
}

void Edge::replaceWing(idx prevF, idx newF, order newOrdInFace) {
  order ord = wingOrder(prevF);
  assert(ordInF(ord) != INVALID);

  setWing(ord, newF, newOrdInFace);
}

}  // namespace Internal
}  // namespace MeshSimpl
