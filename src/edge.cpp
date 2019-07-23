//
// Created by nickl on 5/11/19.
//

#include <tuple>

#include "edge.hpp"
#include "util.hpp"
#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

bool Edge::planCollapse(bool fixBoundary) {
  _q = vertices.q(_vv[0]) + vertices.q(_vv[1]);

  if (fixBoundary) {
    std::array<bool, 2> onBounds{vertices.isBoundary(_vv[0]),
                                 vertices.isBoundary(_vv[1])};
    if (onBounds[0] && onBounds[1]) {
      // the plan is: no plan is needed because it will never by modified
      return false;
    } else if (onBounds[0] != onBounds[1]) {
      // the plan is: new position is the position of the vertex on border
      _center = vertices.position(_vv[onBounds[0] ? 0 : 1]);
      _error = _q.error(_center);
      return true;
    }
  }

  // the plan is: new position leads to the lowest error

  // computes the inverse of matrix A in quadric
  const double aDet = _q.aDeterminant();

  if (aDet != 0) {
    // invertible, find position yielding minimal error
    std::tie(_center, _error) = _q.optimal(aDet);

    // prevent the optimal position from being too far. it is anticipated that
    // such thing happens rarely, when there are coincide faces and the optimal
    // value position might be galaxy away even though any position on their
    // plane will have small enough error. error is kept as it is but we change
    // the collapse center to one of the endpoint so it looks more natural
    vec3d edgeVec = vertices.position(_vv[1]) - vertices.position(_vv[0]);
    vec3d diffVec = _center - vertices.position(_vv[0]);
    for (int i = 0; i < 3; ++i) {
      double lambda = diffVec[i] / edgeVec[i] - 0.5;
      if (std::abs(lambda) > 20) {
        if (lambda > 0) {
          _center = vertices.position(_vv[1]);
        } else {
          _center = vertices.position(_vv[0]);
        }
        break;
      }
    }

    return true;
  }

  // not invertible, choose from endpoints and midpoint
  _center = midpoint(vertices.position(_vv[0]), vertices.position(_vv[1]));
  _error = _q.error(_center);
  for (const idx v : _vv) {
    const double err = _q.error(vertices.position(v));
    if (err < _error) {
      _center = vertices.position(v);
      _error = err;
    }
  }

  return true;
}

void Edge::replaceEndpoint(idx prevV, idx newV) {
  order ord = endpointOrder(prevV);
  _vv[ord] = newV;

  assert(_vv[0] != _vv[1]);
  if (_vv[0] > _vv[1]) {
    std::swap(_vv[0], _vv[1]);
  }
}

bool Edge::dropWing(idx face) {
  order ord = wingOrder(face);
  setWing(ord, 0, INVALID);
  if (ord == 0) swapWings();

  return ordInF(0) != INVALID;
}

}  // namespace Internal
}  // namespace MeshSimpl
