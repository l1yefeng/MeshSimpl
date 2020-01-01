//
// Created by nickl on 5/25/19.
//

#ifndef MESH_SIMPL_VERTICES_HPP
#define MESH_SIMPL_VERTICES_HPP

#include <ext/alloc_traits.h>
#include <utility>
#include <vector>

#include "erasable.hpp"
#include "faces.hpp"
#include "quadric.hpp"
#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

class Vertices : public Erasables {
 private:
  Positions _positions;
  std::vector<Quadric> _quadrics;
  std::vector<bool> _boundary;
  std::vector<bool> _fixed;

 public:
  // Embed positions and allocate space for quadrics
  explicit Vertices(Positions& positions)
      : Erasables(positions.size()),
        _positions(std::move(positions)),
        _quadrics(size()),
        _boundary(size(), false),
        _fixed(size(), false) {}

  // Get/set position of a vertex
  const vec3d& position(idx v) const { return _positions[v]; }
  const vec3d& operator[](idx v) const { return position(v); }
  void setPosition(idx v, const vec3d& pos) { _positions[v] = pos; }

  // Get/set quadric of a vertex
  const Quadric& q(idx v) const { return _quadrics[v]; }
  void increaseQ(idx v, const Quadric& by) { setQ(v, q(v) + by); }
  void setQ(idx v, const Quadric& val) { _quadrics[v] = val; }

  // Get/set if a vertex is on boundary
  bool isBoundary(idx v) const { return _boundary[v]; }
  void setBoundary(idx v, bool b) { _boundary[v] = b; }

  // Get/set if a vertex is fixed
  bool isFixed(idx v) const { return _fixed[v]; }
  void setFixed(idx v, bool b) { _fixed[v] = b; }

  // Erase unreferenced vertices
  void eraseUnref(const Faces& faces) {
    _erased = std::vector<bool>(size(), true);
    for (idx f = 0; f < faces.size(); ++f) {
      if (faces.exists(f)) {
        for (idx v : faces[f]) {
          _erased[v] = false;
        }
      }
    }
  }

  void reduceQByHalf(idx v) { _quadrics[v] *= 0.5; }

  idx duplicate(idx src) {
    idx v = size();
    _erased.push_back(false);
    _positions.push_back(_positions[src]);
    _quadrics.push_back(_quadrics[src]);
    _boundary.push_back(_boundary[src]);
    _fixed.push_back(_fixed[src]);
    return v;
  }

  void compactPositionsAndDie(Positions& positions, Indices& indices);
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_VERTICES_HPP
