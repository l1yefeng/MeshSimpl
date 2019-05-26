//
// Created by nickl on 5/25/19.
//

#ifndef MESH_SIMPL_VERTICES_HPP
#define MESH_SIMPL_VERTICES_HPP

#include "erasable.hpp"
#include "faces.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

class Vertices : public Erasables {
 private:
  Positions _positions;
  std::vector<Quadric> _quadrics;

 public:
  // Embed positions and allocate space for quadrics
  explicit Vertices(Positions& positions)
      : Erasables(positions.size()),
        _positions(std::move(positions)),
        _quadrics(size()) {}

  // Get/set position of a vertex
  const vec3d& position(idx v) const {
    assert(exists(v));
    return _positions[v];
  }
  const vec3d& operator[](idx v) const { return position(v); }
  void setPosition(idx v, const vec3d& pos) {
    assert(exists(v));
    _positions[v] = pos;
  }

  // Get/set quadric of a vertex
  const Quadric& q(idx v) const {
    assert(exists(v));
    return _quadrics[v];
  }
  void increaseQ(idx v, const Quadric& by) { setQ(v, q(v) + by); }
  void setQ(idx v, const Quadric& val) {
    assert(exists(v));
    _quadrics[v] = val;
  }

  // Erase unreferenced vertices
  void eraseUnref(const Faces& faces) {
    erased = std::vector<bool>(size(), true);
    for (idx f = 0; f < size(); ++f) {
      for (idx v : faces[f]) erased[v] = false;
    }
  }

  void compactPositionsAndDie(Positions& positions, Indices& indices) {
    std::vector<std::array<std::vector<idx>, 3>> vertex2face(size());

    // create mapping from v to f
    for (idx f = 0; f < indices.size(); ++f) {
      for (idx i = 0; i < 3; ++i) vertex2face[indices[f][i]][i].push_back(f);
    }

    // get rid of deleted vertices and keep the mapping valid
    for (size_t lo = 0, hi = size() - 1;; ++lo, --hi) {
      while (exists(lo) && lo <= hi) ++lo;
      while (!exists(hi) && lo < hi) --hi;
      if (lo >= hi) {
        _positions.resize(lo);
        vertex2face.resize(lo);
        break;
      }
      std::swap(vertex2face[lo], vertex2face[hi]);
      std::swap(_positions[lo], _positions[hi]);
    }

    // move to output
    positions = std::move(_positions);

    // update indices
    for (idx v = 0; v < vertex2face.size(); ++v) {
      for (order i : {0, 1, 2})
        for (idx f : vertex2face[v][i]) indices[f][i] = v;
    }
  }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_VERTICES_HPP
