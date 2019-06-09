//
// Created by nickl on 6/9/19.
//

#include <stddef.h>
#include <array>
#include <initializer_list>

#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

void Vertices::compactPositionsAndDie(Positions& positions, Indices& indices) {
  std::vector<std::array<std::vector<idx>, 3>> vertex2face(size());

  // create mapping from v to f
  for (idx f = 0; f < indices.size(); ++f) {
    for (idx i = 0; i < 3; ++i) vertex2face[indices[f][i]][i].push_back(f);
  }

  // get rid of deleted vertices and keep the mapping valid
  for (size_t lo = 0, hi = size() - 1; true; ++lo, --hi) {
    while (lo <= hi && exists(lo)) ++lo;
    while (lo < hi && !exists(hi)) --hi;
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

}  // namespace Internal
}  // namespace MeshSimpl
