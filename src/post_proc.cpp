//
// Created by nickl on 1/8/19.
//

#include "post_proc.hpp"

namespace MeshSimpl {
namespace Internal {

void compact_data(V& vertices, F& indices, const Marker& marker) {
  std::vector<std::array<std::vector<idx>, 3>> vertex2face(vertices.size());

  // get rid of all deleted faces
  for (size_t lo = 0, hi = indices.size() - 1;; ++lo, --hi) {
    while (marker.exist_f(lo) && lo <= hi) ++lo;
    while (!marker.exist_f(hi) && lo < hi) --hi;
    if (lo >= hi) {
      indices.resize(lo);
      break;
    }
    std::swap(indices[lo], indices[hi]);
  }

  // create mapping from v to f
  for (idx f = 0; f < indices.size(); ++f)
    for (idx i = 0; i < 3; ++i)
      vertex2face[indices[f][i]][i].push_back(static_cast<idx>(f));

  // get rid of deleted vertices and keep the mapping valid
  for (size_t lo = 0, hi = vertices.size() - 1;; ++lo, --hi) {
    while (marker.exist_v(lo) && lo <= hi) ++lo;
    while (!marker.exist_v(hi) && lo < hi) --hi;
    if (lo >= hi) {
      vertices.resize(lo);
      vertex2face.resize(lo);
      break;
    }
    std::swap(vertex2face[lo], vertex2face[hi]);
    std::swap(vertices[lo], vertices[hi]);
  }

  for (idx v = 0; v < vertex2face.size(); ++v)
    for (idx i = 0; i < 3; ++i)
      for (const idx f : vertex2face[v][i]) indices[f][i] = static_cast<idx>(v);
}

}  // namespace Internal
}  // namespace MeshSimpl
