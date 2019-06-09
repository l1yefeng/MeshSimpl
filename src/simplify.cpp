//
// Created by nickl on 1/8/19.
//

#include <cstddef>
#include <limits>
#include <stdexcept>

#include "collapser.hpp"
#include "edge.hpp"
#include "faces.hpp"
#include "proc.hpp"
#include "qemheap.hpp"
#include "simplify.hpp"
#include "vertices.hpp"

namespace MeshSimpl {

using namespace Internal;

void validateOptions(const SimplifyOptions &options) {
  if (options.strength > 1)
    throw std::invalid_argument("ERROR::INVALID_OPTION: strength > 1");
  if (options.strength < 0)
    throw std::invalid_argument("ERROR::INVALID_OPTION: strength < 0");
}

void simplify(Positions &positions, Indices &indices,
              const SimplifyOptions &options) {
  validateOptions(options);

  const size_t NF = indices.size();
  const size_t nfToDecimate = std::lround(options.strength * NF);

  if (nfToDecimate == 0) return;

  // construct vertices and faces from positions and indices
  // positions and indices are moved and no longer hold data
  Vertices vertices(positions);
  Faces faces(indices);

  // [1] find out information of edges (endpoints, incident faces) and face2edge
  Edges edges;
  buildConnectivity(vertices, faces, edges);

  // [2] compute quadrics of vertices
  computeQuadrics(vertices, faces, edges, options);

  // [3] assigning edge errors using quadrics
  for (auto &edge : edges) edge.planCollapse(options.fixBoundary);

  // [4] create priority queue on quadric error
  QEMHeap heap(edges);
  if (options.fixBoundary) {
    for (auto &edge : edges) {
      if (edge.bothEndsOnBoundary()) heap.remove(&edge);
    }
  }

  int nf = nfToDecimate;
  while (!heap.empty() && nf > 0) {
    // target the least-error edge, if it is what we saw last iteration,
    // it means loop should stop because all remaining edges have been penalized
    Edge *const edge = heap.top();
    if (!edge->exists() || !heap.contains(edge)) {
      heap.pop();
      continue;
    }
    if (edge->error() >= std::numeric_limits<double>::max()) break;

    // [5] collapse the least-error edge until mesh is simplified enough
    Collapser collapser(vertices, faces, heap, edge, options);
    int removed = collapser.collapse();
    nf -= removed;
  }

  vertices.eraseUnref(faces);

  // edges are useless
  // faces and vertices will be used to generate indices and positions
  // then they are useless as well
  faces.compactIndicesAndDie(indices);
  vertices.compactPositionsAndDie(positions, indices);
}

}  // namespace MeshSimpl
