//
// Created by nickl on 1/8/19.
//

#include "simplify.hpp"

#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>

#include "collapser.hpp"
#include "edge.hpp"
#include "faces.hpp"
#include "proc.hpp"
#include "qemheap.hpp"
#include "vertices.hpp"

namespace MeshSimpl {

using namespace Internal;

namespace Internal {
void validateOptions(const SimplifyOptions &options,
                     const Positions &positions) {
  // clang-format off
  if (options.strength > 1)
    throw std::invalid_argument("ERROR::INVALID_OPTION: strength > 1");
  if (options.strength < 0)
    throw std::invalid_argument("ERROR::INVALID_OPTION: strength < 0");
  if (options.borderConstraint < 0)
    throw std::invalid_argument("ERROR::INVALID_OPTION: border constraint < 0");
  if (options.foldOverAngleThreshold > 1 || options.foldOverAngleThreshold < -1)
    throw std::invalid_argument("ERROR::INVALID_OPTION: fold-over angle not between -1 and 1");
  if (options.aspectRatioThreshold > 1.0)
    throw std::invalid_argument("ERROR::INVALID_OPTION: aspect-ratio-threshold cannot exceed 1");
  if (!options.fixedVertices.empty() && options.fixedVertices.size() != positions.size())
    throw std::invalid_argument("ERROR::INVALID_OPTION: fixedVertices is neither empty nor equal with 'positions' in size");
  // clang-format on
}
}  // namespace Internal

void simplify(Positions &positions, Indices &indices,
              const SimplifyOptions &options) {
  validateOptions(options, positions);

  const size_t NF = indices.size();
  const size_t nfToDecimate = std::lround(options.strength * NF);

  if (nfToDecimate == 0) return;

  // construct vertices and faces from positions and indices
  // positions and indices are moved and no longer hold data
  Vertices vertices(positions);
  Faces faces(indices);

  // find out information of edges (endpoints, incident faces) and face2edge
  Edges edges;
  buildConnectivity(vertices, faces, edges);

  // determine each vertex should be fixed or not
  if (options.fixedVertices.empty()) {
    if (options.fixBoundary)
      for (idx v = 0; v < vertices.size(); ++v)
        vertices.setFixed(v, vertices.isBoundary(v));
  } else {
    for (idx v = 0; v < vertices.size(); ++v)
      vertices.setFixed(v, options.fixedVertices[v]);
  }

  // compute quadrics of vertices
  computeQuadrics(vertices, faces, options);

  // assigning edge errors using quadrics
  QEMHeap heap(edges);
  for (size_t e = 0; e < edges.size(); ++e) {
    if (!edges[e].planCollapse()) {
      heap.markRemovedById(e);
    }
  }
  heap.prioritize();

  int nf = nfToDecimate;
  Collapser collapser(vertices, faces, heap, options);
  while (!heap.empty() && nf > 0) {
    // target the least-error edge, if it is what we saw last iteration,
    // it means loop should stop because all remaining edges have been penalized
    Edge *const edge = heap.top();
    if (!edge->exists() || !heap.contains(edge)) {
      heap.pop();
      continue;
    }
    if (edge->error() >= std::numeric_limits<double>::max()) break;

    // collapse the least-error edge until mesh is simplified enough
    int removed = collapser.collapse(edge);
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
