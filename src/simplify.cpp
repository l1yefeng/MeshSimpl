//
// Created by nickl on 1/8/19.
//

#include "simplify.hpp"
#include <cstddef>       // for size_t
#include <limits>        // for numeric_limits
#include <stdexcept>     // for invalid_argument
#include "ecol.hpp"      // for edge_collapse
#include "edge.hpp"      // for Edge
#include "faces.hpp"     // for Faces
#include "pre_proc.hpp"  // for compute_quadrics, construct_edges
#include "qem_heap.hpp"  // for QEMHeap
#include "vertices.hpp"  // for Vertices

namespace MeshSimpl {

using namespace Internal;

void options_validation(const SimplifyOptions &options) {
  if (options.strength > 1)
    throw std::invalid_argument("ERROR::INVALID_OPTION: strength > 1");
  if (options.strength < 0)
    throw std::invalid_argument("ERROR::INVALID_OPTION: strength < 0");
}

void simplify(Positions &positions, Indices &indices,
              const SimplifyOptions &options) {
  options_validation(options);

  // construct vertices and faces from positions and indices
  // positions and indices are moved and no longer hold data
  Vertices vertices(positions);
  Faces faces(indices);
  vertices.eraseUnref(faces);

  const size_t NV = vertices.size();
  const size_t nv_to_decimate = std::lround(options.strength * NV);

  if (nv_to_decimate == 0) return;

  // [1] find out information of edges (endpoints, incident faces) and face2edge
  E edges;
  Edge::embedVertices(vertices);
  construct_edges(vertices, faces, edges);

  // [2] compute quadrics of vertices
  compute_quadrics(vertices, faces, edges, options);

  // [3] assigning edge errors using quadrics
  for (auto &edge : edges) edge.plan_collapse(options.fix_boundary);

  // [4] create priority queue on quadric error
  /*QEMHeap heap(edges, !options.fix_boundary);*/
  QEMHeap heap(edges);
  for (idx e = 0; e < edges.size(); ++e) {
    if (!(options.fix_boundary && edges[e].both_v_on_border())) heap.push(e);
  }
  heap.heapilize();

  size_t nv = nv_to_decimate;
  while (!heap.empty() && nv > 0) {
    // target the least-error edge, if it is what we saw last iteration,
    // it means loop should stop because all remaining edges have been penalized
    Edge *const edge = heap.top();
    if (edge->col_error() >= std::numeric_limits<double>::max()) break;

    // [5] collapse the least-error edge until mesh is simplified enough
    int collapsed = edge_collapse(vertices, faces, heap, *edge, options);

    nv -= collapsed;
  }

  // edges are useless
  // faces and vertices will be used to generate indices and positions
  // then they are useless as well
  faces.compactIndicesAndDie(indices);
  vertices.compactPositionsAndDie(positions, indices);
}

}  // namespace MeshSimpl
