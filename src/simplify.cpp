//
// Created by nickl on 1/8/19.
//

#include "simplify.hpp"
#include <limits>
#include "ecol.hpp"
#include "marker.hpp"
#include "post_proc.hpp"
#include "pre_proc.hpp"

namespace MeshSimpl {

using namespace Internal;

void options_validation(const SimplifyOptions& options) {
  if (options.strength > 1)
    throw std::invalid_argument("ERROR::INVALID_OPTION: strength > 1");
  if (options.strength < 0)
    throw std::invalid_argument("ERROR::INVALID_OPTION: strength < 0");
}

bool is_valid_edge_target(const Edge& edge, const Marker& marker,
                          bool fix_boundary) {
  assert(edge[0] != edge[1]);
  assert(marker.exist_v(edge[0]) && marker.exist_v(edge[1]));
  if (fix_boundary) {
    assert(!edge.both_v_on_border());
  }

  if (edge.on_boundary()) {
    assert(marker.exist_f(edge.face(0)));
    assert(edge.both_v_on_border());
  } else {
    assert(edge.face(0) != edge.face(1));
    assert(marker.exist_f(edge.face(0)) && marker.exist_f(edge.face(1)));
  }
  return true;
}

void simplify(TriMesh& mesh, const SimplifyOptions& options) {
  options_validation(options);

  const size_t NV = mesh.vertices.size();
  const size_t nv_to_decimate =
      NV -
      std::max(static_cast<int>(std::lround((1 - options.strength) * NV)), 3);
  auto& vertices = mesh.vertices;
  auto& indices = mesh.indices;

  if (nv_to_decimate == 0) return;

  // keeps a reference to indices and init face2edge and edges
  Connectivity conn{indices, {}, {}};

  // marker keeps track of face/vertex state: exist or deleted
  Marker marker(NV, indices.size());
  marker.mark_unref_v(indices);

  // [1] find out information of edges (endpoints, incident faces) and face2edge
  construct_edges(vertices, conn);

  // [2] compute quadrics of vertices
  Q quadrics;
  compute_quadrics(vertices, conn, quadrics, options);

  // [3] assigning edge errors using quadrics
  for (auto& edge : conn.edges)
    edge.plan_collapse(vertices, quadrics, options.fix_boundary);

  // [4] create priority queue on quadric error
  QEMHeap heap(conn.edges, !options.fix_boundary);

  size_t nv = nv_to_decimate;
  while (!heap.empty() && nv > 0) {
    // target the least-error edge, if it is what we saw last iteration,
    // it means loop should stop because all remaining edges have been penalized
    const idx target = heap.top();
    if (conn.edges[target].col_error() >= std::numeric_limits<double>::max())
      break;

    const auto& edge = conn.edges[target];

    assert(is_valid_edge_target(edge, marker, options.fix_boundary));

    // [5] collapse the least-error edge until mesh is simplified enough
    order del_ord = edge.v_del_order();
    bool collapsed =
        edge_collapse(vertices, conn, quadrics, heap, target, del_ord, options);
    if (!collapsed) continue;

    // mark adjacent faces deleted
    marker.mark_f(edge.face(0));
    if (!edge.both_v_on_border())
      marker.mark_f(edge.face(1));
    else
      assert(edge.on_boundary());

    // mark one (chosen) endpoint deleted
    marker.mark_v(edge[del_ord]);

    // of course there might have been edges deleted,
    // they were removed from heap during `edge_collapse`

    --nv;
  }

  // edges are pointless from this point on, but need to fix vertices and
  // indices
  compact_data(vertices, indices, marker);
}

}  // namespace MeshSimpl
