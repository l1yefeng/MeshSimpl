//
// Created by nickl on 5/13/19.
//

#include "ring.hpp"
#include <algorithm>     // for sort
#include <cassert>       // for assert
#include "proc.hpp"      // for updateError, isFaceElongated
#include "qemheap.hpp"   // for QEMHeap
#include "vertices.hpp"  // for Vertices

namespace MeshSimpl {
namespace Internal {

bool Ring::checkTopo() {
  std::vector<const Neighbor *> starDel, starKept;
  starDel.reserve(vDelNeighbors.size());
  starKept.reserve(vKeptNeighbors.size());

  for (const auto &nb : vDelNeighbors) starDel.push_back(&nb);
  for (const auto &nb : vKeptNeighbors) starKept.push_back(&nb);

  const auto cmp = [&](const Neighbor *nb0, const Neighbor *nb1) -> bool {
    return nb0->secondV(_faces) < nb1->secondV(_faces);
  };
  std::sort(starDel.begin(), starDel.end(), cmp);
  std::sort(starKept.begin(), starKept.end(), cmp);
  for (auto itd = starDel.begin(), itk = starKept.begin();
       itd != starDel.end() && itk != starKept.end();) {
    if (cmp(*itd, *itk)) {
      ++itd;
    } else if (cmp(*itk, *itd)) {
      ++itk;
    } else {
      return false;
    }
  }

  return true;
}

bool Ring::checkGeom(double foldOverAngle) const {
  for (const auto &nb : vDelNeighbors) {
    if (isFaceFolded(_vertices, _faces, nb.f(), nb.center(), edge.center(),
                     foldOverAngle))
      return false;
  }
  for (const auto &nb : vKeptNeighbors) {
    if (isFaceFolded(_vertices, _faces, nb.f(), nb.center(), edge.center(),
                     foldOverAngle))
      return false;
  }
  return true;
}

bool Ring::checkQuality(double aspectRatio) const {
  for (const auto &nb : vDelNeighbors) {
    if (isFaceElongated(_vertices[vDel], _vertices[nb.firstV(_faces)],
                        _vertices[nb.secondV(_faces)], aspectRatio))
      return false;
  }
  for (const auto &nb : vDelNeighbors) {
    if (isFaceElongated(_vertices[vKept], _vertices[nb.firstV(_faces)],
                        _vertices[nb.secondV(_faces)], aspectRatio))
      return false;
  }
  return true;
}

void InteriorRing::collect() {
  const idx f0 = edge.face(0);
  const idx f1 = edge.face(1);

  Neighbor nb(f0, edge.ordInF(0), ccw);
  while (true) {
    nb.rotate(_faces);
    if (nb.f() == f1) break;
    vDelNeighbors.push_back(nb);
  }

  bool boundaryHit = false;  // useful in case vKept is on boundary
  nb = Neighbor(f1, edge.ordInF(1), ccw);
  while (true) {
    if (nb.secondEdge(_faces)->onBoundary()) {
      if (boundaryHit)
        break;  // while-loop breaks here if one vKept is on boundary
      boundaryHit = true;
      nb = Neighbor(f0, edge.ordInF(0), !ccw);
      continue;
    }

    nb.rotate(_faces);

    // this face is deleted thus unnecessary to check fold-over if f == f0
    if (nb.f() == f0) {
      // while-loop breaks here if there is no business of boundary
      assert(!boundaryHit);
      break;
    }

    vKeptNeighbors.push_back(nb);
  }
}

bool InteriorRing::checkEnv() {
  // special case: there are 2 faces, 3 vertices in current component
  // this happens if input contains such component because
  // this edge collapse implementation avoid generating them
  if (vDelNeighbors.empty()) {
    assert(_faces[edge.face(0)][edge.ordInF(0)] ==
           _faces[edge.face(1)][edge.ordInF(1)]);
    return false;
  }

  // special case: avoid collapsing a tetrahedron into folded faces
  if (vDelNeighbors.size() == 1 && vKeptNeighbors.size() == 1) return false;

  // fine target to collapse
  return true;
}

void InteriorRing::collapse(QEMHeap &heap, bool fixBoundary) {
  const idx f0 = edge.face(0);
  const idx f1 = edge.face(1);

  // update vertex position and quadric
  if (edge.neitherEndOnBoundary()) _vertices.setPosition(vKept, edge.center());
  _vertices.setQ(vKept, edge.q());

  // two kept edges in two deleted faces
  Edge *edgeKept0 = _faces.edgeAcrossFrom(f0, vDel);
  Edge *edgeKept1 = _faces.edgeAcrossFrom(f1, vDel);

  auto it = vDelNeighbors.begin();
  Edge *dirtyEdge = edgeKept0;

  // first face to process: deleted face 0
  _faces.setV(it->f(), it->center(), vKept);
  heap.erase(_faces.side(it->f(), it->j()));

  _faces.setSide(it->f(), it->j(), edgeKept0);
  dirtyEdge->replaceWing(f0, it->f(), it->j());
  _faces.erase(f0);

  _vertices.setBoundary(vDel, _vertices.isBoundary(vKept));

  // every face centered around deleted vertex
  for (++it; it != vDelNeighbors.end(); ++it) {
    _faces.setV(it->f(), it->center(), vKept);
    dirtyEdge = _faces.side(it->f(), it->j());

    dirtyEdge->replaceEndpoint(vDel, vKept);

    updateError(_vertices, heap, dirtyEdge, fixBoundary);
  }

  // the deleted face 1
  --it;
  heap.erase(_faces.side(it->f(), it->i()));
  _faces.setSide(it->f(), it->i(), edgeKept1);
  dirtyEdge = edgeKept1;
  dirtyEdge->replaceWing(f1, it->f(), it->i());
  _faces.erase(f1);

  _vertices.erase(vDel);

  // every edge centered around the kept vertex
  for (auto nb : vKeptNeighbors) {
    dirtyEdge = nb.firstEdge(_faces);
    if (fixBoundary && dirtyEdge->bothEndsOnBoundary()) continue;
    updateError(_vertices, heap, dirtyEdge, fixBoundary);
  }

  if (!(fixBoundary && edgeKept0->bothEndsOnBoundary())) {
    updateError(_vertices, heap, edgeKept0, fixBoundary);
  }
}

void BoundaryRing::collect() {
  const idx f = edge.face(0);

  Neighbor nb(f, edge.ordInF(0), ccw);
  while (true) {
    const Edge *nextEdge = nb.secondEdge(_faces);

    if (nextEdge->onBoundary()) break;

    nb.rotate(_faces);
    vDelNeighbors.push_back(nb);
  }

  nb = Neighbor(f, edge.ordInF(0), !ccw);
  while (true) {
    const Edge *nextEdge = nb.secondEdge(_faces);

    if (nextEdge->onBoundary()) break;

    nb.rotate(_faces);
    vKeptNeighbors.push_back(nb);
  }
}

bool BoundaryRing::checkEnv() {
  // special case: there is ONE triangle in current component
  if (vKeptNeighbors.empty()) return false;

  return true;
}

void BoundaryRing::collapse(QEMHeap &heap, bool fixBoundary) {
  const idx f = edge.face(0);

  // update vertex position and quadric
  _vertices.setPosition(vKept, edge.center());
  _vertices.setQ(vKept, edge.q());

  Edge *edgeKept = _faces.edgeAcrossFrom(f, vDel);
  Edge *dirtyEdge = edgeKept;

  // it is possible that vDelNeighbors is empty
  heap.erase(_faces.edgeAcrossFrom(f, vKept));
  if (vDelNeighbors.empty()) {
    dirtyEdge->dropWing(f);
  } else {
    auto it = vDelNeighbors.begin();
    _faces.setSide(it->f(), it->j(), edgeKept);
    dirtyEdge->replaceWing(f, it->f(), it->j());
  }

  _faces.erase(f);

  // every face centered around deleted vertex
  for (const auto &nb : vDelNeighbors) {
    _faces.setV(nb.f(), nb.center(), vKept);
    dirtyEdge = nb.secondEdge(_faces);

    dirtyEdge->replaceEndpoint(vDel, vKept);

    assert(!dirtyEdge->neitherEndOnBoundary());
    updateError(_vertices, heap, dirtyEdge, false);
  }

  _vertices.erase(vDel);

  updateError(_vertices, heap, edgeKept, false);
  for (auto nb : vKeptNeighbors) {
    dirtyEdge = nb.secondEdge(_faces);
    updateError(_vertices, heap, dirtyEdge, false);
  }
}

}  // namespace Internal
}  // namespace MeshSimpl
