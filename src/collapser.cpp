//
// Created by nickl on 5/31/19.
//

#include <algorithm>
#include <iterator>
#include <tuple>

#include "collapser.hpp"
#include "util.hpp"
#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

void Collapser::collect() {
  for (order i : {0, 1}) {
    idx v = target->endpoint(i);
    if (!vertices.isBoundary(v)) {
      // traverse around like a fan
      Neighbor nb(target, i, v, faces);
      for (nb.rotate(); nb.f() != target->face(1 - i); nb.rotate()) {
        neighbors[i].push_back(nb);
      }
    } else {
      // traverse with two rows stop at boundary
      for (order column : {0, 1}) {
        Neighbor nb(target, column, v, faces);
        while (!nb.secondEdge()->onBoundary()) {
          assert(target->ordInF(1 - column) == INVALID ||
                 nb.f() != target->face(1 - column));
          nb.rotate();
          neighbors[i].push_back(nb);
        }

        // if true, target only has face[0]
        if (target->onBoundary()) break;
      }
    }
  }
}

void Collapser::findCoincideEdges(idx vKept) {
  // comparision used in sorting and finding intersection
  const auto cmp = [](const Neighbor* nb0, const Neighbor* nb1) -> bool {
    return nb0->secondV() < nb1->secondV();
  };

  std::array<std::vector<const Neighbor*>, 2> pointers;
  for (int i : {0, 1}) {
    pointers[i].reserve(neighbors[i].size());
    for (const auto& nb : neighbors[i]) pointers[i].push_back(&nb);
    std::sort(pointers[i].begin(), pointers[i].end(), cmp);
  }

  // finding intersection
  for (auto it0 = pointers[0].begin(), it1 = pointers[1].begin();
       it0 != pointers[0].end() && it1 != pointers[1].end();) {
    if (cmp(*it0, *it1)) {
      ++it0;
    } else if (cmp(*it1, *it0)) {
      ++it1;
    } else {
      nonMani.emplace_back(vKept, (*it0)->secondV(), (*it0)->secondEdge(),
                           (*it1)->secondEdge());

      ++it0;
      ++it1;
    }
  }
}

bool Collapser::cleanup() {
  auto it = nonMani.begin();
  if (it == nonMani.end()) return false;

  // use the first non-manifold edge to separate this mess
  std::vector<std::tuple<Edge*, idx, idx>> edgesReplaceEnd;
  std::vector<std::tuple<idx, order, idx>> facesSetV;

  // e0 will be attached to the e0->face(0) now and the face connects to it
  // on e1 e0 will keep current endpoints e1 will be attached to the
  // e0->face(1) now and the face connects to it on e1 e1 will have
  // endpoints vKept replaced with vDel, vOther with vOtherFork
  Edge* e0 = it->edges[0];
  Edge* e1 = it->edges[1];
  assert(e0->endpoints() == e1->endpoints());

  idx vKept = it->vKept;
  vertices.reduceQByHalf(vKept);
  idx vKeptFork = vertices.duplicate(vKept);
  idx vOther = it->vOther;
  vertices.reduceQByHalf(vOther);
  idx vOtherFork = vertices.duplicate(vOther);
  vertices.setBoundary(vKeptFork, false);

  idx fExch0;

  // select a direction (1, but 0 will work as well), circle around
  // until met the coincided edge and every visited face will be separated
  // to the forked edge to turn the non-manifold into 2-manifold
  order traverseOrd = 0;
  Neighbor nb(e0, traverseOrd, vKept, faces);
  assert(faces.v(nb.f(), nb.center()) == vKept);
  while (true) {
    edgesReplaceEnd.emplace_back(nb.secondEdge(), vKept, vKeptFork);
    facesSetV.emplace_back(nb.f(), nb.center(), vKeptFork);

    visitNonMani(vKept, nb.secondV());

    if (nb.secondEdge() == e1) {
      fExch0 = nb.f();
      break;
    }
    if (nb.secondEdge()->onBoundary()) {
      // switch direction in order to separate edge in one pass
      edgesReplaceEnd.clear();
      facesSetV.clear();
      for (auto& nm : nonMani) {
        if (nm.vKept == vKept) {
          nm.status = 0;
        }
      }

      traverseOrd = 1 - traverseOrd;
      nb.replace(e0, traverseOrd, vKept);
      continue;
    }

    nb.rotate();
    assert(nb.secondEdge() != e0);
  }

  nb.replace(e0, traverseOrd, vOther);
  assert(e0->ordInF(traverseOrd) != INVALID);
  assert(faces.v(nb.f(), nb.center()) == vOther);
  bool hitBoundary = false;
  while (true) {
    edgesReplaceEnd.emplace_back(nb.secondEdge(), vOther, vOtherFork);
    facesSetV.emplace_back(nb.f(), nb.center(), vOtherFork);

    if (hitBoundary) assert(nb.secondEdge() != e1);

    if (nb.secondEdge() == e1) {
      assert(fExch0 == nb.f());
      assert(!hitBoundary);
      break;
    }
    if (nb.secondEdge()->onBoundary()) {
      if (!hitBoundary) {
        nb.replace(e1, e1->wingOrder(fExch0), vOther);
        hitBoundary = true;
        continue;
      } else {
        edgesReplaceEnd.emplace_back(e1, vOther, vOtherFork);
        break;
      }
    }

    nb.rotate();
    assert(nb.secondEdge() != e0);
  }

  if (vertices.isBoundary(vOther)) {
    vertices.setBoundary(vOtherFork, hitBoundary);
    vertices.setBoundary(vOther, !hitBoundary);
  }

  // exchange wing between e0 and e1. when completed, e0 will have
  // two wings remain on surface while e1 will have two wings
  // detached (endpoints are forked)
  order fExch1Ord = 1 - e1->wingOrder(fExch0);
  faces.setSide(e0->face(traverseOrd), e0->ordInF(traverseOrd), e1);
  if (!e1->onBoundary()) {
    faces.setSide(e1->face(fExch1Ord), e1->ordInF(fExch1Ord), e0);
    idx fExch1 = e1->face(fExch1Ord);
    order e1OrdInFExch = e1->ordInF(fExch1Ord);
    e1->setWing(fExch1Ord, e0->face(traverseOrd), e0->ordInF(traverseOrd));
    e0->setWing(traverseOrd, fExch1, e1OrdInFExch);
  } else {
    assert(fExch1Ord == 1);
    e1->setWing(fExch1Ord, e0->face(traverseOrd), e0->ordInF(traverseOrd));
    bool edgeValid = e0->dropWing(e0->face(traverseOrd));
    if (!edgeValid) {
      e0->erase();
    }
  }

  for (auto& ere : edgesReplaceEnd) {
    std::get<0>(ere)->replaceEndpoint(std::get<1>(ere), std::get<2>(ere));
    dirtyEdges.push_back(std::get<0>(ere));
  }
  for (auto& fsv : facesSetV)
    faces.setV(std::get<0>(fsv), std::get<1>(fsv), std::get<2>(fsv));

  updateNonManiGroup(vKept, vKeptFork);

  if (options.fixBoundary) {
    for (auto& ere : edgesReplaceEnd) {
      heap.unmarkRemoved(std::get<0>(ere));
    }
  }

  return true;
}

int Collapser::collapse(Edge* edge) {
  target = edge;
  collect();

  order delOrd = vertices.isBoundary(target->endpoint(0)) ? 1 : 0;

  bool neck = target->bothEndsOnBoundary() && !target->onBoundary();

  // check cause of topo change
  idx vDel = target->endpoint(delOrd);
  idx vKept = target->endpoint(1 - delOrd);
  findCoincideEdges(vKept);

  // reject now, before any modification that changes topology is applied
  if (!options.topologyModifiable) {
    // collapse will create non-manifold edges
    if (!nonMani.empty()) {
      return reject();
    }

    // reject if it will erase/separate the component
    if (neck) {
      std::array<bool, 2> single{false, false};
      for (order i : {0, 1}) {
        idx f = target->face(i);
        order ord = target->ordInF(i);
        single[i] = faces.side(f, next(ord))->onBoundary() &&
                    faces.side(f, prev(ord))->onBoundary();
      }

      if (single[0] == single[1]) {
        return reject();
      }
    }

    // reject if what remains is a tetrahedron
    if (target->bothEndsOnBoundary() && neighbors[0].size() == 1 &&
        neighbors[1].size() == 1) {
      return reject();
    }

    // reject if what remains are two faces folded back-to-back
    if (!vertices.isBoundary(vDel) && neighbors[delOrd].empty()) {
      return reject();
    }
  }

  // reject if topology preserves but some face will be flipped
  if (nonMani.empty()) {
    if (!checkGeom()) {
      return reject();
    }
  }

  // special case: two faces folded (#f=2, #v=3)
  if (!vertices.isBoundary(vDel) && neighbors[delOrd].empty()) {
    idx f0 = target->face(0);
    idx f1 = target->face(1);
    for (order ord : {0, 1, 2}) {
      Edge* edg = faces.side(f0, ord);
      assert(edg->face(0) + edg->face(1) == f0 + f1);

      edg->erase();
    }
    eraseF(f0);
    eraseF(f1);
    return accept();
  }

  // there is topo change or not, collapse the target now. cleanup afterwords
  // update vertex data
  vertices.setPosition(vKept, target->center());
  vertices.setQ(vKept, target->q());

  // replace face corner
  for (auto& nb : neighbors[delOrd]) {
    faces.setV(nb.f(), nb.center(), vKept);
  }

  std::array<std::vector<Edge*>, 2> initDirtyEdges;
  // collect edges who need update around endpoint 0 and 1
  for (int i : {0, 1}) {
    if (!vertices.isBoundary(target->endpoint(i))) {
      auto it = neighbors[i].begin();
      initDirtyEdges[i].push_back(it->firstEdge());
      for (; it != neighbors[i].end(); ++it) {
        initDirtyEdges[i].push_back(it->secondEdge());
      }
    } else {
      for (int column : {0, 1}) {
        initDirtyEdges[i].push_back(faces.edgeAcrossFrom(
            target->face(column), target->endpoint(1 - i)));
        if (target->onBoundary()) break;
      }
      for (auto& nb : neighbors[i]) {
        initDirtyEdges[i].push_back(nb.secondEdge());
      }
    }
  }

  // replace edge endpoint
  for (auto& edg : initDirtyEdges[delOrd]) {
    edg->replaceEndpoint(vDel, vKept);
  }

  // update error of edges
  for (auto& ide : initDirtyEdges) {
    std::move(ide.begin(), ide.end(), std::back_inserter(dirtyEdges));
  }

  // take away face 0 and 1
  std::array<bool, 2> edgeValid = {true, true};
  std::array<Edge*, 2> edgeKept = {nullptr, nullptr};
  for (int i : {0, 1}) {
    idx fDel = target->face(i);
    edgeKept[i] = faces.edgeAcrossFrom(fDel, vDel);
    Edge* edgeDel = faces.edgeAcrossFrom(fDel, vKept);

    if (!edgeDel->onBoundary()) {
      order wingOrd = 1 - edgeDel->wingOrder(fDel);
      idx fDirty = edgeDel->face(wingOrd);
      idx edgeDelOrdInFDirty = edgeDel->ordInF(wingOrd);
      edgeKept[i]->setWing(edgeKept[i]->wingOrder(fDel), fDirty,
                           edgeDelOrdInFDirty);
      faces.setSide(fDirty, edgeDelOrdInFDirty, edgeKept[i]);
    } else {
      edgeValid[i] = edgeKept[i]->dropWing(fDel);
      if (!edgeValid[i]) {
        edgeKept[i]->erase();
      }
    }
    edgeDel->erase();

    eraseF(fDel);

    if (target->onBoundary()) break;
  }

  target->erase();

  // special case: component is separated
  if (neck && edgeValid[0] && edgeValid[1]) {
    Edge* seed = edgeKept[1];
    vertices.reduceQByHalf(vKept);
    idx vFork = vertices.duplicate(vKept);
    std::vector<Neighbor> dirtyNeighbors;

    for (int column : {0, 1}) {
      Neighbor nb(seed, column, vKept, faces);
      while (true) {
        dirtyNeighbors.push_back(nb);
        visitNonMani(vKept, nb.secondV());
        if (nb.secondEdge()->onBoundary()) break;
        nb.rotate();
      }
      if (seed->onBoundary()) break;
    }

    seed->replaceEndpoint(vKept, vFork);
    dirtyEdges.push_back(seed);
    for (auto& nb : dirtyNeighbors) {
      faces.setV(nb.f(), nb.center(), vFork);
      nb.secondEdge()->replaceEndpoint(vKept, vFork);
      dirtyEdges.push_back(nb.secondEdge());
    }

    if (options.fixBoundary) {
      heap.unmarkRemoved(seed);
      for (auto& nb : dirtyNeighbors) heap.unmarkRemoved(nb.secondEdge());
    }

    updateNonManiGroup(vKept, vFork);
  }

  while (cleanup())
    ;

  for (auto& dirty : dirtyEdges) {
    double errorPrev = dirty->error();
    if (dirty->planCollapse(options.fixBoundary)) {
      heap.fix(dirty, errorPrev);
    } else {
      heap.markRemoved(dirty);
    }
  }

  return accept();
}

void Collapser::updateNonManiGroup(idx vKept, idx vFork) {
  for (auto it = nonMani.begin(); it != nonMani.end();) {
    if (it->vKept == vKept) {
      switch (it->status) {
        case 1:
          it = nonMani.erase(it);
          break;
        case 0:
          ++it;
          break;
        case 2:
          it->status = 0;
          it->vKept = vFork;
          ++it;
          break;
        default:
          assert(false);
      }
    } else {
      ++it;
    }
  }
}

void Collapser::visitNonMani(idx vKept, idx vOther) {
  for (auto& nm : nonMani) {
    if (nm.vKept == vKept && nm.vOther == vOther) {
      ++nm.status;
      assert(nm.status == 1 || nm.status == 2);
      return;
    }
  }
}

}  // namespace Internal
}  // namespace MeshSimpl
