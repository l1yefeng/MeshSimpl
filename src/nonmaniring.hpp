//
// Created by nickl on 5/31/19.
//

#ifndef MESH_SIMPL_NONMANIRING_HPP
#define MESH_SIMPL_NONMANIRING_HPP

#include <algorithm>
#include <iostream>
#include <map>
#include <set>

#include "edge.hpp"
#include "neighbor.hpp"
#include "proc.hpp"
#include "qemheap.hpp"
#include "vertices.hpp"

namespace MeshSimpl {
namespace Internal {

class NonManiRing {
 private:
  Vertices& vertices;
  Faces& faces;
  QEMHeap& heap;
  Edge* target;
  const SimplifyOptions& options;

  typedef std::map<idx, std::array<Edge*, 2>> NonManiInfo;

  void collect(std::vector<Neighbor>& vDelNeighbors,
               std::vector<Neighbor>& vKeptNeighbors) {
    const idx f0 = target->face(0);
    const idx f1 = target->face(1);
    const idx vDel = target->endpoint(target->delEndpointOrder());
    const idx vKept = target->endpoint(1 - target->delEndpointOrder());

    // from f0 to f1 around vDel
    Neighbor nb(target, 0, vDel, faces);
    while (true) {
      nb.rotate();
      if (nb.f() == f1) break;
      vDelNeighbors.push_back(nb);
    }

    if (target->neitherEndOnBoundary()) {
      // from f1 to f0 around vKept
      nb.replace(target, 1, vKept);
      while (true) {
        nb.rotate();
        if (nb.f() == f0) break;
        vKeptNeighbors.push_back(nb);
      }
    } else {
      // from f1 to boundary around vKept
      nb.replace(target, 1, vKept);
      while (!nb.secondEdge()->onBoundary()) {
        nb.rotate();
        vKeptNeighbors.push_back(nb);
      }

      // from f0 to boundary around vKept (switch direction)
      nb.replace(target, 0, vKept);
      while (!nb.secondEdge()->onBoundary()) {
        nb.rotate();
        vKeptNeighbors.push_back(nb);
      }
    }
  }

  void findCoincideEdges(const std::vector<Neighbor>& vDelNeighbors,
                         const std::vector<Neighbor>& vKeptNeighbors,
                         NonManiInfo& nonMani) {
    std::vector<const Neighbor*> starDel, starKept;
    starDel.reserve(vDelNeighbors.size());
    starKept.reserve(vKeptNeighbors.size());

    for (const auto& nb : vDelNeighbors) starDel.push_back(&nb);
    for (const auto& nb : vKeptNeighbors) starKept.push_back(&nb);

    const auto cmp = [&](const Neighbor* nb0, const Neighbor* nb1) -> bool {
      return nb0->secondV() < nb1->secondV();
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
        std::array<Edge*, 2> coincided = {(*itd)->secondEdge(),
                                          (*itk)->secondEdge()};
        nonMani.emplace((*itd)->secondV(), coincided);

        ++itd;
        ++itk;
      }
    }
  }

  bool checkGeom(const std::vector<Neighbor>& vDelNeighbors,
                 const std::vector<Neighbor>& vKeptNeighbors) const {
    for (const auto& nb : vDelNeighbors) {
      if (isFaceFolded(vertices, faces, nb.f(), nb.center(), target->center(),
                       options.foldOverAngleThreshold))
        return false;
    }
    for (const auto& nb : vKeptNeighbors) {
      if (isFaceFolded(vertices, faces, nb.f(), nb.center(), target->center(),
                       options.foldOverAngleThreshold))
        return false;
    }
    return true;
  }

  void cleanup(std::map<idx, NonManiInfo>& nonManiGroup) {
    assert(!nonManiGroup.empty());
    idx vKept = nonManiGroup.begin()->first;
    NonManiInfo& nonMani = nonManiGroup.begin()->second;
    assert(!nonMani.empty());

    // use the first non-manifold edge to separate this mess
    std::vector<std::tuple<Edge*, idx, idx>> edgesReplaceEnd;
    std::vector<std::tuple<idx, order, idx>> facesSetV;

    std::map<idx, int> visited;
    for (auto& nm : nonMani) {
      visited.emplace(nm.first, 0);
    }

    auto& nm = *nonMani.begin();
    // e0 will be attached to the e0->face(0) now and the face connects to it
    // on e1 e0 will keep current endpoints e1 will be attached to the
    // e0->face(1) now and the face connects to it on e1 e1 will have
    // endpoints vKept replaced with vDel, vOther with vOtherFork
    Edge* e0 = nm.second[0];
    Edge* e1 = nm.second[1];
    assert(e0->endpoints() == e1->endpoints());

    idx vKeptFork = vertices.duplicateV(vKept);
    idx vOther = nm.first;
    idx vOtherFork = vertices.duplicateV(vOther);

    // run = false => around vKept (replaced with vDel)
    // run = true  => around the other (replaced with vOtherFork)
    idx fSave;
    std::vector<Neighbor> aroundVKept, aroundVOther;

    Neighbor nb(e0, 1, vKept, faces);
    assert(faces.v(nb.f(), nb.center()) == vKept);
    while (true) {
      edgesReplaceEnd.emplace_back(nb.secondEdge(), vKept, vKeptFork);
      facesSetV.emplace_back(nb.f(), nb.center(), vKeptFork);

      auto itVisited = visited.find(nb.secondV());
      if (itVisited != visited.end()) {
        itVisited->second += 1;
      }

      if (nb.secondEdge() == e1) {
        fSave = nb.f();
        break;
      }
      nb.rotate();
      assert(nb.secondEdge() != e0);
    }

    for (auto it = nonMani.begin(), last = nonMani.end(); it != last;) {
      switch (visited[it->first]) {
        case 1:
          it = nonMani.erase(it);
          break;
        case 0:
          ++it;
          break;
        case 2: {
          auto insertRes = nonManiGroup.emplace(vKeptFork, NonManiInfo{});
          insertRes.first->second.emplace(it->first, it->second);
          it = nonMani.erase(it);
          break;
        }
        default:
          assert(false);
      }
    }

    if (nonMani.empty()) {
      nonManiGroup.erase(vKept);
    }

    nb.replace(e0, 1, vOther);
    assert(faces.v(nb.f(), nb.center()) == vOther);
    while (true) {
      edgesReplaceEnd.emplace_back(nb.secondEdge(), vOther, vOtherFork);
      facesSetV.emplace_back(nb.f(), nb.center(), vOtherFork);

      if (nb.secondEdge() == e1) {
        assert(fSave == nb.f());
        break;
      }
      nb.rotate();
      assert(nb.secondEdge() != e0);
    }

    // e0 has correct f0 but incorrect f1 (currently attached to e1)
    // e0->face(0) has correct sides
    // e0->face(1) currently has side e0 but should be replaced with e1
    // e1 has correct fSave but incorrect the other fSave2
    // fSave has correct sides
    // fSave2 currently has side e1 but should be replaced with e0
    faces.setSide(e0->face(1), e0->ordInF(1), e1);
    order fSaveOrd = e1->wingOrder(fSave);
    faces.setSide(e1->face(1 - fSaveOrd), e1->ordInF(1 - fSaveOrd), e0);
    idx fSave2 = e1->face(1 - fSaveOrd);
    order fSave2Ord = e1->ordInF(1 - fSaveOrd);
    e1->setWing(1 - fSaveOrd, e0->face(1), e0->ordInF(1));
    e0->setWing(1, fSave2, fSave2Ord);

    for (auto& ere : edgesReplaceEnd)
      std::get<0>(ere)->replaceEndpoint(std::get<1>(ere), std::get<2>(ere));
    for (auto& fsv : facesSetV)
      faces.setV(std::get<0>(fsv), std::get<1>(fsv), std::get<2>(fsv));
  }

 public:
  NonManiRing(Vertices& vertices, Faces& faces, QEMHeap& heap, Edge* target,
              const SimplifyOptions& options)
      : vertices(vertices),
        faces(faces),
        heap(heap),
        target(target),
        options(options) {}

  int collapse() {
    std::vector<Neighbor> vDelNeighbors, vKeptNeighbors;
    collect(vDelNeighbors, vKeptNeighbors);

    // special case: 2 faces and 3 vertices in current component
    if (vDelNeighbors.empty()) {
      assert(vKeptNeighbors.empty());

      idx f0 = target->face(0);
      idx f1 = target->face(1);
      for (order ord : {0, 1, 2}) {
        Edge* edge = faces.side(f0, ord);
        assert(edge->face(0) + edge->face(1) == f0 + f1);

        edge->erase();
        vertices.erase(faces.v(f0, ord));
      }

      faces.erase(f0);
      faces.erase(f1);

      return 3;
    }

    // special case: tetrahedron
    if (vDelNeighbors.size() == 1 && vKeptNeighbors.size() == 1) {
      const auto& nbd = vDelNeighbors.front();
      const auto& nbk = vKeptNeighbors.front();
      assert(faces.side(nbd.f(), nbd.center()) ==
             faces.side(nbk.f(), nbk.center()));
      for (order i : {0, 1}) {
        vertices.erase(target->endpoint(i));
        vertices.erase(faces.v(target->face(i), target->ordInF(i)));
        faces.erase(target->face(i));
      }
      nbd.firstEdge()->erase();
      nbd.secondEdge()->erase();
      nbk.firstEdge()->erase();
      nbk.secondEdge()->erase();
      target->erase();
      faces.side(nbd.f(), nbd.center())->erase();
      faces.erase(nbd.f());
      faces.erase(nbk.f());

      return 4;
    }

    // check cause of topo change
    order vDelOrd = target->delEndpointOrder();
    idx vDel = target->endpoint(vDelOrd);
    idx vKept = target->endpoint(1 - vDelOrd);
    NonManiInfo nonMani;
    findCoincideEdges(vDelNeighbors, vKeptNeighbors, nonMani);

    if (nonMani.empty()) {
      if (!checkGeom(vDelNeighbors, vKeptNeighbors)) {
        heap.penalize(target);
        return 0;
      }
    }

    // there is topo change or not, collapse the target now. cleanup afterwords
    std::vector<Edge*> dirtyEdges;

    vertices.setPosition(vKept, target->center());
    vertices.setQ(vKept, target->q());

    auto it = vDelNeighbors.begin();
    dirtyEdges.push_back(it->firstEdge());
    for (; it != vDelNeighbors.end(); ++it) {
      faces.setV(it->f(), it->center(), vKept);
      Edge* dirty = it->secondEdge();
      dirty->replaceEndpoint(vDel, vKept);
      dirtyEdges.push_back(dirty);
    }

    if (target->oneEndOnBoundary()) {
      for (order i : {0, 1})
        dirtyEdges.push_back(faces.edgeAcrossFrom(target->face(i), vDel));
      for (auto& nb : vKeptNeighbors) dirtyEdges.push_back(nb.secondEdge());
    } else {
      it = vKeptNeighbors.begin();
      dirtyEdges.push_back(it->firstEdge());
      for (; it != vKeptNeighbors.end(); ++it) {
        dirtyEdges.push_back(it->secondEdge());
      }
    }

    for (auto& dirty : dirtyEdges) {
      if (dirty->bothEndsOnBoundary()) continue;
      double errorPrev = dirty->error();
      dirty->planCollapse(options.fixBoundary);
      heap.fix(dirty, errorPrev);
    }

    // take away f0
    idx f = target->face(0);
    Edge* edgeK = faces.edgeAcrossFrom(f, vDel);
    order fOrdInEdgeK = edgeK->wingOrder(f);
    edgeK->setWing(fOrdInEdgeK, vDelNeighbors.front().f(),
                   vDelNeighbors.front().secondVOrd());
    vDelNeighbors.front().firstEdge()->erase();
    faces.setSide(vDelNeighbors.front().f(), vDelNeighbors.front().secondVOrd(),
                  edgeK);

    // take away f1
    f = target->face(1);
    edgeK = faces.edgeAcrossFrom(f, vDel);
    fOrdInEdgeK = edgeK->wingOrder(f);
    edgeK->setWing(fOrdInEdgeK, vDelNeighbors.back().f(),
                   vDelNeighbors.back().firstVOrd());
    vDelNeighbors.back().secondEdge()->erase();
    faces.setSide(vDelNeighbors.back().f(), vDelNeighbors.back().firstVOrd(),
                  edgeK);

    target->erase();
    for (order i : {0, 1}) faces.erase(target->face(i));

    vertices.erase(vDel);
    if (nonMani.empty()) {
      return 1;
    }

    int forked;
    std::map<idx, NonManiInfo> nonManiGroup = {{vKept, nonMani}};
    for (forked = 0; !nonManiGroup.empty(); forked++) {
      cleanup(nonManiGroup);
    }

    return 1 - 2 * forked;
  }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_NONMANIRING_HPP
