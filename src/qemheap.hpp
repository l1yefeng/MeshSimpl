//
// Created by nickl on 1/8/19.
//

#ifndef MESH_SIMPL_QEMHEAP_HPP
#define MESH_SIMPL_QEMHEAP_HPP

#include <cstddef>
#include <vector>

#include "edge.hpp"
#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

// Reference: https://algs4.cs.princeton.edu/24pq/MinPQ.java
class QEMHeap {
 public:
  // Construct a min-binary-heap with edge ecol errors as keys;
  // store a reference of the list of edges and store all handles
  explicit QEMHeap(Edges &edges);

  void prioritize() {
    for (size_t k = n / 2; k >= 1; --k) sink(k);
    assert(isMinHeap());
  }

  // Returns the edge id with minimum ecol error in heap. It is possible that
  // the output is an erased edge or it has been removed from heap
  Edge *top() const { return &edges[keys[1]]; }

  // Remove the top edge from heap
  void pop();

  // Fix the priority of an edge after the error value is modified;
  // Param `errorPrev` is used to determine the direction of priority change
  void fix(const Edge *edge, double errorPrev);

  // Suppress this edge until it is, if ever, updated next time
  void penalize(Edge *edge);

  // Returns true if heap is empty
  bool empty() const { return size() == 0; }

  // Returns the size of the heap
  size_t size() const { return n; };

  // Returns true if the edge exists in heap and is not marked as removed
  bool contains(const Edge *edge) const;

  // Mark removed but does not touch actual heap data
  void markRemoved(const Edge *edge);
  void markRemovedById(idx e);
  void unmarkRemoved(const Edge *edge) { removed[id(edge)] = false; }

 private:
  std::vector<idx> keys;        // binary heap array, indexed from 1
  Edges &edges;                 // a reference to `edges`
  std::vector<size_t> handles;  // handles[e] is the position of e in keys
  size_t n;                     // = heap.size() = keys.size() - 1
  std::vector<bool> removed;    // erased from heap

  // Compare function: larger error --> lower priority
  bool greater(size_t i, size_t j) const;

  // Helper function: the sole function that modifies handles and keys data
  void exchange(size_t i, size_t j);

  // For assertion purposes
  bool isMinHeap() const { return isMinHeap(1); }

  bool isMinHeap(size_t k) const;

  // Adjust the priority of node k: direction is towards higher priority
  void swim(size_t k);

  // Adjust the priority of node k: direction is towards lower priority
  void sink(size_t k);

  idx id(const Edge *edge) const { return edge - edges.data(); }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_QEMHEAP_HPP
