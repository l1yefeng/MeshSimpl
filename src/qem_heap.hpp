//
// Created by nickl on 1/8/19.
//

#ifndef MESH_SIMPL_QEM_HEAP_HPP
#define MESH_SIMPL_QEM_HEAP_HPP

#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

// Reference: https://algs4.cs.princeton.edu/24pq/MinPQ.java
class QEMHeap {
 public:
  // Construct a min-binary-heap with edge ecol errors as keys;
  // store a reference of the list of edges and store all handles
  explicit QEMHeap(E &edges, bool include_boundary);

  // Returns the edge id with minimum ecol error
  Edge *top() const { return &edges[keys[1]]; }

  // Remove the top edge from heap
  void pop();

  // Fix the priority of an edge after the error value is modified;
  // Param `error_prev` is used to determine the direction of priority change
  void fix(const Edge *ptr, double error_prev);

  // Suppress this edge until it is, if ever, updated next time
  void penalize(Edge *edge);

  void erase(const Edge *ptr);

  // Returns true if heap is empty
  bool empty() const { return size() == 0; }

  // Returns the size of the heap
  size_t size() const { return n; };

 private:
  std::vector<idx> keys;        // binary heap array, indexed from 1
  E &edges;                     // a reference to `edges`
  std::vector<size_t> handles;  // handles[e] is the position of edge in keys
  size_t n;                     // = heap.size() = keys.size() - 1

  // Compare function: larger error --> lower priority
  bool greater(size_t i, size_t j) const;

  // Helper function: the sole function that modifies handles and keys data
  void exchange(size_t i, size_t j);

  // For assertion purposes
  bool is_min_heap() const { return is_min_heap(1); }

  bool is_min_heap(size_t k) const;

  // Adjust the priority of node k: direction is towards higher priority
  void swim(size_t k);

  // Adjust the priority of node k: direction is towards lower priority
  void sink(size_t k);
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_QEM_HEAP_HPP
