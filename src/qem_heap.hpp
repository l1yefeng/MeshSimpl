//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_QEM_HEAP_H
#define LIB_MESH_SIMPL_QEM_HEAP_H

#include "types.hpp"

namespace MeshSimpl {

namespace Internal {

// Reference: https://algs4.cs.princeton.edu/24pq/MinPQ.java
class QEMHeap {
public:
    // Construct a min-binary-heap with edge ecol errors as keys;
    // store a reference of the list of edges and store all handles
    explicit QEMHeap(E& edges, bool include_boundary);

    // Returns the edge id with minimum ecol error
    idx top() const { return keys[1]; }

    // Remove the top edge from heap
    void pop();

    // Fix the priority of an edge after the error value is modified;
    // Param sinking should be true if the node will have lower priority after fix-up
    void fix(const Edge* ptr, double error_prev);
    void fix_all() { prioritize_all(); }

    // Suppress this edge until it is, if ever, updated next time
    void penalize(idx e);

    void erase(idx e);
    void erase_by_ptr(const Edge* ptr) { erase(static_cast<idx>(ptr - edges.data())); }

    // Returns true if heap is empty
    bool empty() const { return n == 0; }

    // Returns the size of the heap, which should be
    // the number of edges with boundary_v != BOTH
    size_t size() const { return n; };

    // Returns the keys iterator at the begin / end of active edges
    std::vector<idx>::const_iterator begin() const { return keys.begin() + 1; }
    std::vector<idx>::const_iterator end() const { return keys.end(); }

private:
    std::vector<idx> keys;
    E& edges;
    std::vector<size_t> handles;
    size_t n;

    // Compare function: larger error --> lower priority
    bool greater(size_t i, size_t j) const {
        return edges[keys[i]].error > edges[keys[j]].error;
    }

    // Helper function: the sole function that modifies handles and keys data
    void exchange(size_t i, size_t j);

    // For assertion purposes
    bool is_min_heap() const { return is_min_heap(1); }
    bool is_min_heap(size_t k) const;

    // Adjust the priority of node k: direction is towards higher priority
    void swim(size_t k);

    // Adjust the priority of node k: direction is towards lower priority
    void sink(size_t k);

    // Prioritize keys
    void prioritize_all();
};

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_QEM_HEAP_H