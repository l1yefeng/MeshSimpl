//
// Created by nickl on 5/25/19.
//

#ifndef MESH_SIMPL_ERASABLE_HPP
#define MESH_SIMPL_ERASABLE_HPP

#include <cassert>
#include <vector>
#include "types.hpp"

namespace MeshSimpl {
namespace Internal {

class Erasable {
 protected:
  bool erased;

  void checkExistence() const { assert(!erased); }

 public:
  void erase() {
    checkExistence();
    erased = true;
  }

  bool exists() const { return !erased; }
};

class Erasables {
 protected:
  std::vector<bool> erased;

 public:
  explicit Erasables(size_t sz) : erased(sz, false) {}

  void erase(idx i) {
    assert(exists(i));
    erased[i] = true;
  }

  bool exists(idx i) const {
    assert(i < size());
    return !erased[i];
  }

  size_t size() const { return erased.size(); }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_ERASABLE_HPP
