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
  bool _erased;

  void checkExistence() const { assert(!_erased); }

 public:
  void erase() {
    checkExistence();
    _erased = true;
  }

  bool exists() const { return !_erased; }
};

class Erasables {
 protected:
  std::vector<bool> _erased;

 public:
  explicit Erasables(size_t sz) : _erased(sz, false) {}

  void erase(idx i) {
    assert(exists(i));
    _erased[i] = true;
  }

  bool exists(idx i) const {
    assert(i < size());
    return !_erased[i];
  }

  size_t size() const { return _erased.size(); }
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_ERASABLE_HPP
