//
// Created by nickl on 7/23/19.
//

#ifndef MESH_SIMPL_QUADRIC_HPP
#define MESH_SIMPL_QUADRIC_HPP

#include <cassert>
#include <array>

#include "types.hpp"
#include "util.hpp"

namespace MeshSimpl {
namespace Internal {

class Quadric {
 public:
  Quadric() = default;

  // face defined by normal * v + d = 0
  Quadric(const vec3d normal, double d)
      : _value({normal[0] * normal[0], normal[0] * normal[1],
                normal[0] * normal[2], normal[1] * normal[1],
                normal[1] * normal[2], normal[2] * normal[2], normal[0] * d,
                normal[1] * d, normal[2] * d, d * d}) {}

  double error(const vec3d& pos) const {
    return dot({dot({_value[0], _value[1], _value[2]}, pos),
                dot({_value[1], _value[3], _value[4]}, pos),
                dot({_value[2], _value[4], _value[5]}, pos)},
               pos) +
           dot({_value[6], _value[7], _value[8]}, pos) * 2 + _value[9];
  }

  double aDeterminant() const {
    return _value[0] * (_value[3] * _value[5] - _value[4] * _value[4]) -
           _value[1] * (_value[1] * _value[5] - _value[4] * _value[2]) +
           _value[2] * (_value[1] * _value[4] - _value[3] * _value[2]);
  }

  std::pair<vec3d, double> optimal(const double aDet) const {
    assert(aDet != 0);
    const double aDetInv = 1.0 / aDet;
    const std::array<double, 6> aInv{
        (_value[3] * _value[5] - _value[4] * _value[4]) * aDetInv,
        (_value[2] * _value[4] - _value[1] * _value[5]) * aDetInv,
        (_value[1] * _value[4] - _value[2] * _value[3]) * aDetInv,
        (_value[0] * _value[5] - _value[2] * _value[2]) * aDetInv,
        (_value[1] * _value[2] - _value[0] * _value[4]) * aDetInv,
        (_value[0] * _value[3] - _value[1] * _value[1]) * aDetInv,
    };

    const vec3d b{_value[6],_value[7],_value[8]};
    const double c = _value[9];

    vec3d center = {-dot({aInv[0], aInv[1], aInv[2]}, b),
               -dot({aInv[1], aInv[3], aInv[4]}, b),
               -dot({aInv[2], aInv[4], aInv[5]}, b)};
    const double err = dot(b, center) + c;
    return {center, err};
  }

  Quadric& operator+=(const Quadric& q) {
    for (auto i = 0; i < _value.size(); ++i) _value[i] += q._value[i];
    return *this;
  }
  Quadric& operator*=(double s) {
    for (auto& v : _value) v *= s;
    return *this;
  }
  Quadric operator+(const Quadric& q) const {
    auto res = *this;
    res += q;
    return res;
  }

 private:
  std::array<double, 10> _value;
};

}  // namespace Internal
}  // namespace MeshSimpl

#endif  // MESH_SIMPL_QUADRIC_HPP
