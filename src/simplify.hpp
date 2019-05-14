//
// Created by nickl on 1/8/19.
//

#ifndef MESH_SIMPL_SIMPLIFY_HPP
#define MESH_SIMPL_SIMPLIFY_HPP

#include "types.hpp"

namespace MeshSimpl {

void options_validation(const SimplifyOptions& options);

void simplify(std::vector<vec3d>& vertices, std::vector<vec3i>& indices,
              const SimplifyOptions& options = {});

}  // namespace MeshSimpl

#endif  // MESH_SIMPL_SIMPLIFY_HPP
