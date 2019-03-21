//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_SIMPLIFY_HPP
#define LIB_MESH_SIMPL_SIMPLIFY_HPP

#include "types.hpp"

namespace MeshSimpl {

void options_validation(const SimplifyOptions& options);

// Mesh simplification main method.
// TODO: Ignoring the 4th and the following values (if exist) in vertices.
std::pair<V, F> simplify(const V& vertices, const F& indices,
                         const SimplifyOptions& options = {});

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_SIMPLIFY_HPP
