//
// Created by nickl on 1/8/19.
//

#ifndef MESH_SIMPL_SIMPLIFY_HPP
#define MESH_SIMPL_SIMPLIFY_HPP

#include "types.hpp"

namespace MeshSimpl {

namespace Internal {
void validateOptions(const SimplifyOptions& options);
}  // namespace Internal

// Simplify the mesh defined by vertex positions and face indices with an
// optional configuration. Change is written in-place.
void simplify(Positions& positions, Indices& indices,
              const SimplifyOptions& options = {});

}  // namespace MeshSimpl

#endif  // MESH_SIMPL_SIMPLIFY_HPP
