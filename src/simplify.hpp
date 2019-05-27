//
// Created by nickl on 1/8/19.
//

#ifndef MESH_SIMPL_SIMPLIFY_HPP
#define MESH_SIMPL_SIMPLIFY_HPP

#include "types.hpp"

namespace MeshSimpl {

void validateOptions(const SimplifyOptions& options);

void simplify(Positions& positions, Indices& indices,
              const SimplifyOptions& options = {});

}  // namespace MeshSimpl

#endif  // MESH_SIMPL_SIMPLIFY_HPP
