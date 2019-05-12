//
// Created by nickl on 1/8/19.
//

#ifndef LIB_MESH_SIMPL_SIMPLIFY_HPP
#define LIB_MESH_SIMPL_SIMPLIFY_HPP

#include "types.hpp"

namespace MeshSimpl {

void options_validation(const SimplifyOptions& options);

void simplify(TriMesh& mesh, const SimplifyOptions& options = {});

// Mesh simplification main method.
bool simplify(const TriMesh& input, TriMesh& output,
              const SimplifyOptions& options = {});

}  // namespace MeshSimpl

#endif  // LIB_MESH_SIMPL_SIMPLIFY_HPP
