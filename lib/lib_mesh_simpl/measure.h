//
// Created by nickl on 1/10/19.
//

#ifndef LIB_MESH_SIMPL_MEASURE_H
#define LIB_MESH_SIMPL_MEASURE_H

#include <cassert>
#include <chrono>
#include <map>

namespace MeshSimpl {

namespace Internal {

class Measure {
public:
    void start(const std::string& key) { before.emplace(key, std::chrono::steady_clock::now()); }
    long stop(const std::string& key) {
        const auto after = std::chrono::steady_clock::now();
        auto search = before.find(key);
        assert(search != before.end());
        const auto ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(after - search->second);
        return ms.count();
    }

private:
    std::map<std::string, std::chrono::steady_clock::time_point> before;
};

} // namespace Internal

} // namespace MeshSimpl

#endif // LIB_MESH_SIMPL_MEASURE_H
