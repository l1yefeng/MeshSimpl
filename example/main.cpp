#include "clipp.h"
#include "simplify.hpp"
#include <chrono>
#include <cstdlib>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace clipp;

int main(int argc, char* argv[]) {
    string in, out;
    bool fix_boundary;
    string weighting;
    float strength;

    auto cli = (value("input", in).doc("input .obj file"),
                option("-o", "--output") & word("path", out).doc("output file path"),
                option("-f", "--fix-boundary")
                    .set(fix_boundary)
                    .doc("do not move vertices on boundary"),
                option("-w", "--weighting") &
                    word("strategy", weighting).doc("one of [area, uniform]"),
                option("-s", "--strength") &
                    number("ratio", strength).doc("0.8 means remove 80% vertices"));

    if (!parse(argc, argv, cli)) {
        cout << make_man_page(cli, argv[0]);
        return 1;
    }

    MeshSimpl::SimplifyOptions options;
    options.strength = strength;
    options.fix_boundary = fix_boundary;
    if (weighting == "area")
        options.weighting = MeshSimpl::BY_AREA;

    // read obj file
    vector<vector<double>> vertices;
    vector<vector<unsigned int>> indices;
    igl::readOBJ(in, vertices, indices);

    const auto before = chrono::steady_clock::now();

    try {
        auto res = MeshSimpl::simplify(vertices, indices, options);
        const auto after = chrono::steady_clock::now();
        const long duration =
            chrono::duration_cast<chrono::milliseconds>(after - before).count();
        std::cout << "[INFO] Simplification completed (" << duration << " milliseconds)"
                  << std::endl;

        Eigen::MatrixXd V(res.first.size(), 3);
        for (int i = 0; i < V.rows(); ++i)
            V.row(i) << res.first[i][0], res.first[i][1], res.first[i][2];
        Eigen::MatrixXi F(res.second.size(), 3);
        for (int i = 0; i < F.rows(); ++i)
            F.row(i) << res.second[i][0], res.second[i][1], res.second[i][2];

        if (out.empty()) {
            // view
            igl::opengl::glfw::Viewer viewer;
            viewer.data().set_mesh(V, F);
            viewer.launch();
        } else {
            // write
            igl::writeOBJ(out, V, F);
        }
    } catch (char const* exception) {
        cerr << exception << endl;
        return 1;
    }
}
