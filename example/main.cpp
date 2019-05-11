#include "clipp.h"
#include "simplify.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace clipp;

void load_obj(const string& filename, MeshSimpl::TriMesh& mesh);

void write_to_obj(const string& filename, const MeshSimpl::TriMesh& mesh);

int main(int argc, char* argv[]) {
    string in, out;
    MeshSimpl::SimplifyOptions options;
    string weighting;

    auto cli =
        ((value("input", in)) % "input .obj file path",
         (value("output", out)) % "output file path",
         (option("-f", "--fix-boundary").set(options.fix_boundary)) %
             "do not move vertices on boundary",
         (option("-w", "--weighting") & value("strategy", weighting)) %
             "one of { uniform, by-area, by-area-inv }",
         (option("-s", "--strength") & number("ratio", options.strength)) %
             "0.8 means remove 80% vertices",
         (option("--border-constraint") &
          number("constant", options.border_constraint)) %
             "default is 2, assign larger constant to make border more "
             "reluctant to shrink",
         (option("--fold-over") &
          number("angle", options.fold_over_angle_threshold)) %
             "default is cos(160), change of angle of faces cannot be larger "
             "than this angle",
         (option("--triangle-quality") &
          number("ratio", options.aspect_ratio_at_least)) %
             "default is 0.02, aspect_ratio cannot be smaller than this value");

    if (!parse(argc, argv, cli)) {
        cout << make_man_page(cli, argv[0]);
        return 1;
    }

    if (weighting == "uniform")
        options.weighting = MeshSimpl::UNIFORM;
    else if (weighting == "by-area")
        options.weighting = MeshSimpl::BY_AREA;
    else if (weighting == "by-area-inv")
        options.weighting = MeshSimpl::BY_AREA_INV;
    else if (!weighting.empty()) {
        cout << make_man_page(cli, argv[0]);
        return 1;
    }

    // read obj file
    MeshSimpl::TriMesh mesh;
    load_obj(in, mesh);
    cout << "Loaded mesh (#V = " << mesh.vertices.size()
         << "; #F = " << mesh.indices.size() << ") from " << in << endl;

    // simplify
    const auto before = chrono::steady_clock::now();
    MeshSimpl::simplify(mesh, options);
    const auto after = chrono::steady_clock::now();
    const long duration =
        chrono::duration_cast<chrono::milliseconds>(after - before).count();
    cout << "Simplification completed (" << duration << " ms)" << endl;

    // write to obj file
    write_to_obj(out, mesh);
    cout << "Wrote mesh (#V = " << mesh.vertices.size()
         << "; #F = " << mesh.indices.size() << ") to " << out << endl;

    return 0;
}

void load_obj(const string& filename, MeshSimpl::TriMesh& mesh) {
    ifstream ifs(filename);

    mesh.vertices.resize(0);
    mesh.indices.resize(0);

    string lead;
    for (string line; getline(ifs, line);) {
        if (line[0] == '#')
            continue;

        stringstream ss(line);
        ss >> lead;

        if (lead == "v") {
            double x, y, z;
            ss >> x >> y >> z;
            mesh.vertices.push_back({x, y, z});
        } else if (lead == "f") {
            MeshSimpl::vec3i face;
            for (int i = 0; i < 3; ++i) {
                ss >> face[i];
                char ch;
                do {
                    ch = ss.get();
                } while (!ss.eof() && ch != ' ');
            }
            for (int i = 0; i < 3; ++i)
                --face[i];
            mesh.indices.push_back(face);
        }
    }

    ifs.close();
}
void write_to_obj(const string& filename, const MeshSimpl::TriMesh& mesh) {
    ofstream ofs(filename);

    const int precision = 9;

    int coord_col_len = precision + 4;
    double coord = 0.0;
    for (const auto& v : mesh.vertices)
        coord = max({coord, abs(v[0]), abs(v[1]), abs(v[2])});
    for (unsigned sz = floor(coord); sz > 0; sz /= 10, ++coord_col_len)
        ;

    int v_col_len = 1;
    for (unsigned sz = mesh.vertices.size(); sz > 0; sz /= 10, ++v_col_len)
        ;

    ofs << "#" << endl
        << "# plain triangle mesh" << endl
        << "# vertex count: " << mesh.vertices.size() << endl
        << "# face count:   " << mesh.indices.size() << endl
        << "#" << endl
        << endl;

    ofs << fixed << setprecision(precision);
    for (const auto& v : mesh.vertices) {
        ofs << "v";
        for (const double x : v)
            ofs << right << setw(coord_col_len) << x;
        ofs << endl;
    }
    ofs << endl;
    for (const auto& f : mesh.indices) {
        ofs << "f";
        for (const auto v : f)
            ofs << right << setw(v_col_len) << v + 1;
        ofs << endl;
    }
    ofs << endl;

    ofs.close();
}
