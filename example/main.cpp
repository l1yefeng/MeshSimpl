#include <algorithm>  // for max
#include <array>      // for array
#include <chrono>     // for milliseconds, duration_cast, operator-, stea...
#include <cmath>      // for floor
#include <fstream>
#include <iomanip>   // for operator<<, setw, setprecision
#include <iostream>  // for operator<<, basic_ostream, endl, basic_ostre...
#include <string>    // for string, operator==, operator<<, getline, ope...
#include <vector>    // for vector

#include "clipp.h"       // for operator%, operator,, group, parameter, option
#include "simplify.hpp"  // for simplify
#include "types.hpp"     // for vec3d, vec3i, SimplifyOptions, BY_AREA, BY_A...

using namespace std;
using namespace clipp;
using namespace MeshSimpl;

void load_obj(const string& filename, vector<vec3d>& vertices,
              vector<vec3i>& indices);

void write_to_obj(const string& filename, const vector<vec3d>& vertices,
                  const vector<vec3i>& indices);

int main(int argc, char* argv[]) {
  string in, out;
  SimplifyOptions options;

  auto cli =
      ((value("input", in)) % "input .obj file path",
       (value("output", out)) % "output file path",
       (option("-f", "--fix-boundary").set(options.fixBoundary)) %
           "do not move vertices on boundary",
       (option("-t", "--topo-modifiable").set(options.topologyModifiable)) %
           "permit change of topology during simplification",
       (option("-w", "--weight-by-area").set(options.weightByArea)) %
           "quadrics are scaled by triangle area",
       (option("-s", "--strength") & number("ratio", options.strength)) %
           "0.8 means remove 80% vertices",
       (option("--border-constraint") &
        number("constant", options.borderConstraint)) %
           ("assign larger constant to make border more reluctant to shrink "
            "(default to " +
            to_string(options.borderConstraint) + ")"));

  if (!parse(argc, argv, cli)) {
    cout << make_man_page(cli, argv[0]);
    return 1;
  }

  vector<vec3d> vertices;
  vector<vec3i> indices;

  // read obj file
  load_obj(in, vertices, indices);
  cout << "Loaded mesh (#V = " << vertices.size() << "; #F = " << indices.size()
       << ") from " << in << endl;

  // simplify
  const auto before = chrono::steady_clock::now();
  simplify(vertices, indices, options);
  const auto after = chrono::steady_clock::now();
  const long duration =
      chrono::duration_cast<chrono::milliseconds>(after - before).count();
  cout << "Simplification completed (" << duration << " ms)" << endl;

  // write to obj file
  write_to_obj(out, vertices, indices);
  cout << "Wrote mesh (#V = " << vertices.size() << "; #F = " << indices.size()
       << ") to " << out << endl;

  return 0;
}

void load_obj(const string& filename, vector<vec3d>& vertices,
              vector<vec3i>& indices) {
  ifstream ifs(filename);

  vertices.resize(0);
  indices.resize(0);

  string lead;
  for (string line; getline(ifs, line);) {
    if (line[0] == '#') continue;

    lead = "";
    stringstream ss(line);
    ss >> lead;

    if (lead == "v") {
      double x, y, z;
      ss >> x >> y >> z;
      vertices.push_back({x, y, z});
    } else if (lead == "f") {
      vec3i face;
      for (int i = 0; i < 3; ++i) {
        ss >> face[i];
        char ch;
        do {
          ch = ss.get();
        } while (!ss.eof() && ch != ' ');
      }
      for (int i = 0; i < 3; ++i) --face[i];
      indices.push_back(face);
    }
  }

  ifs.close();
}

void write_to_obj(const string& filename, const vector<vec3d>& vertices,
                  const vector<vec3i>& indices) {
  ofstream ofs(filename);

  const int precision = 9;

  int coord_col_len = precision + 4;
  double coord = 0.0;
  for (const auto& v : vertices)
    coord = max({coord, abs(v[0]), abs(v[1]), abs(v[2])});
  for (unsigned sz = floor(coord); sz > 0; sz /= 10, ++coord_col_len)
    ;

  int v_col_len = 1;
  for (unsigned sz = vertices.size(); sz > 0; sz /= 10, ++v_col_len)
    ;

  ofs << "#" << endl
      << "# plain triangle mesh" << endl
      << "# vertex count: " << vertices.size() << endl
      << "# face count:   " << indices.size() << endl
      << "#" << endl
      << endl;

  ofs << fixed << setprecision(precision);
  for (const auto& v : vertices) {
    ofs << "v";
    for (const double x : v) ofs << right << setw(coord_col_len) << x;
    ofs << endl;
  }
  ofs << endl;
  for (const auto& f : indices) {
    ofs << "f";
    for (const auto v : f) ofs << right << setw(v_col_len) << v + 1;
    ofs << endl;
  }
  ofs << endl;

  ofs.close();
}
