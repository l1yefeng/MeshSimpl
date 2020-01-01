#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <simplify.hpp>
#include <string>
#include <types.hpp>
#include <vector>

#include "clipp.h"

using namespace std;
using namespace clipp;
using namespace MeshSimpl;

void load_obj(const string& filename, vector<vec3d>& positions,
              vector<vec3i>& indices);

void write_to_obj(const string& filename, const vector<vec3d>& positions,
                  const vector<vec3i>& indices);

void specifyFixedVertices(const string& filename, vector<bool>& fixed);

int main(int argc, char* argv[]) {
  string in, out, fixedVerticesFile;
  SimplifyOptions options;

  auto cli = (
      // clang-format off
      (value("input", in))
       % "input .obj file path",
      (value("output", out))
       % "output file path",
      (option("-f", "--fix-boundary").set(options.fixBoundary))
       % "do not move vertices on boundary",
      (option("-t", "--topo-modifiable").set(options.topologyModifiable))
       % "permit change of topology during simplification",
      (option("-w", "--weight-by-area").set(options.weightByArea))
       % "quadrics are scaled by triangle area",
      (option("-s", "--strength") & number("ratio", options.strength))
       % "0.8 means remove 80% vertices",
      (option("--border-constraint") & number("constant", options.borderConstraint))
       % ("assign larger constant to make border more reluctant to shrink (default to " + to_string(options.borderConstraint) + ")"),
      (option("--aspect-ratio") & number("ratio", options.aspectRatioThreshold))
       % ("faces with aspect ratio larger than 1/ratio won't be created; assign non-positive value to disable the checking (default to " + to_string(options.borderConstraint) + ")"),
      (option("--fixed-vertices") & value("file", fixedVerticesFile))
       % "a file with a vertex number on each line, included vertices will be fixed during the simplification"
      // clang-format on
       );


  if (!parse(argc, argv, cli)) {
    cout << make_man_page(cli, argv[0]);
    return 1;
  }

  vector<vec3d> positions;
  vector<vec3i> indices;

  // read obj file
  load_obj(in, positions, indices);
  cout << "Loaded mesh (#V = " << positions.size()
       << "; #F = " << indices.size() << ") from " << in << endl;

  if (!fixedVerticesFile.empty()) {
    // fix vertices according to given file
    options.fixedVertices.resize(positions.size());

    specifyFixedVertices(fixedVerticesFile, options.fixedVertices);
  }

  // simplify
  const auto before = chrono::steady_clock::now();
  simplify(positions, indices, options);
  const auto after = chrono::steady_clock::now();
  const long duration =
      chrono::duration_cast<chrono::milliseconds>(after - before).count();
  cout << "Simplification completed (" << duration << " ms)" << endl;

  // write to obj file
  write_to_obj(out, positions, indices);
  cout << "Wrote mesh (#V = " << positions.size() << "; #F = " << indices.size()
       << ") to " << out << endl;

  return 0;
}

void load_obj(const string& filename, vector<vec3d>& positions,
              vector<vec3i>& indices) {
  ifstream ifs(filename);

  positions.resize(0);
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
      positions.push_back({x, y, z});
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

void write_to_obj(const string& filename, const vector<vec3d>& positions,
                  const vector<vec3i>& indices) {
  ofstream ofs(filename);

  const int precision = 9;

  int coord_col_len = precision + 4;
  double coord = 0.0;
  for (const auto& v : positions)
    coord = max({coord, abs(v[0]), abs(v[1]), abs(v[2])});
  for (unsigned sz = floor(coord); sz > 0; sz /= 10, ++coord_col_len)
    ;

  int v_col_len = 1;
  for (unsigned sz = positions.size(); sz > 0; sz /= 10, ++v_col_len)
    ;

  ofs << "#" << endl
      << "# plain triangle mesh" << endl
      << "# vertex count: " << positions.size() << endl
      << "# face count:   " << indices.size() << endl
      << "#" << endl
      << endl;

  ofs << fixed << setprecision(precision);
  for (const auto& v : positions) {
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

void specifyFixedVertices(const string& filename, vector<bool>& fixed) {
  vector<int> vids;
  ifstream ifs(filename);
  copy(istream_iterator<int>(ifs), istream_iterator<int>(),
       back_inserter(vids));
  ifs.close();
  for_each(vids.begin(), vids.end(), [&fixed](int v) { fixed[v - 1] = true; });
}
