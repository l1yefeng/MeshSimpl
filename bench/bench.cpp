//
// Created by nickl on 5/27/19.
//

#include <benchmark/benchmark.h>
#include <fstream>
#include <sstream>
#include "edge.hpp"
#include "proc.hpp"
#include "qemheap.hpp"
#include "simplify.hpp"
#include "vertices.hpp"

using namespace std;
using namespace MeshSimpl;
using namespace MeshSimpl::Internal;

static const string file = "/home/nickl/Desktop/before.obj";

static Positions positions;
static Indices indices;

static void load_obj(const string& filename) {
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

static void BM_Simplify(benchmark::State& state) {
  Positions P = positions;
  Indices I = indices;
  SimplifyOptions options;
  options.strength = 0.99;

  for (auto _ : state) simplify(P, I, options);
}

BENCHMARK(BM_Simplify);

static void BM_PreProcessing(benchmark::State& state) {
  Positions P = positions;
  Indices I = indices;

  SimplifyOptions options;
  options.strength = 0.99;

  for (auto _ : state) {
    Vertices vertices(P);
    Faces faces(I);
    vertices.eraseUnref(faces);

    // [1] find out information of edges (endpoints, incident faces) and
    // face2edge
    Edges edges;
    buildConnectivity(vertices, faces, edges);

    // [2] compute quadrics of vertices
    computeQuadrics(vertices, faces, edges, options);

    // [3] assigning edge errors using quadrics
    for (auto& edge : edges) edge.planCollapse(options.fixBoundary);

    // [4] create priority queue on quadric error

    QEMHeap heap(edges);
    for (idx e = 0; e < edges.size(); ++e) {
      if (!(options.fixBoundary && edges[e].bothEndsOnBoundary())) heap.push(e);
    }
    heap.heapilize();
  }
}

BENCHMARK(BM_PreProcessing);

int main(int argc, char** argv) {
  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;

  load_obj(file);

  ::benchmark::RunSpecifiedBenchmarks();
}