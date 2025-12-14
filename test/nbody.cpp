#include "../src/fmmtree.h"
#include "../src/point.h"
#include "../src/vector.h"
#include "kernels.h"

#include <iostream>
#include <random>
#include <vector>

template <class Kernel>
std::vector<double> reference(const std::vector<Point> &sources) {
  int num_sources = sources.size();
  std::vector<double> potentials(num_sources);

  for (int i = 0; i < num_sources; i++) {
    for (int j = 0; j < num_sources; j++) {
      if (i == j) {
        continue;
      }

      potentials[i] += Kernel::potential(sources[j], sources[i].position);
    }
  }
  return potentials;
}

int main() {
  int num_sources = 10000;
  std::vector<Point> sources;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for (int i = 0; i < num_sources; i++) {
    Point rand_point(Vector2(dist(gen), dist(gen)), dist(gen));
    sources.push_back(rand_point);
  }

  int height = ceil(std::log(num_sources) / std::log(4));
  int p = 5;
  NaiveFmmTree<GravityKernel> fmm_tree(p, sources, height);

  std::vector<double> fmm_potentials = fmm_tree.evaluateSources();
  std::vector<double> reference_potentials = reference<GravityKernel>(sources);

  double max_error = 0.0;
  double sum_error = 0.0;
  for (int i = 0; i < num_sources; i++) {
    double error = std::abs(fmm_potentials[i] - reference_potentials[i]) /
                   std::abs(reference_potentials[i]);
    max_error = std::max(max_error, error);
    sum_error += error;
  }

  std::cout << "Max relative error: " << max_error << std::endl;
  std::cout << "Average relative error: " << sum_error / num_sources
            << std::endl;

  return 0;
}
