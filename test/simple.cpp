#include "../src/fmmtree.h"
#include "../src/point.h"
#include "../src/vector.h"
#include "kernels.h"

#include <iostream>
#include <random>
#include <vector>

template <class Kernel>
double reference(const std::vector<Point> &sources, const Vector2 &point) {
  double result = 0.0;

  for (const Point &source : sources) {
    result += Kernel::potential(source, point);
  }

  return result;
}

int main() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  int num_sources = 10000;
  int height = ceil(std::log(num_sources) / std::log(4));
  int p = 5;

  int num_experiments = 20;
  double max_relative_error = 0.0;
  double sum_relative_error = 0.0;
  for (int i = 0; i < num_experiments; i++) {
    std::vector<Point> sources;
    for (int i = 0; i < num_sources; i++) {
      Point rand_point(Vector2(dist(gen), dist(gen)), dist(gen));
      sources.push_back(rand_point);
    }

    NaiveFmmTree<GravityKernel> fmm_tree(p, sources, height);

    Vector2 test_point(0.5, 0.5);
    double fmm_result = fmm_tree.evaluate(test_point);
    double reference_result = reference<GravityKernel>(sources, test_point);

    double relative_error =
        std::abs(fmm_result - reference_result) / std::abs(reference_result);
    max_relative_error = std::max(max_relative_error, relative_error);
    sum_relative_error += relative_error;
  }
  std::cout << "Max relative error: " << max_relative_error << std::endl;
  std::cout << "Average relative error: "
            << sum_relative_error / num_experiments << std::endl;

  return 0;
}
