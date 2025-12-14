#pragma once

#include "point.h"
#include "tables.h"
#include "vector.h"

#include <complex>
#include <vector>

using Complex = std::complex<double>;

struct MultipoleExpansion {
  int p;
  Complex center;
  std::vector<Complex> coeffs;

  BinomialTable binomial_table;
  ExponentialTable exp_table;

  MultipoleExpansion(int p, Vector2 center)
      : p(p), center(center.x, center.y), coeffs(p + 1) {}

  MultipoleExpansion(Vector2 center, const std::vector<Complex> &coeffs) {
    this->coeffs = coeffs;
    this->p = coeffs.size() - 1;
    this->center = Complex(center.x, center.y);
  }

  MultipoleExpansion &operator+=(const MultipoleExpansion &other);
  void clear() { coeffs.assign(p + 1, 0.0); }

  double evaluate(Vector2 point) const;

  void buildExpansion(const std::vector<Point> &sources);
  MultipoleExpansion M2M(const Complex &shift);
};
