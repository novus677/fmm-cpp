#pragma once

#include "multipole.h"
#include "point.h"
#include "tables.h"
#include "vector.h"

#include <complex>
#include <vector>

using Complex = std::complex<double>;

struct LocalExpansion {
  int p;
  Complex center;
  std::vector<Complex> coeffs;

  InverseExponentialTable inv_exp_table;
  BinomialTable binomial_table;

  LocalExpansion(int p, Vector2 center)
      : p(p), center(center.x, center.y), coeffs(p + 1) {}

  LocalExpansion(Vector2 center, const std::vector<Complex> &coeffs) {
    this->coeffs = coeffs;
    this->p = coeffs.size() - 1;
    this->center = Complex(center.x, center.y);
  }

  LocalExpansion &operator+=(const LocalExpansion &other);
  void clear() { coeffs.assign(p + 1, 0.0); }

  double evaluate(Vector2 point) const;

  void M2L(const MultipoleExpansion &multipole);
  LocalExpansion L2L(const Complex &shift);
};
