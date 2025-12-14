#include "local.h"

LocalExpansion &LocalExpansion::operator+=(const LocalExpansion &other) {
  if (p != other.p || center != other.center) {
    throw std::runtime_error("Cannot add incompatible local expansions");
  }

  for (int k = 0; k <= p; k++) {
    coeffs[k] += other.coeffs[k];
  }
  return *this;
}

double LocalExpansion::evaluate(Vector2 point) const {
  Complex z(point.x, point.y);
  z -= center;

  Complex result = 0;
  Complex z_power = 1;

  for (int l = 0; l <= p; l++) {
    result += coeffs[l] * z_power;
    z_power *= z;
  }

  return result.real();
}

void LocalExpansion::M2L(const MultipoleExpansion &multipole) {
  // Lemma 2.2.2
  Complex shift = multipole.center - center;

  coeffs[0] = multipole.coeffs[0] * std::log(-shift);
  for (int k = 1; k <= p; k++) {
    int sign = (k % 2 == 0) ? 1 : -1;
    coeffs[0] += static_cast<double>(sign) * multipole.coeffs[k] *
                 inv_exp_table.inv_exp(shift, k);
  }

  for (int l = 1; l <= p; l++) {
    coeffs[l] = -multipole.coeffs[0] / static_cast<double>(l);
    for (int k = 1; k <= p; k++) {
      int sign = (k % 2 == 0) ? 1 : -1;
      coeffs[l] +=
          static_cast<double>(sign) * multipole.coeffs[k] *
          inv_exp_table.inv_exp(shift, k) *
          static_cast<double>(binomial_table.binomial(l + k - 1, k - 1));
    }
    coeffs[l] *= inv_exp_table.inv_exp(shift, l);
  }
}

LocalExpansion LocalExpansion::L2L(const Complex &shift) {
  // Lemma 2.2.3
  Complex new_center = center - shift;
  std::vector<Complex> new_coeffs(coeffs);

  for (int j = 0; j < p; j++) {
    for (int k = p - j - 1; k < p; k++) {
      new_coeffs[k] -= shift * new_coeffs[k + 1];
    }
  }

  return LocalExpansion(Vector2(new_center.real(), new_center.imag()),
                        new_coeffs);
}
