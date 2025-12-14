#include "multipole.h"

MultipoleExpansion &
MultipoleExpansion::operator+=(const MultipoleExpansion &other) {
  if (p != other.p || center != other.center) {
    throw std::runtime_error("Cannot add incompatible multipole expansions");
  }

  for (int k = 0; k <= p; k++) {
    coeffs[k] += other.coeffs[k];
  }
  return *this;
}

double MultipoleExpansion::evaluate(Vector2 point) const {
  Complex z(point.x, point.y);
  z -= center;

  Complex result = coeffs[0].real() * std::log(z);

  Complex z_inv_power = 1.0 / z;
  for (int k = 1; k <= p; k++) {
    result += coeffs[k] * z_inv_power;
    z_inv_power /= z;
  }

  return result.real();
}

void MultipoleExpansion::buildExpansion(const std::vector<Point> &sources) {
  // Theorem 2.1.1
  for (const Point &source : sources) {
    Complex z(source.position.x, source.position.y);
    z -= center;
    coeffs[0] += source.strength;

    Complex z_power = z;
    for (int k = 1; k <= p; k++) {
      coeffs[k] -= source.strength * z_power / static_cast<double>(k);
      z_power *= z;
    }
  }
}

MultipoleExpansion MultipoleExpansion::M2M(const Complex &shift) {
  // Lemma 2.2.1
  Complex new_center = center - shift;
  std::vector<Complex> new_coeffs(p + 1);

  new_coeffs[0] = coeffs[0].real();
  for (int l = 1; l <= p; l++) {
    new_coeffs[l] =
        -coeffs[0].real() * exp_table.exp(shift, l) / static_cast<double>(l);
    for (int k = 1; k <= l; k++) {
      new_coeffs[l] +=
          coeffs[k] * exp_table.exp(shift, l - k) *
          static_cast<double>(binomial_table.binomial(l - 1, k - 1));
    }
  }

  return MultipoleExpansion(Vector2(new_center.real(), new_center.imag()),
                            new_coeffs);
}
