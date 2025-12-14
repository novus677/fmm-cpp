#include "../src/local.h"
#include "../src/multipole.h"
#include "../src/point.h"
#include "../src/vector.h"

class GravityKernel {
public:
  static double potential(const Point &source, Vector2 point) {
    Vector2 delta = source.position - point;
    if (delta == Vector2::zeros()) {
      return 0.0;
    }
    return source.strength * std::log(delta.norm());
  }

  using Multipole = MultipoleExpansion;
  using Local = LocalExpansion;
};
