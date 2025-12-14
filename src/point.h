#pragma once

#include "vector.h"

#include <cmath>
#include <limits>
#include <vector>

struct Point {
  Vector2 position;
  double strength;

  Point(Vector2 position, double strength)
      : position(position), strength(strength) {}
};

struct Box2 {
  Vector2 center;
  double half_side;

  Box2(Vector2 center, double half_side)
      : center(center), half_side(half_side) {}

  bool contains(const Point &p) const {
    return std::abs(p.position.x - center.x) <= half_side &&
           std::abs(p.position.y - center.y) <= half_side;
  }
};

inline Box2 computeBoundingBox(const std::vector<Point> &points) {
  if (points.empty()) {
    return Box2{Vector2{0.0, 0.0}, 0.0};
  }

  double min_x = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::lowest();

  for (const Point &p : points) {
    min_x = std::min(min_x, p.position.x);
    max_x = std::max(max_x, p.position.x);
    min_y = std::min(min_y, p.position.y);
    max_y = std::max(max_y, p.position.y);
  }

  double center_x = 0.5 * (min_x + max_x);
  double center_y = 0.5 * (min_y + max_y);
  double half_side = 0.5 * std::max(max_x - min_x, max_y - min_y);

  constexpr double padding = 1e-5;
  half_side *= (1.0 + padding);

  return Box2{Vector2{center_x, center_y}, half_side};
}

inline bool adjacent(const Box2 &box1, const Box2 &box2) {
  double dx = std::abs(box1.center.x - box2.center.x);
  double dy = std::abs(box1.center.y - box2.center.y);

  double side_sum = box1.half_side + box2.half_side;

  const double tolerance = std::min(box1.half_side, box2.half_side) * 1e-5;
  return dx <= side_sum + tolerance && dy <= side_sum + tolerance;
}

// Helper functions for quadtree operations
inline uint32_t getQuadrant(const Box2 &box, Vector2 position) {
  // Quadrant layout:
  // 2 | 3
  // --|--
  // 0 | 1
  uint32_t quadrant = 0;
  if (position.x >= box.center.x) {
    quadrant |= 1;
  }
  if (position.y >= box.center.y) {
    quadrant |= 2;
  }
  return quadrant;
}

inline Box2 getChildBox(const Box2 &parent, uint32_t quadrant) {
  switch (quadrant) {
  case 0: {
    Vector2 center = Vector2{parent.center.x - 0.5 * parent.half_side,
                             parent.center.y - 0.5 * parent.half_side};
    return Box2{center, 0.5 * parent.half_side};
  }
  case 1: {
    Vector2 center = Vector2{parent.center.x + 0.5 * parent.half_side,
                             parent.center.y - 0.5 * parent.half_side};
    return Box2{center, 0.5 * parent.half_side};
  }
  case 2: {
    Vector2 center = Vector2{parent.center.x - 0.5 * parent.half_side,
                             parent.center.y + 0.5 * parent.half_side};
    return Box2{center, 0.5 * parent.half_side};
  }
  case 3: {
    Vector2 center = Vector2{parent.center.x + 0.5 * parent.half_side,
                             parent.center.y + 0.5 * parent.half_side};
    return Box2{center, 0.5 * parent.half_side};
  }
  default: {
    throw std::runtime_error("Invalid quadrant");
  }
  }
}
