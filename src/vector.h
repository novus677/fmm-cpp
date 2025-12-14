#pragma once

#include <cmath>

struct Vector2 {
  double x, y;

  Vector2(double x, double y) : x(x), y(y) {}

  Vector2 operator+(const Vector2 &other) const {
    return Vector2(x + other.x, y + other.y);
  }

  Vector2 operator-(const Vector2 &other) const {
    return Vector2(x - other.x, y - other.y);
  }

  Vector2 operator*(const double &scalar) const {
    return Vector2(x * scalar, y * scalar);
  }

  Vector2 operator/(const double &scalar) const {
    return Vector2(x / scalar, y / scalar);
  }

  Vector2 &operator+=(const Vector2 &other) {
    x += other.x;
    y += other.y;
    return *this;
  }

  Vector2 &operator-=(const Vector2 &other) {
    x -= other.x;
    y -= other.y;
    return *this;
  }

  Vector2 &operator*=(const double &scalar) {
    x *= scalar;
    y *= scalar;
    return *this;
  }

  Vector2 &operator/=(const double &scalar) {
    x /= scalar;
    y /= scalar;
    return *this;
  }

  bool operator==(const Vector2 &other) const {
    return x == other.x && y == other.y;
  }
  bool operator!=(const Vector2 &other) const { return !(*this == other); }

  double dot(const Vector2 &other) const { return x * other.x + y * other.y; }
  double norm2() const { return x * x + y * y; }
  double norm() const { return std::sqrt(norm2()); }

  static Vector2 zeros() { return Vector2(0.0, 0.0); }
  static Vector2 ones() { return Vector2(1.0, 1.0); }
};
