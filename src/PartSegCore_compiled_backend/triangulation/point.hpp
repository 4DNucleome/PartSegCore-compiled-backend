//
// Created by Grzegorz Bokota on 11.10.24.
//

#ifndef PARTSEGCORE_POINT_H
#define PARTSEGCORE_POINT_H

#include <algorithm>
#include <sstream>
#include <vector>

namespace partsegcore::point {

struct Vector;

/* Point class with x and y coordinates */
struct Point {
  typedef float coordinate_t;

  coordinate_t x;
  coordinate_t y;
  bool operator==(const Point &p) const { return x == p.x && y == p.y; }
  bool operator!=(const Point &p) const { return !(*this == p); }
  Point(coordinate_t x, coordinate_t y) : x(x), y(y) {}
  Point() = default;

  bool operator<(const Point &p) const {
    if (this->y == p.y) {
      return this->x < p.x;
    }
    return this->y < p.y;
  }

  // add operator
  Vector operator+(const Point &p) const;

  // subtract operator
  Vector operator-(const Point &p) const;

  Point operator+(const Vector &v) const;

  //  Point operator/(coordinate_t f) const {
  //      return {this->x / f, this->y / f};
  //  }
  //
  //  Point operator*(coordinate_t f) const {
  //      return {this->x * f, this->y * f};
  //  }

  // Overload the << operator for Point
  friend std::ostream &operator<<(std::ostream &os, const Point &point) {
    os << "(x=" << point.x << ", y=" << point.y << ")";
    return os;
  }

  struct PointHash {
    std::size_t operator()(const Point &p) const {
      return std::hash<Point::coordinate_t>()(p.x) ^
             (std::hash<Point::coordinate_t>()(p.y) << 1);
    }
  };
};

struct Vector {
  typedef Point::coordinate_t coordinate_t;
  coordinate_t x;
  coordinate_t y;
  Vector(coordinate_t x, coordinate_t y) : x(x), y(y) {}
  Vector() = default;
  explicit Vector(const Point &p) : x(p.x), y(p.y) {}

  Vector operator+(const Vector &v) const {
    return {this->x + v.x, this->y + v.y};
  }
  Vector operator-(const Vector &v) const {
    return {this->x - v.x, this->y - v.y};
  }
  Vector operator/(coordinate_t f) const { return {this->x / f, this->y / f}; }
  Vector operator*(coordinate_t f) const { return {this->x * f, this->y * f}; }
  Vector operator-() const { return {-this->x, -this->y}; }
};

inline Vector Point::operator+(const Point &p) const {
  return {this->x + p.x, this->y + p.y};
}

inline Vector Point::operator-(const Point &p) const {
  return {this->x - p.x, this->y - p.y};
}

inline Point Point::operator+(const Vector &v) const {
  return {this->x + v.x, this->y + v.y};
}

/*Struct to represent edge of polygon with points ordered*/
struct Segment {
  Point top{};
  Point bottom{};
  Segment(Point p1, Point p2) {
    if (p1 == p2) {
      std::ostringstream oss;
      oss << "Segment cannot have two identical points: " << p1;
      throw std::invalid_argument(oss.str());
    }
    if (p1 < p2) {
      bottom = p1;
      top = p2;
    } else {
      bottom = p2;
      top = p1;
    }
  }

  [[nodiscard]] bool is_horizontal() const { return bottom.y == top.y; }
  [[nodiscard]] bool is_vertical() const { return bottom.x == top.x; }

  Segment() = default;

  bool operator<(const Segment &s) const {
    if (this->bottom == s.bottom) {
      return this->top < s.top;
    }
    return this->bottom < s.bottom;
  }

  bool operator==(const Segment &s) const {
    return this->bottom == s.bottom && this->top == s.top;
  }

  bool operator!=(const Segment &s) const { return !(*this == s); }

  /**
   * Computes the x-coordinate of a point on the line segment at a given
   * y-coordinate. This function assumes that the line segment is defined
   * between the `bottom` and `top` points.
   *
   * If the line segment is horizontal (`bottom.y == top.y`), it returns the
   * x-coordinate of `bottom`. Otherwise, it computes the x-coordinate using
   * linear interpolation.
   *
   * @param y The y-coordinate at which to find the corresponding x-coordinate
   * on the line segment.
   * @return The x-coordinate of the point on the line segment at the specified
   * y-coordinate.
   */
  [[nodiscard]] Point::coordinate_t point_on_line_x(
      Point::coordinate_t y) const {
    if (bottom.y == top.y) {
      return bottom.x;
    }
    return bottom.x +
           (y - bottom.y) * ((top.x - bottom.x) / (top.y - bottom.y));
  }

  /**
   * t = \frac{(x_1 - x_2)(x_3 - x_2) + (y_1 - y_2)(y_3 - y_2)}{(x_3 - x_2)^2 +
   *(y_3 - y_2)^2}
   **/
  [[nodiscard]] double point_projection_factor(Point p) const {
    return ((p.x - this->top.x) * (this->bottom.x - this->top.x) +
            (p.y - this->top.y) * (this->bottom.y - this->top.y)) /
           ((this->top.x - this->bottom.x) * (this->top.x - this->bottom.x) +
            (this->top.y - this->bottom.y) * (this->top.y - this->bottom.y));
  }

  [[nodiscard]] bool point_on_line(Point p) const {
    if (this->is_horizontal()) {
      return (this->bottom.x <= p.x && p.x <= this->top.x);
    }
    if (this->is_vertical()) {
      return (this->bottom.y <= p.y && p.y <= this->top.y);
    }
    auto x_cord = this->point_on_line_x(p.y);
    return (this->bottom.x <= x_cord && x_cord <= this->top.x);
  }

  /**
   * Overloads the << operator for the Segment structure, enabling
   * segments to be directly inserted into output streams.
   *
   * This operator outputs a Segment in the format:
   * [bottom=<bottom_point>, top=<top_point>]
   *
   * @param os The output stream to which the Segment is inserted.
   * @param segment The Segment instance to be inserted into the stream.
   * @return The output stream after the Segment has been inserted.
   */
  friend std::ostream &operator<<(std::ostream &os, const Segment &segment) {
    os << "[bottom=" << segment.bottom << ", top=" << segment.top << "]";
    return os;
  }

  /**
   * A hash functor for the Segment structure.
   *
   * This functor computes a hash value for a given Segment instance.
   * The hash is determined by combining the hash values of the bottom
   * and top points of the segment.
   *
   * It is required for the Segment structure to be used as a key in
   * unordered containers, such as unordered_map and unordered_set.
   */
  struct SegmentHash {
    std::size_t operator()(const Segment &segment) const {
      std::size_t h1 = Point::PointHash()(segment.bottom);
      std::size_t h2 = Point::PointHash()(segment.top);
      return h1 ^ (h2 << 1);
    }
  };
};
Point centroid(const std::vector<Point> &point_list) {
  if (point_list.empty()) {
    return {0, 0};
  }
  Point res(0, 0);
  for (auto &point : point_list) {
    res.x += point.x;
    res.y += point.y;
  }
  res.x /= float(point_list.size());
  res.y /= float(point_list.size());
  return res;
}
}  // namespace partsegcore::point

// overload of hash function for
// an unordered map and set
namespace std {
template <>
struct hash<partsegcore::point::Point> {
  size_t operator()(const partsegcore::point::Point &point) const noexcept {
    return partsegcore::point::Point::PointHash()(point);
  }
};

template <>
struct hash<partsegcore::point::Segment> {
  size_t operator()(const partsegcore::point::Segment &segment) const noexcept {
    return partsegcore::point::Segment::SegmentHash()(segment);
  }
};

}  // namespace std

#endif  // PARTSEGCORE_POINT_H
