//
// Created by Grzegorz Bokota on 11.10.24.
//

#ifndef PARTSEGCORE_POINT_H
#define PARTSEGCORE_POINT_H

#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace partsegcore {
namespace point {

/* Point class with x and y coordinates */
struct Point {
  float x;
  float y;
  bool operator==(const Point &p) const { return x == p.x && y == p.y; }
  bool operator!=(const Point &p) const { return !(*this == p); }
  Point(float x, float y) : x(x), y(y) {}
  Point() = default;

  bool operator<(const Point &p) const {
    if (this->x == p.x) {
      return this->y < p.y;
    }
    return this->x < p.x;
  }

  struct PointHash {
    std::size_t operator()(const Point &p) const {
      return std::hash<float>()(p.x) * 31 + std::hash<float>()(p.y);
    }
  };
};

/*Struct to represent edge of polygon with points ordered*/
struct Segment {
  Point left{};
  Point right{};
  Segment(Point p1, Point p2) {
    if (p1 < p2) {
      left = p1;
      right = p2;
    } else {
      left = p2;
      right = p1;
    }
  }
  Segment() = default;
};
}  // namespace point
}  // namespace partsegcore

// overload of hash function for
// unordered map and set
namespace std {
template <>
struct hash<partsegcore::point::Point> {
  std::size_t operator()(const partsegcore::point::Point &point) const {
    return partsegcore::point::Point::PointHash()(point);
  }
};
}  // namespace std

#endif  // PARTSEGCORE_POINT_H
