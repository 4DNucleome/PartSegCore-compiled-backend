//
// Created by Grzegorz Bokota on 11.10.24.
//

#ifndef PARTSEGCORE_POINT_H
#define PARTSEGCORE_POINT_H

#include <algorithm>
#include <iostream>
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
      return this->y > p.y;
    }
    return this->x < p.x;
  }

  // Overload the << operator for Point
  friend std::ostream &operator<<(std::ostream &os, const Point &point) {
    os << "(x=" << point.x << ", y=" << point.y << ")";
    return os;
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

  bool operator<(const Segment &s) const {
    if (this->left == s.left) {
      return this->right < s.right;
    }
    return this->left < s.left;
  }

  bool operator==(const Segment &s) const {
    return this->left == s.left && this->right == s.right;
  }

  // Overload the << operator for Segment
  friend std::ostream &operator<<(std::ostream &os, const Segment &segment) {
    os << "[" << segment.left << " -- " << segment.right << "]";
    return os;
  }

  struct SegmentHash {
    std::size_t operator()(const Segment &segment) const {
      std::size_t h1 = Point::PointHash()(segment.left);
      std::size_t h2 = Point::PointHash()(segment.right);
      return h1 ^ (h2 << 1);
    }
  };
};
}  // namespace point
}  // namespace partsegcore

// overload of hash function for
// unordered map and set
namespace std {
template <>
struct hash<partsegcore::point::Point> {
  size_t operator()(const partsegcore::point::Point &point) const {
    return partsegcore::point::Point::PointHash()(point);
  }
};

template <>
struct hash<partsegcore::point::Segment> {
  size_t operator()(const partsegcore::point::Segment &segment) const {
    return partsegcore::point::Segment::SegmentHash()(segment);
  }
};

}  // namespace std

#endif  // PARTSEGCORE_POINT_H
