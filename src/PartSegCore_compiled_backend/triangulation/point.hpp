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
    if (this->y == p.y) {
      return this->x < p.x;
    }
    return this->y < p.y;
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
  Point bottom{};
  Point top{};
  Segment(Point p1, Point p2) {
    if (p1 < p2) {
      bottom = p1;
      top = p2;
    } else {
      bottom = p2;
      top = p1;
    }
  }

  bool is_horizontal() const { return bottom.y == top.y; }

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

  // Overload the << operator for Segment
  friend std::ostream &operator<<(std::ostream &os, const Segment &segment) {
    os << "[bottom=" << segment.bottom << ", top=" << segment.top << "]";
    return os;
  }

  struct SegmentHash {
    std::size_t operator()(const Segment &segment) const {
      std::size_t h1 = Point::PointHash()(segment.bottom);
      std::size_t h2 = Point::PointHash()(segment.top);
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
