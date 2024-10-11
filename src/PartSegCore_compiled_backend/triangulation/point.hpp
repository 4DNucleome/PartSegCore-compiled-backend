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
  Point(float x, float y) : x(x), y(y) {}
  Point() = default;

  bool operator<(const Point &p) const {
    if (this->x == p.x) {
      return this->y < p.y;
    }
    return this->x < p.x;
  }
};

struct PointHash {
  std::size_t operator()(const Point &p) const {
    return std::hash<float>()(p.x) * 31 + std::hash<float>()(p.y);
  }
};

bool point_eq(const Point &p, const Point &q) {
  return p.x == q.x && p.y == q.y;
}

bool cmp_pair_point(const std::pair<Point, int> &p,
                    const std::pair<Point, int> &q) {
  return p.first < q.first;
}

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

bool cmp_point(const Point &p, const Point &q) { return p < q; }
}  // namespace point
}  // namespace partsegcore

#endif  // PARTSEGCORE_POINT_H
