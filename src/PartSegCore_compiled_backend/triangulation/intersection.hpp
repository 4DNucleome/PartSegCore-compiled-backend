//
// Created by Grzegorz Bokota on 11.10.24.
// Based on
// https://www.geeksforgeeks.org/given-a-set-of-line-segments-find-if-any-two-segments-intersect/
// and
// Mark Berg, Otfried Cheong, Marc Kreveld, Mark Overmars.
// Computational Geometry: Algorithms and Applications.
// 3rd ed., Springer Berlin, Heidelberg, 2008.
// ISBN: 978-3-540-77974-2.
// https://link.springer.com/book/10.1007/978-3-540-77974-2
// https://doi.org/10.1007/978-3-540-77974-2
//

#ifndef PARTSEGCORE_INTERSECTION_H
#define PARTSEGCORE_INTERSECTION_H

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "point.hpp"

namespace partsegcore {
namespace intersection {

struct Event {
  point::Point p;
  int index;
  bool is_left;

  Event(float x, float y, int index, bool is_left)
      : p(x, y), index(index), is_left(is_left) {}
  Event(const point::Point &p, int index, bool is_left)
      : p(p), index(index), is_left(is_left) {}
  Event() = default;

  bool operator<(const Event &e) const {
    if (p == e.p) {
      if (is_left == e.is_left) {
        return index < e.index;
      }
      return is_left > e.is_left;
    }
    return p < e.p;
  }
};

struct EventData {
  std::vector<std::size_t> begins;
  std::vector<std::size_t> ends;
};

struct PairHash {
  std::size_t operator()(const std::pair<int, int> &p) const {
    return std::hash<int>()(p.first) * 31 + std::hash<int>()(p.second);
  }
};

typedef std::map<point::Point, EventData> IntersectionEvents;

/**
 * Checks whether point q lies on the line segment defined by points p and r.
 *
 * This function determines if a given point q lies on the line segment
 * that connects points p and r. It checks if the x and y coordinates
 * of q are within the bounding box formed by p and r.
 * It works, because we use it when orientation is 0, so we know that
 * if q is in bounding box, it lies on the segment.
 *
 * @param p The first endpoint of the line segment.
 * @param q The point to check.
 * @param r The second endpoint of the line segment.
 * @return True if point q lies on the segment pr, false otherwise.
 */
bool _on_segment(const point::Point &p, const point::Point &q,
                 const point::Point &r) {
  if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
      q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y))
    return true;
  return false;
}

/**
 * Determines the orientation of the triplet (p, q, r).
 *
 * @param p The first point.
 * @param q The second point.
 * @param r The third point.
 *
 * @return 0 if p, q and r are collinear.
 *         1 if the triplet (p, q, r) is in a clockwise orientation.
 *         2 if the triplet (p, q, r) is in a counterclockwise orientation.
 */
int _orientation(const point::Point &p, const point::Point &q,
                 const point::Point &r) {
  float val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
  if (val == 0) return 0;
  return (val > 0) ? 1 : 2;
}

/**
 * Determines if two line segments intersect.
 *
 * This function checks whether two given line segments, s1 and s2, intersect.
 * It makes use of orientation and on-segment tests to evaluate intersection
 * based on the positions of the endpoints of the segments.
 *
 * @param s1 The first line segment, represented by two endpoints.
 * @param s2 The second line segment, represented by two endpoints.
 * @return True if the segments intersect, false otherwise.
 */
bool _do_intersect(const point::Segment &s1, const point::Segment &s2) {
  const point::Point &p1 = s1.bottom;
  const point::Point &q1 = s1.top;
  const point::Point &p2 = s2.bottom;
  const point::Point &q2 = s2.top;

  int o1 = _orientation(p1, q1, p2);
  int o2 = _orientation(p1, q1, q2);
  int o3 = _orientation(p2, q2, p1);
  int o4 = _orientation(p2, q2, q1);

  if (o1 != o2 && o3 != o4) return true;

  if (o1 == 0 && _on_segment(p1, p2, q1)) return true;
  if (o2 == 0 && _on_segment(p1, q2, q1)) return true;
  if (o3 == 0 && _on_segment(p2, p1, q2)) return true;
  if (o4 == 0 && _on_segment(p2, q1, q2)) return true;

  return false;
}

/**
 * @brief Checks if two segments share an endpoint.
 *
 * This function determines whether two segments, each defined by
 * two endpoints, share any endpoint. Specifically, it checks if
 * the bottom or top endpoint of the first segment is equal to the
 * bottom or top endpoint of the second segment.
 *
 * @param s1 The first segment.
 * @param s2 The second segment.
 * @return true if the segments share at least one endpoint, false otherwise.
 */
inline bool _share_endpoint(const point::Segment &s1,
                            const point::Segment &s2) {
  return s1.bottom == s2.bottom || s1.bottom == s2.top || s1.top == s2.bottom ||
         s1.top == s2.top;
}

std::set<Event>::iterator pred(std::set<Event> &s,
                               std::set<Event>::iterator it) {
  if (it == s.begin()) {
    return s.end();
  }
  return --it;
}

std::set<Event>::iterator succ(std::set<Event> &s,
                               std::set<Event>::iterator it) {
  return ++it;
}

std::unordered_set<std::pair<int, int>, PairHash> _find_intersections(
    const std::vector<point::Segment> &segments) {
  std::unordered_set<std::pair<int, int>, PairHash> intersections;
  IntersectionEvents intersection_events;
  std::vector<Event> events;
  std::set<Event> active;
  events.reserve(2 * segments.size());
  for (std::size_t i = 0; i < segments.size(); i++) {
    intersection_events[segments[i].bottom].begins.push_back(i);
    intersection_events[segments[i].top].ends.push_back(i);
    events.emplace_back(segments[i].bottom, i, true);
    events.emplace_back(segments[i].top, i, false);
  }

  for (auto &event : intersection_events) {
    for (auto &begin : event.second.begins) {
      for (auto &end : event.second.ends) {
        if (begin != end) {
          std::cout << "Bing" << std::endl;
        }
      }
    }
  }

  std::sort(events.begin(), events.end());

  //     Print events
  std::cout << "Events: " << std::endl;
  for (auto &event : events) {
    std::cout << "Event: x=" << event.p.x << " y=" << event.p.y
              << " index=" << event.index
              << " bottom=" << (event.is_left ? "true" : "false") << std::endl;
  }
  std::cout << "End of events" << std::endl << std::flush;

  for (auto &event : events) {
    if (event.index == 4) {
      std::cout << "Event 4 " << event.is_left << std::endl;
      for (auto el : active) {
        std::cout << "Active: x=" << el.p.x << " y=" << el.p.y
                  << " index=" << el.index
                  << " bottom=" << (el.is_left ? "true" : "false") << std::endl;
      }
    }

    if (event.is_left) {
      auto next = active.lower_bound(event);
      auto prev = pred(active, next);
      // we use while, because more than two segments can intersect at the same
      // point
      while (next != active.end() &&
             _do_intersect(segments[event.index], segments[next->index])) {
        // to not lost some intersection, but exclude edges, that share endpoint
        if (!_share_endpoint(segments[event.index], segments[next->index])) {
          if (event.index < next->index) {
            intersections.emplace(event.index, next->index);
          } else {
            intersections.emplace(next->index, event.index);
          }
        }
        next = succ(active, next);
      }
      // we use while, because more than two segments can intersect at the same
      // point
      while (prev != active.end() &&
             _do_intersect(segments[event.index], segments[prev->index])) {
        // to not lost some intersection, but exclude edges, that share endpoint
        if (!_share_endpoint(segments[event.index], segments[prev->index])) {
          if (event.index < prev->index) {
            intersections.emplace(event.index, prev->index);
          } else {
            intersections.emplace(prev->index, event.index);
          }
        }
        prev = pred(active, prev);
      }
      active.insert(event);
    } else {
      auto it =
          active.find(Event(segments[event.index].bottom, event.index, true));
      auto next = succ(active, it);
      auto prev = pred(active, it);
      if (next != active.end() && prev != active.end() &&
          _do_intersect(segments[next->index], segments[prev->index])) {
        if (!_share_endpoint(segments[next->index], segments[prev->index])) {
          if (next->index < prev->index) {
            intersections.emplace(next->index, prev->index);
          } else {
            intersections.emplace(prev->index, next->index);
          }
        }
      }
      active.erase(it);
    }
  }
  return intersections;
}

/**
 * @brief Finds the intersection point of two line segments, if it exists.
 *
 * This function calculates the intersection point of two given line segments.
 * Each segment is defined by two endpoints. If the segments do not intersect,
 * the function returns the point (0, 0).
 *
 * @param s1 The first line segment.
 * @param s2 The second line segment.
 * @return The intersection point of the two segments, or (0, 0) if they do not
 * intersect.
 */
point::Point _find_intersection(const point::Segment &s1,
                                const point::Segment &s2) {
  float a1, b1, c1, a2, b2, c2, det, x, y;
  a1 = s1.top.y - s1.bottom.y;
  b1 = s1.bottom.x - s1.top.x;
  c1 = a1 * s1.bottom.x + b1 * s1.bottom.y;
  a2 = s2.top.y - s2.bottom.y;
  b2 = s2.bottom.x - s2.top.x;
  c2 = a2 * s2.bottom.x + b2 * s2.bottom.y;
  det = a1 * b2 - a2 * b1;
  if (det == 0) return {0, 0};
  x = (b2 * c1 - b1 * c2) / det;
  y = (a1 * c2 - a2 * c1) / det;
  return {x, y};
}
}  // namespace intersection
}  // namespace partsegcore

#endif  // PARTSEGCORE_INTERSECTION_H
