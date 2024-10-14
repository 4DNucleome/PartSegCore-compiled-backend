#ifndef PARTSEGCORE_TRIANGULATE_H
#define PARTSEGCORE_TRIANGULATE_H

#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "intersection.hpp"
#include "point.hpp"

namespace partsegcore {
namespace triangulation {

enum PointType { NORMAL, SPLIT, MERGE, INTERSECTION };

struct Interval {
  point::Point last_seen;
  Interval() = default;
  explicit Interval(const point::Point &p) : last_seen(p) {};
};

/* Comparator for segments
 * To determine if segment is left or right of other segment
 *
 * It assumes that segments are not intersecting
 */

struct SegmentLeftRightComparator {
  bool operator()(const point::Segment &s1, const point::Segment &s2) const {
    return false;
  }
};

struct OrderedPolygon {
  point::Point top;
  point::Point bottom;
  std::vector<point::Point> left;
  std::vector<point::Point> right;
};

struct Triangle {
  std::size_t x;
  std::size_t y;
  std::size_t z;
  Triangle(std::size_t x, std::size_t y, std::size_t z) : x(x), y(y), z(z) {};
  Triangle() = default;
};

typedef std::size_t EdgeIndex;

struct PointEdges {
  EdgeIndex edge_index;
  point::Point opposite_point;
  PointEdges(EdgeIndex edge_index, point::Point opposite_point)
      : edge_index(edge_index), opposite_point(opposite_point) {}
  bool operator<(const PointEdges &e) const {
    return opposite_point < e.opposite_point;
  }
};

typedef std::unordered_map<point::Point, std::vector<PointEdges>> PointToEdges;
typedef std::map<point::Segment, Interval *, SegmentLeftRightComparator>
    SegmentToLine;

/**
 * Checks if a given polygon is convex.
 *
 * This function takes a polygon represented as a vector of points and
 * determines if it is convex. A convex polygon is one where all the interior
 * angles are less than 180 degrees, meaning that the edges never turn back on
 * themselves.
 *
 * @param polygon A vector of points representing the vertices of the polygon in
 * order.
 * @return True if the polygon is convex, false otherwise.
 */
bool _is_convex(const std::vector<point::Point> &polygon) {
  int orientation = 0;
  int triangle_orientation;
  for (std::size_t i = 0; i < polygon.size() - 2; i++) {
    triangle_orientation =
        intersection::_orientation(polygon[i], polygon[i + 1], polygon[i + 2]);
    if (triangle_orientation == 0) continue;
    if (orientation == 0)
      orientation = triangle_orientation;
    else if (orientation != triangle_orientation)
      return false;
  }
  triangle_orientation = intersection::_orientation(
      polygon[polygon.size() - 2], polygon[polygon.size() - 1], polygon[0]);
  if (triangle_orientation != 0 && triangle_orientation != orientation)
    return false;
  triangle_orientation = intersection::_orientation(polygon[polygon.size() - 1],
                                                    polygon[0], polygon[1]);
  if (triangle_orientation != 0 && triangle_orientation != orientation)
    return false;
  return true;
}

/**
 * Divides a convex polygon into triangles using the fan triangulation method.
 * Starting from the first point in the polygon, each triangle is formed by
 * connecting the first point with two consecutive points from the polygon list.
 *
 * @param polygon A vector of Point objects representing the vertices of the
 * polygon.
 * @return A vector of Triangle objects where each triangle is defined by
 * indices corresponding to the vertices in the input polygon.
 */
std::vector<Triangle> _triangle_convex_polygon(
    const std::vector<point::Point> &polygon) {
  std::vector<Triangle> result;
  for (std::size_t i = 1; i < polygon.size() - 1; i++) {
    if (intersection::_orientation(polygon[0], polygon[i], polygon[i + 1]) !=
        0) {
      result.emplace_back(0, i, i + 1);
    }
  }
  return result;
}

std::vector<point::Segment> calc_edges(
    const std::vector<point::Point> &polygon) {
  std::vector<point::Segment> edges;
  edges.reserve(polygon.size());
  for (std::size_t i = 0; i < polygon.size() - 1; i++) {
    edges.emplace_back(polygon[i], polygon[i + 1]);
  }
  edges.emplace_back(polygon[polygon.size() - 1], polygon[0]);
  return edges;
}

std::vector<point::Point> find_intersection_points(
    const std::vector<point::Point> &polygon) {
  /* find all edge intersections and add mid-points for all such intersection
   * place*/
  auto edges = calc_edges(polygon);

  auto intersections = intersection::_find_intersections(edges);
  if (intersections.empty()) return polygon;
  std::unordered_map<std::size_t, std::vector<point::Point>>
      intersections_points;
  for (const auto &intersection : intersections) {
    auto inter_point = intersection::_find_intersection(
        edges[intersection.first], edges[intersection.second]);
    intersections_points[intersection.first].push_back(inter_point);
    intersections_points[intersection.second].push_back(inter_point);
  }
  std::size_t points_count = polygon.size();
  for (auto &intersections_point : intersections_points) {
    points_count += intersections_point.second.size() - 1;
    intersections_point.second.push_back(
        edges[intersections_point.first].right);
    intersections_point.second.push_back(edges[intersections_point.first].left);
    std::sort(intersections_point.second.begin(),
              intersections_point.second.end());
  }

  std::vector<point::Point> new_polygon;
  new_polygon.reserve(points_count);
  for (std::size_t i = 0; i < polygon.size(); i++) {
    auto point = polygon[i];
    new_polygon.push_back(point);
    if (intersections_points.count(i)) {
      auto new_points = intersections_points[i];
      if (new_points[0] == point) {
        for (std::size_t j = 1; j < new_points.size() - 1; j++) {
          new_polygon.push_back(new_points[j]);
        }
      } else {
        for (std::size_t j = new_points.size() - 2; j > 0; j++) {
          new_polygon.push_back(new_points[j]);
        }
      }
    }
  }
  return new_polygon;
}

/*
    Calculate point type.
    If there is more than two edges adjusted to point, it is intersection point.
    If there are two adjusted edges, it could be one of split, merge and normal
    point. If both adjusted edges have opposite end before given point p, this
   is merge point. If both adjusted edges have opposite end after given point p,
   split point. Otherwise it is normal point.
    */
PointType get_point_type(point::Point p, PointToEdges &point_to_edges) {
  if (point_to_edges.at(p).size() != 2) return PointType::INTERSECTION;
  const auto &edges = point_to_edges.at(p);
  if (edges[0].opposite_point < p && edges[1].opposite_point < p)
    return PointType::MERGE;
  if (p < edges[0].opposite_point && p < edges[1].opposite_point)
    return PointType::SPLIT;
  return PointType::NORMAL;
}

/*
    Get map from point to list of edges which contains this point.
    Also sort each list by point order.
    */
PointToEdges get_points_edges(std::vector<point::Segment> &edges) {
  PointToEdges point_to_edges;
  for (std::size_t i = 0; i < edges.size(); i++) {
    point_to_edges[edges[i].left].emplace_back(i, edges[i].right);
    point_to_edges[edges[i].right].emplace_back(i, edges[i].left);
  }
  for (auto &point_to_edge : point_to_edges) {
    std::sort(point_to_edge.second.begin(), point_to_edge.second.end());
  }
  return point_to_edges;
}

void _process_normal_point(
    const point::Point &p, const std::vector<point::Segment> &edges,
    const PointToEdges &point_to_edges,
    std::map<point::Segment, Interval *> &segment_to_line) {
  const point::Segment &edge_prev =
      edges[point_to_edges.at(p).at(0).edge_index];
  const point::Segment &edge_next =
      edges[point_to_edges.at(p).at(1).edge_index];
  segment_to_line[edge_next] = segment_to_line.at(edge_prev);
  segment_to_line.at(edge_prev)->last_seen = p;
  segment_to_line.erase(edge_prev);
}

void _process_split_point(
    const point::Point &p, const std::vector<point::Segment> &edges,
    const PointToEdges &point_to_edges,
    std::map<point::Segment, Interval *> &segment_to_line) {
  const point::Segment &edge_left =
      edges[point_to_edges.at(p).at(0).edge_index];
  const point::Segment &edge_right =
      edges[point_to_edges.at(p).at(1).edge_index];
}

/* process merge point
 * When merge point is found, we need to merge two intervals into one
 *
 */
void _process_merge_point(
    const point::Point &p, const std::vector<point::Segment> &edges,
    const PointToEdges &point_to_edges,
    std::map<point::Segment, Interval *> &segment_to_line) {
  const point::Segment &edge_left =
      edges[point_to_edges.at(p).at(0).edge_index];
  const point::Segment &edge_right =
      edges[point_to_edges.at(p).at(1).edge_index];
}

///* this is processing of point that have more than 2 edges adjusted */
// void _process_intersection_point(
//     const point::Point &p, std::unordered_map<point::Segment, Interval>
//     &segment_to_line, std::unordered_set<point::Segment>
//     &sweeping_line_intersect) {
//   for (auto &edge : point_to_edges[p]) {
//     sweeping_line_intersect.insert(edges[edge.edge_index]);
//   }
//   std::vector<Interval> intervals;
//   for (auto &edge : point_to_edges[p]) {
//     intervals.push_back(Interval(p, edges[edge.edge_index],
//     edges[edge.edge_index]));
//   }
// }

/*
    This is implementation of sweeping line triangulation of polygon
    Its assumes that there is no edge intersections, but may be a point with
   more than 2 edges. described on this lecture:
    https://www.youtube.com/playlist?list=PLtTatrCwXHzEqzJMaTUFgqoCNllgwk4DH
    */
std::pair<std::vector<Triangle>, std::vector<point::Point>>
sweeping_line_triangulation(const std::vector<point::Point> &polygon) {
  std::vector<Triangle> result;
  point::Segment *edge_prev, *edge_next;
  auto edges = calc_edges(polygon);
  PointToEdges point_to_edges = get_points_edges(edges);

  std::vector<point::Point> sorted_points = polygon;
  // copy to avoid modification of original vector
  std::sort(sorted_points.begin(), sorted_points.end());
  std::vector<OrderedPolygon> ordered_polygon_li;
  std::map<point::Segment, Interval *> segment_to_line;
  std::vector<Interval> intervals;
  ordered_polygon_li.emplace_back();
  Interval line;
  for (auto &sorted_point : sorted_points) {
    auto point_type = get_point_type(sorted_point, point_to_edges);
    switch (point_type) {
      case PointType::NORMAL:
        // change edge adjusted to current sweeping line
        _process_normal_point(sorted_point, edges, point_to_edges,
                              segment_to_line);

        break;
      case PointType::SPLIT:
        // split sweeping line on two lines
        // add edge sor cutting polygon on two parts
        _process_split_point(sorted_point, edges, point_to_edges,
                             segment_to_line);
        break;
      case PointType::MERGE:
        // merge two sweeping lines to one
        // save point as point to start new line for SPLIT point case
        _process_merge_point(sorted_point, edges, point_to_edges,
                             segment_to_line);
        break;
      case PointType::INTERSECTION:
        // this is merge and split point at same time
        // this is not described in original algorithm
        // but we need it to handle self intersecting polygons
        // Remember about more than 4 edges case
        break;
    }
  }
  return std::make_pair(result, polygon);
}

std::pair<std::vector<Triangle>, std::vector<point::Point>>
_triangulate_polygon(const std::vector<point::Point> &polygon) {
  if (polygon.size() < 3)
    return std::make_pair(std::vector<Triangle>(), polygon);
  if (polygon.size() == 3)
    return std::make_pair(std::vector<Triangle>({Triangle(0, 1, 2)}), polygon);
  if (polygon.size() == 4) {
    if (partsegcore::intersection::_orientation(polygon[0], polygon[1],
                                                polygon[2]) !=
        partsegcore::intersection::_orientation(polygon[0], polygon[3],
                                                polygon[2]))
      return std::make_pair(
          std::vector<Triangle>({Triangle(0, 1, 2), Triangle(0, 3, 2)}),
          polygon);
  }

  if (_is_convex(polygon))
    return std::make_pair(_triangle_convex_polygon(polygon), polygon);

  // Implement sweeping line algorithm for triangulation
  // described on this lecture:
  // https://www.youtube.com/playlist?list=PLtTatrCwXHzEqzJMaTUFgqoCNllgwk4DH
  //
  return sweeping_line_triangulation(find_intersection_points(polygon));
}
}  // namespace triangulation
}  // namespace partsegcore

#endif  // PARTSEGCORE_TRIANGULATE_H
