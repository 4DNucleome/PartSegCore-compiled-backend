#ifndef PARTSEGCORE_TRIANGULATE_H
#define PARTSEGCORE_TRIANGULATE_H

#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "intersection.hpp"
#include "point.hpp"

namespace partsegcore {
namespace triangulation {

enum PointType { NORMAL, SPLIT, MERGE, INTERSECTION };

struct Line {
  point::Segment left, right;
  Line() = default;
  Line(const point::Segment &left, const point::Segment &right)
      : left(left), right(right) {};
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

typedef std::unordered_map<point::Point, std::vector<PointEdges>,
                           point::PointHash>
    PointToEdges;

bool cmp_point_edges(const PointEdges &p, const PointEdges &q) { return p < q; }

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
              intersections_point.second.end(), point::cmp_point);
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
    std::sort(point_to_edge.second.begin(), point_to_edge.second.end(),
              cmp_point_edges);
  }
  return point_to_edges;
}

/*
    This is implementation of sweeping line triangulation of polygon
    Its assumes that there is no edge intersections, but may be a point with
   more than 2 edges. described on this lecture:
    https://www.youtube.com/playlist?list=PLtTatrCwXHzEqzJMaTUFgqoCNllgwk4DH
    */
std::pair<std::vector<Triangle>, std::vector<point::Point>>
sweeping_line_triangulation(const std::vector<point::Point> &polygon) {
  std::vector<Triangle> result;
  auto edges = calc_edges(polygon);
  PointToEdges point_to_edges = get_points_edges(edges);

  std::vector<point::Point> sorted_points = polygon;
  // copy to avoid modification of original vector
  std::sort(sorted_points.begin(), sorted_points.end(), point::cmp_point);
  std::vector<OrderedPolygon> ordered_polygon_li;
  ordered_polygon_li.emplace_back();
  Line line;
  for (auto &sorted_point : sorted_points) {
    auto point_type = get_point_type(sorted_point, point_to_edges);
    switch (point_type) {
      case PointType::NORMAL:
        // change edge adjusted to current sweeping line
        break;
      case PointType::SPLIT:
        // split sweeping line on two lines
        // add edge sor cutting polygon on two parts
        line = Line(
            point::Segment(sorted_point,
                           point_to_edges.at(sorted_point)[0].opposite_point),
            point::Segment(sorted_point,
                           point_to_edges.at(sorted_point)[1].opposite_point));
        break;
      case PointType::MERGE:
        // merge two sweeping lines to one
        // save point as point to start new line for SPLIT pointcase
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
