#ifndef PARTSEGCORE_TRIANGULATE_H
#define PARTSEGCORE_TRIANGULATE_H

#include <algorithm>
#include <map>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "debug_util.hpp"
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

/**
 * Comparator function to determine the relative positioning of two given
 * segments.
 *
 * Compares two segments to identify if the first segment (`s1`) is to the left
 * and below the second segment (`s2`) when moving left to right. The function
 * assumes that the segments do not intersect.
 *
 * TODO Investigate if it could be implemented in a more efficient way.
 *
 * @param s1 The first segment to compare.
 * @param s2 The second segment to compare.
 *
 * @return `true` if `s1` is to the left and below `s2`, `false` otherwise.
 */
bool left_to_right(const point::Segment &s1, const point::Segment &s2) {
  if (s1.is_horizontal() && s2.is_horizontal()) {
    return s1.bottom.x < s2.bottom.x;
  }
  if (s1.is_horizontal()) {
    int i1 = intersection::_orientation(s2.bottom, s2.top, s1.top);
    int i2 = intersection::_orientation(s2.bottom, s2.top, s1.bottom);
    return (i1 == 2 || i2 == 2);
  }
  if (s2.is_horizontal()) {
    int i1 = intersection::_orientation(s1.bottom, s1.top, s2.top);
    int i2 = intersection::_orientation(s1.bottom, s1.top, s2.bottom);
    return (i1 == 1 || i2 == 1);
  }

  if ((s1.top.y < s2.top.y)) {
    if (s1.top == s2.top) {
      return intersection::_orientation(s2.top, s2.bottom, s1.bottom) == 1;
    }
    return intersection::_orientation(s2.top, s2.bottom, s1.top) == 1;
  }
  if (s1.top == s2.top) {
    return intersection::_orientation(s1.top, s1.bottom, s2.bottom) == 2;
  }
  return intersection::_orientation(s1.top, s1.bottom, s2.top) == 2;
}

struct SegmentLeftRightComparator {
  bool operator()(const point::Segment &s1, const point::Segment &s2) const {
    return left_to_right(s1, s2);
  }
};

struct MonotonePolygon {
  point::Point top{};
  point::Point bottom{};
  std::vector<point::Point> left;
  std::vector<point::Point> right;

  MonotonePolygon() = default;
  MonotonePolygon(point::Point top, point::Point bottom,
                  std::vector<point::Point> left,
                  std::vector<point::Point> right)
      : top(top),
        bottom(bottom),
        left(std::move(left)),
        right(std::move(right)) {}
};

struct Triangle {
  std::size_t x;
  std::size_t y;
  std::size_t z;
  Triangle(std::size_t x, std::size_t y, std::size_t z) : x(x), y(y), z(z) {};
  Triangle() = default;
};

struct PointTriangle {
  point::Point p1;
  point::Point p2;
  point::Point p3;
  PointTriangle(point::Point p1, point::Point p2, point::Point p3)
      : p1(p1), p2(p2), p3(p3) {};
  PointTriangle() = default;
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

/**
 * Construct triangles using the current point and edges opposite to it.
 *
 * This function takes a current point and a stack of points, then builds
 * triangles by connecting the current point to adjoining points in the stack.
 * Each triangle is added to the result vector as a `PointTriangle` object.
 *
 * Steps performed:
 * 1. Iterate over the stack to form triangles with the current point.
 * 2. Replace the first element of the stack with the last element.
 * 3. Set the second element of the stack to the current point.
 * 4. Remove all elements from the stack except the first two.
 *
 * @param stack A vector of `point::Point` objects representing the stack of
 * points.
 * @param result A vector to store the resultant `PointTriangle` objects.
 * @param current_point The current `point::Point` used to form triangles with
 * points in the stack.
 */
void _build_triangles_opposite_edge(std::vector<point::Point> &stack,
                                    std::vector<PointTriangle> &result,
                                    point::Point current_point) {
  for (std::size_t i = 0; i < stack.size() - 1; i++) {
    result.emplace_back(current_point, stack[i], stack[i + 1]);
  }
  stack[0] = stack.back();
  stack[1] = current_point;
  // remove all elements except two first
  stack.erase(stack.begin() + 2, stack.end());
}

/**
 * Constructs triangles from the current edge of the y-monotone and updates the
 * stack and result vectors.
 *
 * This function processes a stack of points and an incoming point to form
 * triangles that are part of the processed area. The triangles are determined
 * based on the expected orientation and are added to the result vector. The
 * stack is then updated to reflect the current status of the processed area.
 *
 * @param point_stack A reference to a std::vector of point::Point representing
 * the current state of the stack.
 * @param triangles A reference to a std::vector of PointTriangle where the
 * generated triangles will be stored.
 * @param incoming_point The incoming point::Point to be considered for new
 * triangles and updating the stack.
 * @param expected_orientation The expected orientation (clockwise or
 * counterclockwise) that will guide the triangle formation.
 */
void _build_triangles_current_edge(std::vector<point::Point> &stack,
                                   std::vector<PointTriangle> &result,
                                   point::Point current_point,
                                   int expected_orientation) {
  auto it1 = stack.rbegin();
  auto it2 = stack.rbegin() + 1;
  while (it2 != stack.rend() &&
         intersection::_orientation(*it2, *it1, current_point) ==
             expected_orientation) {
    result.emplace_back(current_point, *it1, *it2);
    it1++;
    it2++;
  }
  stack.erase(it1.base(), stack.end());
  stack.push_back(current_point);
}

/**
 * @enum Side
 * @brief An enumeration to represent different sides.
 *
 * This enumeration defines constants for various sides,
 * in particular:
 * - TOP: Represents the top side.
 * - LEFT: Represents the left side.
 * - RIGHT: Represents the right side.
 */
enum Side {
  TOP = 0,
  LEFT = 2,
  RIGHT = 1,
};

/**
 * Triangulates a given monotone polygon.
 *
 * The function takes a monotone polygon as input and returns a vector
 * of triangles that represent the triangulation of the polygon. The
 * polygon must be y-monotone, meaning that every line segment
 * parallel to the y-axis intersects the polygon at most twice.
 *
 * The algorithm works by sorting points of the polygon, merging the
 * left and right chains, and then using a stack-based approach to
 * recursively build triangles.
 *
 * @param polygon The monotone polygon to be triangulated.
 * @return A vector of PointTriangle representing the triangles of the
 *         triangulated polygon.
 */
std::vector<PointTriangle> triangulate_monotone_polygon(
    const MonotonePolygon &polygon) {
  std::vector<PointTriangle> result;
  std::size_t left_index = 0;
  std::size_t right_index = 0;
  std::vector<point::Point> stack;
  std::vector<std::pair<point::Point, Side>> points;

  points.reserve(polygon.left.size() + polygon.right.size() + 2);
  points.emplace_back(polygon.top, Side::TOP);
  while (left_index < polygon.left.size() &&
         right_index < polygon.right.size()) {
    if (polygon.left[left_index] < polygon.right[right_index]) {
      points.emplace_back(polygon.right[right_index], Side::RIGHT);
      right_index++;
    } else {
      points.emplace_back(polygon.left[left_index], Side::LEFT);
      left_index++;
    }
  }
  while (left_index < polygon.left.size()) {
    points.emplace_back(polygon.left[left_index], Side::LEFT);
    left_index++;
  }
  while (right_index < polygon.right.size()) {
    points.emplace_back(polygon.right[right_index], Side::RIGHT);
    right_index++;
  }
  points.emplace_back(polygon.bottom, Side::TOP);

  stack.push_back(points[0].first);
  stack.push_back(points[1].first);
  Side side = points[1].second;

  for (std::size_t i = 2; i < points.size(); i++) {
    if (side == points[i].second) {
      _build_triangles_current_edge(stack, result, points[i].first, side);
    } else {
      _build_triangles_opposite_edge(stack, result, points[i].first);
    }
    side = points[i].second;
  }
  return result;
}

/**
 * Calculates the edges of a polygon represented as a sequence of points.
 *
 * This function takes a vector of points representing the vertices of a polygon
 * and returns a vector of segments representing the edges of the polygon.
 * The edges are created by connecting each point to the next, with the final
 * point connecting back to the first point to close the polygon.
 *
 * @param polygon A vector of points representing the vertices of the polygon.
 * @return A vector of segments representing the edges of the polygon.
 */
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

/**
 * @brief Finds intersection points in a polygon and adds mid-points for all
 * intersections.
 *
 * The function takes a vector of points defining a polygon and finds all edge
 * intersections. It then adds mid-points for all such intersections.
 *
 * @param polygon The polygon defined by a vector of Point objects.
 * @return A new vector of Point objects representing the polygon with added
 * intersection points.
 */
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
    intersections_point.second.push_back(edges[intersections_point.first].top);
    intersections_point.second.push_back(
        edges[intersections_point.first].bottom);
    std::sort(intersections_point.second.begin(),
              intersections_point.second.end());
  }

  std::vector<point::Point> new_polygon;
  new_polygon.reserve(points_count);
  for (std::size_t i = 0; i < polygon.size(); i++) {
    auto point = polygon[i];
    if (new_polygon[new_polygon.size() - 1] != point)
      new_polygon.push_back(point);
    if (intersections_points.count(i)) {
      auto new_points = intersections_points[i];
      if (new_points[0] == point) {
        for (std::size_t j = 1; j < new_points.size() - 1; j++) {
          if (new_polygon[new_polygon.size() - 1] != new_points[j])
            new_polygon.push_back(new_points[j]);
        }
      } else {
        for (std::size_t j = new_points.size() - 2; j > 0; j++) {
          if (new_polygon[new_polygon.size() - 1] != new_points[j])
            new_polygon.push_back(new_points[j]);
        }
      }
    }
  }
  return new_polygon;
}

/**
 * Determines the type of a given point based on its adjacent edges.
 *
 * If the point has more than two edges connected to it, it is considered an
 * intersection point. If it has exactly two edges, further checks categorize
 * it as one of split, merge, or normal points.
 *
 * - If both adjacent edges have their opposite ends before the given point `p`,
 *   the point is classified as a merge point.
 * - If both adjacent edges have their opposite ends after the given point `p`,
 *   the point is classified as a split point.
 * - Otherwise, it is categorized as a normal point.
 *
 * @param p The point to classify.
 * @param point_to_edges A mapping from points to their adjacent edges.
 * @return The type of the point as determined by its adjacent edges.
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

/**
 * Get a mapping from points to the list of edges that contain those points.
 * Each list of edges is sorted by the point order within the edges.
 *
 * @param edges A vector of Segment objects representing the edges of a polygon.
 * @return A PointToEdges map where each key is a point and the corresponding
 * value is a vector of edges containing that point.
 */
PointToEdges get_points_edges(std::vector<point::Segment> &edges) {
  PointToEdges point_to_edges;
  for (std::size_t i = 0; i < edges.size(); i++) {
    point_to_edges[edges[i].bottom].emplace_back(i, edges[i].top);
    point_to_edges[edges[i].top].emplace_back(i, edges[i].bottom);
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
    This is an implementation of sweeping line triangulation of polygon
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
  std::vector<MonotonePolygon> ordered_polygon_li;
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
