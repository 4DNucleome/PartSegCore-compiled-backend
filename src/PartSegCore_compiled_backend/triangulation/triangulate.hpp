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

enum PointType { NORMAL, SPLIT, MERGE, INTERSECTION, EMPTY };

struct MonotonePolygon {
  point::Point top{};
  point::Point bottom{};
  std::vector<point::Point> left;
  std::vector<point::Point> right;

  MonotonePolygon() = default;
  explicit MonotonePolygon(point::Point top) : top(top) {}
  MonotonePolygon(point::Point top, point::Point bottom,
                  std::vector<point::Point> left,
                  std::vector<point::Point> right)
      : top(top),
        bottom(bottom),
        left(std::move(left)),
        right(std::move(right)) {}
};

struct Interval {
  point::Point last_seen{};
  point::Segment left_segment;
  point::Segment right_segment;
  std::vector<MonotonePolygon *> polygons_list{};
  Interval() = default;
  explicit Interval(const point::Point &p, const point::Segment &left,
                    const point::Segment &right)
      : last_seen(p), left_segment(left), right_segment(right) {};
  explicit Interval(const point::Point &p, const point::Segment &left,
                    const point::Segment &right, MonotonePolygon *polygon)
      : last_seen(p),
        left_segment(left),
        right_segment(right),
        polygons_list({polygon}) {};

  void replace_segment(const point::Segment &old_segment,
                       const point::Segment &new_segment) {
    if (left_segment == old_segment) {
      left_segment = new_segment;
      return;
    } else if (right_segment == old_segment) {
      right_segment = new_segment;
      return;
    }
    throw std::runtime_error("Segment not found in interval");
  };
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

struct PointMonotonePolygon {
  MonotonePolygon *left_polygon = nullptr;
  MonotonePolygon *right_polygon = nullptr;
  PointMonotonePolygon() = default;
  explicit PointMonotonePolygon(MonotonePolygon *left_polygon,
                                MonotonePolygon *right_polygon)
      : left_polygon(left_polygon), right_polygon(right_polygon) {};
};

struct Triangle {
  std::size_t x;
  std::size_t y;
  std::size_t z;
  Triangle(std::size_t x, std::size_t y, std::size_t z) : x(x), y(y), z(z) {};
  Triangle() = default;
};

struct PointTriangle {
  point::Point p1{};
  point::Point p2{};
  point::Point p3{};

  PointTriangle(point::Point p1, point::Point p2, point::Point p3) {
    if ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x) < 0) {
      this->p1 = p3;
      this->p2 = p2;
      this->p3 = p1;
    } else {
      this->p1 = p1;
      this->p2 = p2;
      this->p3 = p3;
    }
  };
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
 * Get a mapping from points to the list of edges that contain those points.
 * Each list of edges is sorted by the point order within the edges.
 *
 * @param edges A vector of Segment objects representing the edges of a polygon.
 * @return A PointToEdges map where each key is a point and the corresponding
 * value is a vector of edges containing that point.
 */
PointToEdges get_points_edges(const std::vector<point::Segment> &edges) {
  PointToEdges point_to_edges;
  for (std::size_t i = 0; i < edges.size(); i++) {
    point_to_edges[edges[i].bottom].emplace_back(i, edges[i].top);
    point_to_edges[edges[i].top].emplace_back(i, edges[i].bottom);
  }
  for (auto &point_to_edge : point_to_edges) {
    std::sort(point_to_edge.second.begin(), point_to_edge.second.end(),
              [](const PointEdges &a, const PointEdges &b) {
                return b.opposite_point < a.opposite_point;
              });
  }
  return point_to_edges;
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
  if (point_to_edges.at(p).empty()) return PointType::EMPTY;
  if (point_to_edges.at(p).size() != 2) return PointType::INTERSECTION;

  const auto &edges = point_to_edges.at(p);
  if (edges[0].opposite_point < p && edges[1].opposite_point < p)
    return PointType::SPLIT;
  if (p < edges[0].opposite_point && p < edges[1].opposite_point)
    return PointType::MERGE;
  return PointType::NORMAL;
}

struct MonotonePolygonBuilder {
  std::map<point::Segment, Interval *> segment_to_line{};
  std::vector<point::Segment> edges{};
  PointToEdges point_to_edges{};
  std::vector<MonotonePolygon> monotone_polygons{};

  MonotonePolygonBuilder() = default;
  explicit MonotonePolygonBuilder(const std::vector<point::Segment> &edges)
      : edges(edges) {
    this->point_to_edges = get_points_edges(edges);
  }

  std::pair<const point::Segment &, const point::Segment &>
  get_left_right_edges(const point::Point &p) {
    auto point_info = point_to_edges.at(p);
    auto fst_idx = point_info[0].edge_index;
    auto snd_idx = point_info[1].edge_index;
    if (edges[fst_idx].top == edges[snd_idx].top) {
      if (intersection::_orientation(edges[fst_idx].bottom, edges[fst_idx].top,
                                     edges[snd_idx].bottom) == 2) {
        return {edges[snd_idx], edges[fst_idx]};
      }
    }
    if (intersection::_orientation(edges[fst_idx].top, edges[fst_idx].bottom,
                                   edges[snd_idx].top) == 1) {
      return {edges[snd_idx], edges[fst_idx]};
    }

    return {edges[fst_idx], edges[snd_idx]};
  }

  /**
   * Processes the end point of a segment within a monotone polygon.
   *
   * @param p The end point to be processed.
   *
   * The function retrieves two segments (edges) associated with the
   * given point and checks the monotone polygon information of the
   * top points of these segments. Depending on whether the end point
   * belongs to one or two monotone polygons, it updates the monotone
   * polygon information and stores completed polygons.
   *
   * It also handles the cleanup of relevant mappings and data structures
   * when processing is complete.
   */

  void process_end_point(const point::Point &p) {
    auto [edge_left, edge_right] = get_left_right_edges(p);

    Interval *interval = segment_to_line.at(edge_left);

    for (auto &polygon : interval->polygons_list) {
      polygon->bottom = p;
      monotone_polygons.push_back(*polygon);
      delete polygon;
    }

    segment_to_line.erase(edge_left);
    segment_to_line.erase(edge_right);
    delete interval;
  }

  void process_merge_point(const point::Point &p) {
    auto [edge_left, edge_right] = get_left_right_edges(p);

    if (segment_to_line.at(edge_left) != segment_to_line.at(edge_right)) {
      // merge two intervals into one
      Interval *left_interval = segment_to_line.at(edge_left);
      Interval *right_interval = segment_to_line.at(edge_right);
      segment_to_line.erase(edge_right);
      left_interval->right_segment = right_interval->right_segment;
      segment_to_line[right_interval->right_segment] = left_interval;
      left_interval->last_seen = p;
      left_interval->polygons_list.back()->right.push_back(p);
      right_interval->polygons_list.front()->left.push_back(p);
      for (auto &polygon : right_interval->polygons_list) {
        left_interval->polygons_list.push_back(polygon);
      }
      delete right_interval;
    } else {
      // This is the end point
      this->process_end_point(p);
    }
  };

  /**
   * Processes a normal point within a monotone polygon.
   *
   * @param p The normal point to be processed.
   *
   * The function identifies two edges associated with the given point and
   * retrieves the corresponding segments. It then finds the interval and
   * associated monotone polygons based on the edges. If the last seen point in
   * the interval is a merge point, it updates and stores the completed monotone
   * polygon. The function also updates the right and left polygons and handles
   * segment replacements along with interval mappings, ensuring the data
   * structures are correctly updated.
   */
  void process_normal_point(const point::Point &p) {
    const point::Segment &edge_top =
        edges[point_to_edges.at(p).at(0).edge_index];
    const point::Segment &edge_bottom =
        edges[point_to_edges.at(p).at(1).edge_index];
    if (segment_to_line.count(edge_top) == 0) {
      throw std::runtime_error("Segment not found in the map");
    }
    Interval *interval = segment_to_line.at(edge_top);

    if (interval->polygons_list.size() > 1) {
      // end all the polygons, except 1
      if (edge_top == interval->right_segment) {
        for (auto i = 1; i < interval->polygons_list.size(); i++) {
          interval->polygons_list[i]->bottom = p;
          monotone_polygons.push_back(*interval->polygons_list[i]);
          delete interval->polygons_list[i];
        }

      } else {
        for (auto i = 0; i < interval->polygons_list.size() - 1; i++) {
          interval->polygons_list[i]->bottom = p;
          monotone_polygons.push_back(*interval->polygons_list[i]);
          delete interval->polygons_list[i];
        }
        interval->polygons_list[0] = interval->polygons_list.back();
      }
      interval->polygons_list.erase(interval->polygons_list.begin() + 1,
                                    interval->polygons_list.end());
    }

    if (edge_top == interval->right_segment) {
      interval->polygons_list[0]->right.push_back(p);
    } else {
      interval->polygons_list[0]->left.push_back(p);
    }

    segment_to_line[edge_bottom] = interval;
    interval->last_seen = p;
    interval->replace_segment(edge_top, edge_bottom);
    segment_to_line.erase(edge_top);
  };

  void process_start_point(const point::Point &p) {
    auto [edge_left, edge_right] = get_left_right_edges(p);

    auto *new_polygon = new MonotonePolygon(p);
    auto *interval = new Interval(p, edge_left, edge_right, new_polygon);
    segment_to_line[edge_left] = interval;
    segment_to_line[edge_right] = interval;
  }

  void process_split_point(const point::Point &p) {
    auto [edge_left, edge_right] = get_left_right_edges(p);

    // We need to find to which interval the point belongs.
    // If the point does not belong to any interval, we treat
    // it as a start point and create a new interval

    for (auto &segment_interval : segment_to_line) {
      if (segment_interval.first == segment_interval.second->right_segment) {
        // scan interval only once
        continue;
      }
      // check if the current point is inside the quadrangle defined by edges of
      // the interval
      Interval *interval = segment_interval.second;
      if (interval->left_segment.point_on_line(p.y) < p.x &&
          interval->right_segment.point_on_line(p.y) > p.x) {
        // the point is inside the interval

        // update the sweep line
        auto right_segment = interval->right_segment;
        interval->right_segment = edge_left;
        interval->last_seen = p;
        segment_to_line[edge_left] = interval;

        auto *new_interval = new Interval(p, edge_right, right_segment);
        segment_to_line[edge_right] = new_interval;
        segment_to_line[right_segment] = new_interval;

        if (interval->polygons_list.size() == 1) {
          MonotonePolygon *new_polygon = nullptr;
          if (interval->polygons_list[0]->right.empty()) {
            new_polygon = new MonotonePolygon(interval->polygons_list[0]->top);
          } else {
            new_polygon =
                new MonotonePolygon(interval->polygons_list[0]->right.back());
          }
          new_polygon->left.push_back(p);
          new_interval->polygons_list.push_back(new_polygon);
          interval->polygons_list[0]->right.push_back(p);
        }

        if (interval->polygons_list.size() >= 2) {
          interval->polygons_list[0]->right.push_back(p);
          interval->polygons_list[interval->polygons_list.size() - 1]
              ->left.push_back(p);
          for (auto i = 1; i < interval->polygons_list.size() - 1; i++) {
            interval->polygons_list[i]->bottom = p;
            monotone_polygons.push_back(*interval->polygons_list[i]);
            delete interval->polygons_list[i];
          }
          new_interval->polygons_list.push_back(interval->polygons_list.back());
          interval->polygons_list.erase(interval->polygons_list.begin() + 1,
                                        interval->polygons_list.end());
        }
        return;
      }
    }
    // the point is not inside any interval
    this->process_start_point(p);
  };

  void process_intersection_point(const point::Point &p) {
    throw std::runtime_error("Intersection points are not supported");
  };
};

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
  // remove all elements except two firsts
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
 * - TOP_OR_BOTTOM: Represents the top side.
 * - LEFT: Represents the left side.
 * - RIGHT: Represents the right side.
 */
enum Side {
  TOP_OR_BOTTOM = 0,
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

  result.reserve(polygon.left.size() + polygon.right.size());

  points.reserve(polygon.left.size() + polygon.right.size() + 2);
  points.emplace_back(polygon.top, Side::TOP_OR_BOTTOM);
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
  points.emplace_back(polygon.bottom, Side::TOP_OR_BOTTOM);

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
 * Calculates the edges of polygons from a list of polygons, provided as
 * a list of points for each polygon.
 *
 * Each polygon is represented by a vector of points, where each point is
 * defined as an instance of `point::Point`. The function iterates through
 * each point in each polygon and creates edges (segments) between
 * consecutive points as well as between the last point and the first point
 * of each polygon.
 *
 * @param polygon_list A vector of polygons, with each polygon represented
 *                     as a vector of `point::Point` instances.
 *
 * @return A vector of `point::Segment` instances representing the edges of
 *         all polygons in the input list.
 */
std::vector<point::Segment> calc_edges(
    const std::vector<std::vector<point::Point>> &polygon_list) {
  std::vector<point::Segment> edges;
  std::size_t points_count = 0;
  for (const auto &polygon : polygon_list) {
    points_count += polygon.size();
  }
  edges.reserve(points_count);
  for (const auto &polygon : polygon_list) {
    for (std::size_t i = 0; i < polygon.size() - 1; i++) {
      edges.emplace_back(polygon[i], polygon[i + 1]);
    }
    edges.emplace_back(polygon[polygon.size() - 1], polygon[0]);
  }
  return edges;
}

/**
 * Calculates and returns a list of unique (deduplicated) edges from a list of
 * polygons.
 *
 * This function receives a list of polygons, where each polygon is defined by a
 * series of points. It identifies the edges of the polygons and returns those
 * edges that are unique (i.e., edges that appear an odd number of times across
 * all polygons). Any edge that appears an even number of times is considered to
 * be a duplicate and is removed.
 *
 * @param polygon_list A reference to a vector containing vectors of points,
 * where each inner vector represents a polygon.
 * @return A vector of unique edges (of type point::Segment) that are not
 * duplicated across the polygons.
 */
std::vector<point::Segment> calc_dedup_edges(
    const std::vector<std::vector<point::Point>> &polygon_list) {
  std::set<point::Segment> edges_set;
  point::Segment edge;

  for (const auto &polygon : polygon_list) {
    for (std::size_t i = 0; i < polygon.size() - 1; i++) {
      edge = point::Segment(polygon[i], polygon[i + 1]);
      if (edges_set.count(edge) == 0) {
        edges_set.insert(edge);
      } else {
        edges_set.erase(edge);
      }
    }
    edge = point::Segment(polygon[polygon.size() - 1], polygon[0]);
    if (edges_set.count(edge) == 0) {
      edges_set.insert(edge);
    } else {
      edges_set.erase(edge);
    }
  }
  return {edges_set.begin(), edges_set.end()};
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
std::vector<std::vector<point::Point>> find_intersection_points(
    const std::vector<std::vector<point::Point>> &polygon_list) {
  /* find all edge intersections and add mid-points for all such intersection
   * places*/
  auto edges = calc_edges(polygon_list);

  auto intersections = intersection::_find_intersections(edges);
  if (intersections.empty()) return polygon_list;
  std::unordered_map<std::size_t, std::vector<point::Point>>
      intersections_points;
  for (const auto &intersection : intersections) {
    auto inter_point = intersection::_find_intersection(
        edges[intersection.first], edges[intersection.second]);
    intersections_points[intersection.first].push_back(inter_point);
    intersections_points[intersection.second].push_back(inter_point);
  };
  for (auto &intersections_point : intersections_points) {
    //    points_count += intersections_point.second.size() - 1;
    intersections_point.second.push_back(edges[intersections_point.first].top);
    intersections_point.second.push_back(
        edges[intersections_point.first].bottom);
    std::sort(intersections_point.second.begin(),
              intersections_point.second.end());
  }

  std::vector<std::vector<point::Point>> new_polygons_list;

  for (const auto &polygon : polygon_list) {
    std::vector<point::Point> new_polygon;
    new_polygon.reserve(polygon.size() * 2);
    new_polygon.push_back(polygon[0]);
    for (std::size_t i = 0; i < polygon.size(); i++) {
      auto point = polygon[i];
      if (new_polygon[new_polygon.size() - 1] != point)
        new_polygon.push_back(point);
      if (intersections_points.count(i)) {
        auto &new_points = intersections_points[i];
        if (new_points[0] == point) {
          for (auto it = new_points.begin() + 1; it != new_points.end(); it++) {
            if (new_polygon.back() != *it) new_polygon.push_back(*it);
          }
        } else {
          for (auto it = new_points.rbegin() + 1; it != new_points.rend();
               it++) {
            if (new_polygon.back() != *it) new_polygon.push_back(*it);
          }
        }
      }
    }
    new_polygons_list.push_back(new_polygon);
  }
  return new_polygons_list;
}

std::vector<point::Point> find_intersection_points(
    const std::vector<point::Point> &polygon) {
  auto new_polygon = find_intersection_points(
      std::vector<std::vector<point::Point>>({polygon}));
  return new_polygon[0];
}

std::vector<point::Point> _sorted_polygons_points(
    const std::vector<std::vector<point::Point>> &polygon_list) {
  std::vector<point::Point> result;
  std::unordered_set<point::Point> visited;
  for (const auto &polygon : polygon_list) {
    for (const auto &point : polygon) {
      if (visited.count(point) == 0) {
        result.push_back(point);
        visited.insert(point);
      }
    }
  }
  std::sort(
      result.begin(), result.end(),
      [](const point::Point &p1, const point::Point &p2) { return p2 < p1; });
  return result;
}

/*
    This is an implementation of sweeping line triangulation of polygon
    Its assumes that there is no edge intersections, but may be a point with
    more than 2 edges. described on this lecture:
    https://www.youtube.com/playlist?list=PLtTatrCwXHzEqzJMaTUFgqoCNllgwk4DH
    */
std::pair<std::vector<Triangle>, std::vector<point::Point>>
sweeping_line_triangulation(
    const std::vector<std::vector<point::Point>> &polygon_list) {
  std::vector<Triangle> result;
  auto edges = calc_dedup_edges(polygon_list);
  MonotonePolygonBuilder builder(edges);

  PointToEdges point_to_edges = get_points_edges(edges);

  std::vector<point::Point> sorted_points =
      _sorted_polygons_points(polygon_list);

  result.reserve(sorted_points.size() - 2);

  for (auto &sorted_point : sorted_points) {
    auto point_type = get_point_type(sorted_point, point_to_edges);
    switch (point_type) {
      case PointType::NORMAL:
        // Change edge adjusted to the current sweeping line.
        // Adds edge to monotone polygon.
        builder.process_normal_point(sorted_point);
        break;
      case PointType::SPLIT:
        // Split a sweeping line on two lines,
        // adds edge sor cutting polygon on two parts.
        builder.process_split_point(sorted_point);
        break;
      case PointType::MERGE:
        // merge two sweeping lines to one
        // merge two intervals into one
        // It may be the end of the polygon, then finish the
        // monotone polygon
        builder.process_merge_point(sorted_point);
        break;
      case PointType::INTERSECTION:
        // this is a merge and split point at the same time
        // this is not described in original algorithm
        // but we need it to handle self intersecting polygons
        // Remember about more than 4 edges case
        builder.process_intersection_point(sorted_point);
        break;
      case PointType::EMPTY:
        // this is a point without edges (removed by deduplication)
        break;
    }
  }
  std::unordered_map<point::Point, std::size_t> point_to_index;
  point_to_index.reserve(sorted_points.size());
  for (std::size_t i = 0; i < sorted_points.size(); i++) {
    point_to_index[sorted_points[i]] = i;
  }
  for (auto &monotone_polygon : builder.monotone_polygons) {
    auto triangles = triangulate_monotone_polygon(monotone_polygon);
    for (auto &triangle : triangles) {
      result.emplace_back(point_to_index[triangle.p1],
                          point_to_index[triangle.p2],
                          point_to_index[triangle.p3]);
    }
  }
  return std::make_pair(result, sorted_points);
}

// calculate the triangulation of a symmetric difference of list of polygons

std::pair<std::vector<Triangle>, std::vector<point::Point>> triangulate_polygon(
    const std::vector<std::vector<point::Point>> &polygon_list) {
  if (polygon_list.empty())
    // empty list
    return std::make_pair(std::vector<Triangle>(), std::vector<point::Point>());
  if (polygon_list.size() == 1) {
    // only one polygon in the list
    std::vector<point::Point> polygon = polygon_list[0];
    if (polygon.size() < 3)
      return std::make_pair(std::vector<Triangle>(), polygon);
    if (polygon.size() == 3)
      return std::make_pair(std::vector<Triangle>({Triangle(0, 1, 2)}),
                            polygon);
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
  }

  // Implement the sweeping line algorithm for triangulation
  // described on this lecture:
  // https://www.youtube.com/playlist?list=PLtTatrCwXHzEqzJMaTUFgqoCNllgwk4DH
  //
  return sweeping_line_triangulation(find_intersection_points(polygon_list));
}

std::pair<std::vector<Triangle>, std::vector<point::Point>> triangulate_polygon(
    const std::vector<point::Point> &polygon_list) {
  return triangulate_polygon(
      std::vector<std::vector<point::Point>>({polygon_list}));
}

}  // namespace triangulation
}  // namespace partsegcore

#endif  // PARTSEGCORE_TRIANGULATE_H
