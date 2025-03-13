#ifndef PARTSEGCORE_TRIANGULATE_H
#define PARTSEGCORE_TRIANGULATE_H

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>  // memory header is required on linux, and not on macos
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "intersection.hpp"
#include "point.hpp"

namespace partsegcore::triangulation {

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
  std::vector<std::unique_ptr<MonotonePolygon>> polygons_list{};
  Interval() = default;
  explicit Interval(const point::Point &p, const point::Segment &left,
                    const point::Segment &right)
      : last_seen(p), left_segment(left), right_segment(right) {};
  explicit Interval(const point::Point &p, const point::Segment &left,
                    const point::Segment &right,
                    std::unique_ptr<MonotonePolygon> &polygon)
      : last_seen(p),
        left_segment(left),
        right_segment(right),
        polygons_list() {
    polygons_list.emplace_back(std::move(polygon));
  };

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
  [[nodiscard]] point::Segment opposite_segment(
      const point::Segment &segment) const {
    if (segment == left_segment) {
      return right_segment;
    }
    if (segment == right_segment) {
      return left_segment;
    }
    throw std::runtime_error("Segment not found in interval");
  };

  // ostream operator for Interval
  friend std::ostream &operator<<(std::ostream &os, const Interval &interval) {
    os << "Last Seen: " << interval.last_seen
       << ", Left Segment: " << interval.left_segment
       << ", Right Segment: " << interval.right_segment
       << ", Polygons count: " << interval.polygons_list.size();
    return os;
  }
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
inline bool left_to_right(const point::Segment &s1, const point::Segment &s2) {
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

inline bool left_right_share_top(const point::Segment &s1,
                                 const point::Segment &s2) {
  return intersection::_orientation(s1.bottom, s1.top, s2.bottom) == 1;
}

inline bool left_right_share_bottom(const point::Segment &s1,
                                    const point::Segment &s2) {
  return intersection::_orientation(s1.top, s1.bottom, s2.top) == 2;
}

struct SegmentLeftRightComparator {
  bool operator()(const point::Segment &s1, const point::Segment &s2) const {
    return left_to_right(s1, s2);
  }
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
    /* Sort the points in clockwise order */
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

struct PointEdgeInfo {
  EdgeIndex edge_index;
  point::Point opposite_point;
  PointEdgeInfo(EdgeIndex edge_index, point::Point opposite_point)
      : edge_index(edge_index), opposite_point(opposite_point) {}
  bool operator<(const PointEdgeInfo &e) const {
    return opposite_point < e.opposite_point;
  }
};

typedef std::unordered_map<point::Point, std::vector<PointEdgeInfo>>
    PointToEdges;

/**
 * Get a mapping from points to the list of edges that contain those points.
 * Each list of edges is sorted by the point order within the edges.
 *
 * @param edges A vector of Segment objects representing the edges of a polygon.
 * @return A PointToEdges map where each key is a point and the corresponding
 * value is a vector of edges containing that point.
 */
inline PointToEdges get_points_edges(const std::vector<point::Segment> &edges) {
  PointToEdges point_to_edges;
  for (std::size_t i = 0; i < edges.size(); i++) {
    point_to_edges[edges[i].bottom].emplace_back(i, edges[i].top);
    point_to_edges[edges[i].top].emplace_back(i, edges[i].bottom);
  }
  for (auto &point_to_edge : point_to_edges) {
    std::sort(point_to_edge.second.begin(), point_to_edge.second.end(),
              [](const PointEdgeInfo &a, const PointEdgeInfo &b) {
                return b.opposite_point < a.opposite_point;
              });
  }
  return point_to_edges;
}

/**
 * Determines the type of given point based on its adjacent edges.
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
inline PointType get_point_type(point::Point p,
                                const PointToEdges &point_to_edges) {
  if (point_to_edges.count(p) == 0) return PointType::EMPTY;
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
  std::map<point::Segment, std::shared_ptr<Interval>> segment_to_line{};
  std::vector<point::Segment> edges{};
  PointToEdges point_to_edges{};
  std::vector<MonotonePolygon> monotone_polygons{};

  MonotonePolygonBuilder() = default;
  explicit MonotonePolygonBuilder(const std::vector<point::Segment> &edges)
      : edges(edges) {
    this->point_to_edges = get_points_edges(edges);
  }

  /**
   * Retrieves the left and right segments (edges) of a monotone polygon
   * for a given point when edges are sharing top point.
   *
   * @param p The point for which the left and right edges are to be retrieved.
   * @return A pair containing the left and right segments respectively.
   *
   * The function looks up the segments associated with the given point `p`
   * in the `point_to_edges` map. Based on their mutual orientation at the
   * bottom point of one of the edges, it determines which segment is on the
   * left and which is on the right. This orientation check is conducted
   * using a helper function to determine the relative ordering of the edges.
   *
   * If the orientation is counter-clockwise (value 2), the first segment is
   * on the right while the second segment is on the left. Otherwise, the
   * segments maintain their original order.
   */
  std::pair<const point::Segment &, const point::Segment &>
  get_left_right_edges_top(const point::Point &p) {
    auto point_info = point_to_edges.at(p);
    auto fst_idx = point_info[0].edge_index;
    auto snd_idx = point_info[1].edge_index;
    if (intersection::_orientation(edges[fst_idx].bottom, edges[fst_idx].top,
                                   edges[snd_idx].bottom) == 2) {
      return {edges[snd_idx], edges[fst_idx]};
    }
    return {edges[fst_idx], edges[snd_idx]};
  }

  /**
   * Retrieves the left and right edges of a monotone polygon
   * for a given point at its bottom.
   *
   * @param p The point for which the left and right edges are to be retrieved.
   * @return A pair containing the left and right segments (edges) in order.
   *
   * This function looks up the edges associated with the point `p` in the
   * `point_to_edges` map. Based on the orientation of the segments
   * at their top and bottom points, it determines which edge is on the left
   * and which is on the right. The orientation is calculated using a helper
   * function to determine the relative position of one segment relative to the
   * other.
   *
   * If the orientation indicates a counter-clockwise arrangement, the function
   * swaps the order of the edges to ensure the left and right ordering.
   */
  std::pair<const point::Segment &, const point::Segment &>
  get_left_right_edges_bottom(const point::Point &p) {
    auto point_info = point_to_edges.at(p);
    auto fst_idx = point_info[0].edge_index;
    auto snd_idx = point_info[1].edge_index;
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

  void process_end_point(const point::Point &p, const point::Segment &edge_left,
                         const point::Segment &edge_right,
                         const std::shared_ptr<Interval> &interval) {
    for (auto &polygon : interval->polygons_list) {
      polygon->bottom = p;
      monotone_polygons.push_back(*polygon);
    }

    segment_to_line.erase(edge_left);
    segment_to_line.erase(edge_right);
  }

  void process_merge_point(const point::Point &p) {
    auto [edge_left, edge_right] = get_left_right_edges_bottom(p);

    auto left_interval = segment_to_line.at(edge_left);
    auto right_interval = segment_to_line.at(edge_right);

    if (left_interval != right_interval) {
#ifdef DEBUG
      if (right_interval->right_segment == edge_right) {
        std::ostringstream oss;
        oss << "The right edge of merge point should be the left edge of the "
               "right interval.\nGot interval: "
            << *right_interval << " and edge: " << edge_right;
        throw std::runtime_error(oss.str());
      }
#endif
      segment_to_line.erase(edge_right);
      segment_to_line.erase(edge_left);
      left_interval->right_segment = right_interval->right_segment;
#ifdef DEBUG
      if (segment_to_line.count(right_interval->right_segment) == 0) {
        throw std::runtime_error("Segment not found in the map2");
      }
#endif
      segment_to_line[right_interval->right_segment] = left_interval;
      left_interval->last_seen = p;
      left_interval->polygons_list.back()->right.push_back(p);
      right_interval->polygons_list.front()->left.push_back(p);
      for (auto &polygon : right_interval->polygons_list) {
        left_interval->polygons_list.push_back(std::move(polygon));
      }
    } else {
      // This is the end point
      this->process_end_point(p, edge_left, edge_right, left_interval);
    }
  };

  /**
   * Processes a normal point within a monotone polygon.
   *
   * @param p The point to be processed.
   * @param edge_top The top edge segment associated with the point.
   * @param edge_bottom The bottom edge segment associated with the point.
   *
   * This function updates the monotone polygon interval associated with
   * `edge_top` to include the new point `p`. If the interval contains
   * multiple polygons, it closes all but one of them by setting their
   * bottom to `p` and moving them to the list of completed monotone polygons.
   *
   * Depending on whether `edge_top` represents the right segment of the
   * interval, the function then adds the point `p` to the appropriate (left
   * or right) list of points in the remaining polygon and replaces `edge_top`
   * with `edge_bottom` in the segment-to-line map.
   *
   * The function throws a `std::runtime_error` if `edge_top` is not found in
   * the `segment_to_line` map.
   */
  void _process_normal_point(const point::Point &p,
                             const point::Segment &edge_top,
                             const point::Segment &edge_bottom) {
    if (segment_to_line.count(edge_top) == 0) {
      throw std::runtime_error("Segment not found in the map");
    }
    auto interval = segment_to_line.at(edge_top);

    if (interval->polygons_list.size() > 1) {
      // end all the polygons, except 1
      if (edge_top == interval->right_segment) {
        for (auto i = 1; i < interval->polygons_list.size(); i++) {
          interval->polygons_list[i]->bottom = p;
          monotone_polygons.push_back(*interval->polygons_list[i]);
        }

      } else {
        for (auto i = 0; i < interval->polygons_list.size() - 1; i++) {
          interval->polygons_list[i]->bottom = p;
          monotone_polygons.push_back(*interval->polygons_list[i]);
        }
        interval->polygons_list[0] = std::move(interval->polygons_list.back());
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
#ifdef DEBUG
    if (segment_to_line.count(interval->left_segment) == 0) {
      throw std::runtime_error("Left segment not found in the map");
    }
    if (segment_to_line.count(interval->right_segment) == 0) {
      throw std::runtime_error("Right segment not found in the map");
    }
#endif
  };

  /**
   * Processes a normal point within a monotone polygon.
   *
   * @param p The point to be processed.
   *
   * This function identifies the two segments (top and bottom edges)
   * associated with the point `p` using its index in the `point_to_edges` map.
   * It then calls the helper method `_process_normal_point`
   * with the identified segments to handle the actual processing logic.
   */
  void process_normal_point(const point::Point &p) {
    const point::Segment &edge_top =
        edges[point_to_edges.at(p).at(0).edge_index];
    const point::Segment &edge_bottom =
        edges[point_to_edges.at(p).at(1).edge_index];

    _process_normal_point(p, edge_top, edge_bottom);
  }

  void process_start_point(const point::Point &p,
                           const point::Segment &edge_left,
                           const point::Segment &edge_right) {
    auto new_polygon = std::make_unique<MonotonePolygon>(p);
    auto interval =
        std::make_shared<Interval>(p, edge_left, edge_right, new_polygon);
    segment_to_line[edge_left] = interval;
    segment_to_line[edge_right] = interval;
  }

  void process_split_point(const point::Point &p) {
    auto [edge_left, edge_right] = get_left_right_edges_top(p);

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
      auto interval = segment_interval.second;
      if (interval->left_segment.point_on_line_x(p.y) < p.x &&
          interval->right_segment.point_on_line_x(p.y) > p.x) {
        // the point is inside the interval

        // update the sweep line
        auto right_segment = interval->right_segment;
        interval->right_segment = edge_left;
        interval->last_seen = p;
        segment_to_line[edge_left] = interval;

        auto new_interval =
            std::make_shared<Interval>(p, edge_right, right_segment);
        segment_to_line[edge_right] = new_interval;
        segment_to_line[right_segment] = new_interval;

        if (interval->polygons_list.size() == 1) {
          std::unique_ptr<MonotonePolygon> new_polygon;
          if (interval->polygons_list[0]->right.empty()) {
            new_polygon = std::make_unique<MonotonePolygon>(
                interval->polygons_list[0]->top);
          } else {
            new_polygon = std::make_unique<MonotonePolygon>(
                interval->polygons_list[0]->right.back());
          }
          new_polygon->left.push_back(p);
          new_interval->polygons_list.emplace_back(std::move(new_polygon));
          interval->polygons_list[0]->right.push_back(p);
        }

        if (interval->polygons_list.size() >= 2) {
          interval->polygons_list[0]->right.push_back(p);
          interval->polygons_list[interval->polygons_list.size() - 1]
              ->left.push_back(p);
          for (auto i = 1; i < interval->polygons_list.size() - 1; i++) {
            interval->polygons_list[i]->bottom = p;
            monotone_polygons.push_back(*interval->polygons_list[i]);
          }
          new_interval->polygons_list.push_back(
              std::move(interval->polygons_list.back()));
          interval->polygons_list.erase(interval->polygons_list.begin() + 1,
                                        interval->polygons_list.end());
        }
        return;
      }
    }
    // the point is not inside any interval
    this->process_start_point(p, edge_left, edge_right);
  };

  void process_intersection_point(const point::Point &p) {
    auto point_info = point_to_edges.at(p);
    std::set<point::Segment> processed_segments;
    std::vector<point::Segment> segments_to_normal_process;
    std::vector<point::Segment>
        top_segments;  // segments above the intersection point
    std::vector<point::Segment>
        bottom_segments;  // segments below the intersection point

    for (auto &edge_info : point_info) {
      auto &edge = edges[edge_info.edge_index];
      if (processed_segments.count(edge) > 0) {
        continue;
      }
      if (segment_to_line.count(edge) > 0) {
        auto interval = segment_to_line.at(edge);
        auto opposite_edge = interval->opposite_segment(edge);
        if (edge.bottom == p && opposite_edge.bottom == p) {
          process_end_point(p, edge, opposite_edge, interval);
          processed_segments.insert(edge);
          processed_segments.insert(opposite_edge);
          continue;
        }
      }
      if (edge.top == p) {
        bottom_segments.push_back(edge);
      } else {
        top_segments.push_back(edge);
      }
      segments_to_normal_process.push_back(edge);
    }

    // after this loop we should have at most 4 segments

    if (bottom_segments.empty() && top_segments.empty()) {
      // all segments were processed
      return;
    }

    std::sort(top_segments.begin(), top_segments.end(),
              left_right_share_bottom);
    std::sort(bottom_segments.begin(), bottom_segments.end(),
              left_right_share_top);

    auto bottom_begin = bottom_segments.begin();
    auto bottom_end = bottom_segments.end();
    auto top_begin = top_segments.begin();
    if (!top_segments.empty()) {
      if (top_segments.front() ==
          segment_to_line.at(top_segments.front())->right_segment) {
        ++bottom_begin;
        ++top_begin;
        _process_normal_point(p, top_segments.front(), bottom_segments.front());
      }
      if (top_begin != top_segments.end() &&
          top_segments.back() ==
              segment_to_line.at(top_segments.back())->left_segment) {
        --bottom_end;
        _process_normal_point(p, top_segments.back(), bottom_segments.back());
      }
    }

    for (auto it = bottom_begin; it < bottom_end; it += 2) {
      process_start_point(p, *it, *(it + 1));
    }
  };
};

/*
 * Check if the polygon, that all angles have the same orientation
 * do not have self-intersections
 *
 * @param begin Iterator to the first point of the polygon
 * @param end Iterator to the end of the polygon
 * @param centroid Centroid of the polygon
 *
 * @return True if the polygon is simple, false otherwise
 */
template <typename Iterator>
bool is_simple_polygon(Iterator begin, Iterator end, point::Point centroid) {
  double start_angle = std::atan2(begin->y - centroid.y, begin->x - centroid.x);
  double prev_angle = 0;
  begin++;
  for (; begin != end; begin++) {
    double angle =
        std::atan2(begin->y - centroid.y, begin->x - centroid.x) - start_angle;
    if (angle < 0) {
      angle += 2 * M_PI;
    }
    if (angle < prev_angle) {
      return false;
    } else {
      prev_angle = angle;
    }
  }
  return true;
}

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
inline bool _is_convex(const std::vector<point::Point> &polygon) {
  if (polygon.size() < 3) return false;
  if (polygon.size() == 3) return true;
  intersection::Orientation orientation = intersection::Orientation::COLLINEAR;
  intersection::Orientation triangle_orientation;
  size_t idx = 0;
  for (; idx < polygon.size() - 2; idx++) {
    triangle_orientation = intersection::_orientation(
        polygon[idx], polygon[(idx + 1)], polygon[(idx + 2)]);
    if (triangle_orientation != intersection::Orientation::COLLINEAR) {
      orientation = triangle_orientation;
      break;
    }
  }
  if (orientation == intersection::Orientation::COLLINEAR) {
    return false;
  }
  for (; idx < polygon.size() - 2; idx++) {
    triangle_orientation = intersection::_orientation(
        polygon[idx], polygon[(idx + 1)], polygon[(idx + 2)]);
    if (triangle_orientation != 0 && triangle_orientation != orientation) {
      return false;
    }
  }
  triangle_orientation = intersection::_orientation(
      polygon[polygon.size() - 2], polygon[polygon.size() - 1], polygon[0]);
  if (triangle_orientation != 0 && triangle_orientation != orientation) {
    return false;
  }
  triangle_orientation = intersection::_orientation(polygon[polygon.size() - 1],
                                                    polygon[0], polygon[1]);
  if (triangle_orientation != 0 && triangle_orientation != orientation) {
    return false;
  }

  point::Point centroid = point::centroid(polygon);

  if (orientation == intersection::Orientation::COUNTERCLOCKWISE) {
    return is_simple_polygon(polygon.begin(), polygon.end(), centroid);
  } else {
    return is_simple_polygon(polygon.rbegin(), polygon.rend(), centroid);
  }
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
inline std::vector<Triangle> _triangle_convex_polygon(
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
inline void _build_triangles_opposite_edge(std::vector<point::Point> &stack,
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
 * @param stack A reference to a std::vector of point::Point representing
 * the current state of the stack.
 * @param result A reference to a std::vector of PointTriangle where the
 * generated triangles will be stored.
 * @param current_point The incoming point::Point to be considered for new
 * triangles and updating the stack.
 * @param expected_orientation The expected orientation (clockwise or
 * counterclockwise) that will guide the triangle formation.
 */
inline void _build_triangles_current_edge(std::vector<point::Point> &stack,
                                          std::vector<PointTriangle> &result,
                                          point::Point current_point,
                                          int expected_orientation) {
  auto it1 = stack.rbegin();
  auto it2 = stack.rbegin() + 1;
  while (it2 != stack.rend() &&
         intersection::_orientation(*it2, *it1, current_point) ==
             expected_orientation) {
    result.emplace_back(current_point, *it1, *it2);
    ++it1;
    ++it2;
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
inline std::vector<PointTriangle> triangulate_monotone_polygon(
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
inline std::vector<point::Segment> calc_edges(
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
inline std::vector<point::Segment> calc_edges(
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
    if (polygon.back() != polygon.front())
      edges.emplace_back(polygon.back(), polygon.front());
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
inline std::vector<point::Segment> calc_dedup_edges(
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
    if (polygon.back() == polygon.front()) continue;
    edge = point::Segment(polygon.back(), polygon.front());
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
 * @param polygon_list The list of polygons defined by a vector of ov vector of
 * Point objects.
 * @return A new vector of Point objects representing the polygon with added
 * intersection points.
 */
inline std::vector<std::vector<point::Point>> find_intersection_points(
    const std::vector<std::vector<point::Point>> &polygon_list) {
  /* find all edge intersections and add mid-points for all such intersection
   * places*/
  auto edges = calc_edges(polygon_list);

  auto intersections = intersection::_find_intersections(edges);
  if (intersections.empty()) return polygon_list;
  std::unordered_map<std::size_t, std::vector<std::pair<double, point::Point>>>
      intersections_points;
  for (const auto &intersection : intersections) {
    auto inter_points = intersection::_find_intersection(
        edges[intersection.first], edges[intersection.second]);
    for (auto inter_point : inter_points) {
      intersections_points[intersection.first].emplace_back(
          edges[intersection.first].point_projection_factor(inter_point),
          inter_point);
      intersections_points[intersection.second].emplace_back(
          edges[intersection.second].point_projection_factor(inter_point),
          inter_point);
    }
  }
  for (auto &intersections_point : intersections_points) {
    //    points_count += intersections_point.second.size() - 1;
    auto edge = edges[intersections_point.first];
    intersections_point.second.emplace_back(-1, edge.top);
    intersections_point.second.emplace_back(2, edge.bottom);
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
      if (new_polygon.back() != point) new_polygon.push_back(point);
      if (intersections_points.count(i)) {
        auto &new_points = intersections_points[i];
        if (new_points[0].second == point) {
          for (auto it = new_points.begin() + 1; it != new_points.end(); ++it) {
            if (new_polygon.back() != it->second)
              new_polygon.push_back(it->second);
          }
        } else {
          for (auto it = new_points.rbegin() + 1; it != new_points.rend();
               ++it) {
            if (new_polygon.back() != it->second)
              new_polygon.push_back(it->second);
          }
        }
      }
    }
    if (new_polygon.size() > 1 && new_polygon.front() == new_polygon.back())
      new_polygon.pop_back();
    new_polygons_list.push_back(new_polygon);
  }
  return new_polygons_list;
}

inline std::vector<point::Point> find_intersection_points(
    const std::vector<point::Point> &polygon) {
  auto new_polygon = find_intersection_points(
      std::vector<std::vector<point::Point>>({polygon}));
  return new_polygon[0];
}

inline std::vector<point::Point> _sorted_polygons_points(
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
    This is an implementation of the polygon with sweeping lines.
    It assumes that there are no edge intersections but may be a point with
    more than 2 edges.
    Described on this lecture:
    https://www.youtube.com/playlist?list=PLtTatrCwXHzEqzJMaTUFgqoCNllgwk4DH
    and in the book:
    de Berg,M. et al. (2008) Computational Geometry Springer Berlin Heidelberg.
    */
inline std::pair<std::vector<Triangle>, std::vector<point::Point>>
sweeping_line_triangulation(
    const std::vector<std::vector<point::Point>> &polygon_list) {
  std::vector<Triangle> result;
  auto edges = calc_dedup_edges(polygon_list);
  MonotonePolygonBuilder builder(edges);

  PointToEdges point_to_edges = get_points_edges(edges);

  std::vector<point::Point> sorted_points =
      _sorted_polygons_points(polygon_list);

  result.reserve(sorted_points.size() - 2);
  int idx = 0;
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
        // merge two sweeping lines into one
        // merge two intervals into one
        // It may be the end of the polygon, then finish the
        // monotone polygon
        builder.process_merge_point(sorted_point);
        break;
      case PointType::INTERSECTION:
        // this is a merge and split point at the same time
        // this is not described in the original algorithm,
        // but we need it to handle self-intersecting polygons
        builder.process_intersection_point(sorted_point);
        break;
      case PointType::EMPTY:
        // this is a point without edges (removed by deduplication)
        break;
    }
    ++idx;
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

// calculate the triangulation out of a symmetric difference out of a list of
// polygons

inline std::pair<std::vector<Triangle>, std::vector<point::Point>>
triangulate_polygon_face(
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
  // described in this lecture:
  // https://www.youtube.com/playlist?list=PLtTatrCwXHzEqzJMaTUFgqoCNllgwk4DH
  //
  return sweeping_line_triangulation(find_intersection_points(polygon_list));
}

inline std::pair<std::vector<Triangle>, std::vector<point::Point>>
triangulate_polygon_face(const std::vector<point::Point> &polygon) {
  // #if DDEBUG
  //     try{
  // #endif
  return triangulate_polygon_face(
      std::vector<std::vector<point::Point>>({polygon}));
  // #if DDEBUG
  //       } catch (const std::exception &e) {
  //       std::cerr << "Polygon: [";
  //       for (const auto &point : polygon) {
  //         std::cerr << "(" << point.x << ", " << point.y  << "), ";
  //       }
  //       std::cerr << "]" << std::endl;
  //       std::cerr << "Error: " << e.what() << std::endl;
  //       throw e;
  //     }
  // #endif
}

inline point::Point::coordinate_t vector_length(point::Point p1,
                                                point::Point p2) {
  return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) +
                   (p1.y - p2.y) * (p1.y - p2.y));
}

struct PathTriangulation {
  std::vector<Triangle> triangles;
  std::vector<point::Point> centers;
  std::vector<point::Vector> offsets;

  void reserve(std::size_t size) {
    triangles.reserve(size);
    centers.reserve(size);
    offsets.reserve(size);
  }

  void fix_triangle_orientation() {
    point::Point p1{}, p2{}, p3{};
    for (auto &triangle : triangles) {
      p1 = centers[triangle.x] + offsets[triangle.x];
      p2 = centers[triangle.y] + offsets[triangle.y];
      p3 = centers[triangle.z] + offsets[triangle.z];
      if ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x) < 0) {
        std::swap(triangle.x, triangle.z);
      }
      //      if (intersection::_orientation(p1, p2, p3) == 1) {
      //        std::swap(triangle.x, triangle.z);
      //      }
    }
  }
};

inline std::pair<point::Point::coordinate_t, point::Point::coordinate_t>
sign_abs(point::Point::coordinate_t x) {
  if (x < 0) {
    return {-1, -x};
  }
  if (x > 0) {
    return {1, x};
  }
  return {0, 0};
}

inline point::Point::coordinate_t add_triangles_for_join(
    PathTriangulation &triangles, point::Point p1, point::Point p2,
    point::Point p3, point::Point::coordinate_t prev_length, double cos_limit,
    bool bevel) {
  std::size_t idx = triangles.offsets.size();
  point::Vector mitter{};
  point::Point::coordinate_t length = vector_length(p2, p3);
  point::Vector p1_p2_diff_norm = (p2 - p1) / prev_length;
  point::Vector p2_p3_diff_norm = (p3 - p2) / length;

  point::Point::coordinate_t cos_angle = p1_p2_diff_norm.x * p2_p3_diff_norm.x +
                                         p1_p2_diff_norm.y * p2_p3_diff_norm.y;
  point::Point::coordinate_t sin_angle = p1_p2_diff_norm.x * p2_p3_diff_norm.y -
                                         p1_p2_diff_norm.y * p2_p3_diff_norm.x;

  triangles.centers.push_back(p2);
  triangles.centers.push_back(p2);

  if (sin_angle == 0) {
    mitter = {p1_p2_diff_norm.y / 2, -p1_p2_diff_norm.x / 2};
  } else {
    point::Point::coordinate_t scale_factor = 1 / sin_angle;
    if (bevel || cos_angle < cos_limit) {
      /* Bevel join
       * There is a need to check if inner vector is not to long
       * See https://github.com/napari/napari/pull/7268#user-content-bevel-cut
       */
      auto [sign, mag] = sign_abs(scale_factor);
      scale_factor =
          sign * (float)0.5 * std::min(mag, std::min(prev_length, length));
    }
    mitter = (p1_p2_diff_norm - p2_p3_diff_norm) * scale_factor * 0.5;
  }
  if (bevel || cos_angle < cos_limit) {
    triangles.centers.push_back(p2);
    triangles.triangles.emplace_back(idx, idx + 1, idx + 2);
    if (sin_angle < 0) {
      triangles.offsets.push_back(mitter);
      triangles.offsets.emplace_back(-p1_p2_diff_norm.y * 0.5,
                                     p1_p2_diff_norm.x * 0.5);
      triangles.offsets.emplace_back(-p2_p3_diff_norm.y * 0.5,
                                     p2_p3_diff_norm.x * 0.5);
      triangles.triangles.emplace_back(idx, idx + 2, idx + 3);
      triangles.triangles.emplace_back(idx + 2, idx + 3, idx + 4);
    } else {
      triangles.offsets.emplace_back(p1_p2_diff_norm.y * 0.5,
                                     -p1_p2_diff_norm.x * 0.5);
      triangles.offsets.push_back(-mitter);
      triangles.offsets.emplace_back(p2_p3_diff_norm.y * 0.5,
                                     -p2_p3_diff_norm.x * 0.5);
      triangles.triangles.emplace_back(idx + 1, idx + 2, idx + 3);
      triangles.triangles.emplace_back(idx + 1, idx + 3, idx + 4);
    }
  } else {
    triangles.offsets.push_back(mitter);
    triangles.offsets.push_back(-mitter);
    triangles.triangles.emplace_back(idx, idx + 1, idx + 2);
    triangles.triangles.emplace_back(idx + 1, idx + 2, idx + 3);
  }

  return length;
}

inline PathTriangulation triangulate_path_edge(
    const std::vector<point::Point> &path, bool closed, float limit,
    bool bevel) {
  if (path.size() < 2)
    return {{{0, 1, 3}, {1, 3, 2}},
            {path[0], path[0], path[0], path[0]},
            {{0, 0}, {0, 0}, {0, 0}, {0, 0}}};
  PathTriangulation result;
  point::Vector norm_diff{};
  result.reserve(path.size() * 3);
  double cos_limit = 1.0 / (limit * limit / 2) - 1.0;
  point::Point::coordinate_t prev_length{};

  if (closed) {
    prev_length = vector_length(path[0], path[path.size() - 1]);
    prev_length =
        add_triangles_for_join(result, path[path.size() - 1], path[0], path[1],
                               prev_length, cos_limit, bevel);
  } else {
    prev_length = vector_length(path[0], path[1]);
    norm_diff = (path[1] - path[0]) / prev_length;
    result.centers.push_back(path[0]);
    result.centers.push_back(path[0]);
    result.offsets.emplace_back(norm_diff.y * 0.5, -norm_diff.x * 0.5);
    result.offsets.push_back(-result.offsets.back());
    result.triangles.emplace_back(0, 1, 2);
    result.triangles.emplace_back(1, 2, 3);
  }
  for (std::size_t i = 1; i < path.size() - 1; i++) {
    prev_length =
        add_triangles_for_join(result, path[i - 1], path[i], path[i + 1],
                               prev_length, cos_limit, bevel);
  }
  if (closed) {
    add_triangles_for_join(result, path[path.size() - 2], path[path.size() - 1],
                           path[0], prev_length, cos_limit, bevel);
    result.centers.push_back(result.centers[0]);
    result.centers.push_back(result.centers[0]);
    result.offsets.push_back(result.offsets[0]);
    result.offsets.push_back(result.offsets[1]);
  } else {
    norm_diff = (path[path.size() - 1] - path[path.size() - 2]) / prev_length;
    result.centers.push_back(path[path.size() - 1]);
    result.centers.push_back(path[path.size() - 1]);
    result.offsets.emplace_back(norm_diff.y * 0.5, -norm_diff.x * 0.5);
    result.offsets.push_back(-result.offsets.back());
  }
  result.fix_triangle_orientation();
  return result;
}

/**
 * Represents an edge in a graph structure used for polygon processing.
 * Each edge contains a reference to its opposite point and a flag to track
 * if it has been visited during graph traversal.
 */
struct GraphEdge {
  point::Point opposite_point;
  bool visited;
  explicit GraphEdge(point::Point p) : opposite_point(p), visited(false) {}
};

/**
 * Represents a node in a graph structure used for polygon processing.
 * Each node contains its edges, a sub-index for traversal tracking,
 * and a visited flag for graph traversal.
 */
struct GraphNode {
  std::vector<GraphEdge> edges;
  std::size_t sub_index;
  bool visited;

  GraphNode() : sub_index(0), visited(false) {}
};

/**
 * Splits a polygon into sub-polygons by identifying and removing edges that
 * appear more than once in the polygon's edge list.
 *
 * This function processes the given polygon and separates it wherever an
 * edge is repeated. It generates a collection of sub-polygons such that each
 * resulting sub-polygon contains unique edges. This operation can help to
 * resolve ambiguities in complex or self-intersecting polygons.
 *
 * Note: Polygons with exactly 2 points are treated as single line segments and
 * are not processed for edge deduplication between their first and last points.
 *
 * @param polygon The input polygon represented as a list of edges.
 *
 * @return A vector of sub-polygons, where each sub-polygon is free of repeated
 * edges.
 */
inline std::vector<std::vector<point::Point>> split_polygon_on_repeated_edges(
    const std::vector<point::Point> &polygon) {
  if (polygon.size() < 3) return {polygon};
  auto edges_dedup = calc_dedup_edges({polygon});
  std::vector<std::vector<point::Point>> result;
  point::Segment segment;

  std::unordered_set edges_set(edges_dedup.begin(), edges_dedup.end());
  std::unordered_map<point::Point, GraphNode> edges_map;
  for (std::size_t i = 0; i < polygon.size() - 1; i++) {
    segment = {polygon[i], polygon[(i + 1)]};
    if (edges_set.count(segment) > 0) {
      edges_map[polygon[i]].edges.emplace_back(polygon[i + 1]);
    }
  }
  segment = {polygon.back(), polygon.front()};
  if (edges_set.count(segment) > 0) {
    edges_map[polygon.back()].edges.emplace_back(polygon.front());
  }
  for (auto &edge : edges_map) {
    if (edge.second.visited) {
      continue;
    }
    edge.second.visited = true;
    std::vector<point::Point> new_polygon;
    new_polygon.push_back(edge.first);
    auto *current_edge = &edge.second;
    while (current_edge->sub_index < current_edge->edges.size()) {
      auto *prev = current_edge;
      auto next_point =
          current_edge->edges[current_edge->sub_index].opposite_point;
      current_edge = &edges_map.at(next_point);
      prev->sub_index++;
      current_edge->visited = true;
      new_polygon.push_back(next_point);
    }
    while (new_polygon.front() == new_polygon.back()) {
      new_polygon.pop_back();
    }
    result.push_back(new_polygon);
  }
  return result;
}

}  // namespace partsegcore::triangulation

#endif  // PARTSEGCORE_TRIANGULATE_H
