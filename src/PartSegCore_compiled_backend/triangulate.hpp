#include <set>
#include <unordered_map>
#include <unordered_set>

struct Point {
  float x;
  float y;
  bool operator==(const Point &p) const { return x == p.x && y == p.y; }
  Point(float x, float y) : x(x), y(y) {}
  Point() {}

  bool operator<(const Point &p) const {
    if (this->x == p.x) {
      return this->y < p.y;
    }
    return this->x < p.x;
  }
};

struct Event {
  float x;
  float y;
  int index;
  bool is_left;

  Event(float x, float y, int index, bool is_left)
      : x(x), y(y), index(index), is_left(is_left) {}
  Event(const Point &p, int index, bool is_left)
      : x(p.x), y(p.y), index(index), is_left(is_left) {}
  Event() {}

  bool operator<(const Event &e) const {
    if (x == e.x) {
      if (y == e.y) {
        if (is_left == e.is_left) {
          return index < e.index;
        }
        return is_left > e.is_left;
      }
      return y < e.y;
    }
    return x < e.x;
  }
};

struct PairHash {
  std::size_t operator()(const std::pair<int, int>& p) const {
    return std::hash<int>()(p.first) * 31 + std::hash<int>()(p.second);
  }
};


struct PointHash {
  std::size_t operator()(const Point& p) const{
    return std::hash<int>()(p.x) * 31 + std::hash<int>()(p.y);
  }
};


struct Segment {
  Point left;
  Point right;
  Segment(Point p1, Point p2) {
    if (p1 < p2) {
      this->left = p1;
      this->right = p2;
    } else {
      this->left = p2;
      this->right = p1;
    }
  }
  Segment() {}
};

struct Triangle {
  int x;
  int y;
  int z;
  Triangle(int x, int y, int z) : x(x), y(y), z(z) {}
  Triangle() {}
};

typedef std::unordered_map<Point, std::vector<std::pair<int, Point>>, PointHash> PointToEdges;

bool point_eq(const Point &p, const Point &q) {
  return p.x == q.x && p.y == q.y;
}

bool cmp_point(const Point &p, const Point &q) { return p < q; }

bool cmp_pair_point(const std::pair<Point, int> &p,
                    const std::pair<Point, int> &q) {
  return p.first < q.first;
}

bool cmp_event(const Event &p, const Event &q) { return p < q; }

bool _on_segment(const Point &p, const Point &q, const Point &r) {
  if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
      q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y))
    return true;
  return false;
}

int _orientation(const Point &p, const Point &q, const Point &r) {
  float val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
  if (val == 0)
    return 0;
  return (val > 0) ? 1 : 2;
}

bool _do_intersect(const Segment &s1, const Segment &s2) {
  const Point &p1 = s1.left;
  const Point &q1 = s1.right;
  const Point &p2 = s2.left;
  const Point &q2 = s2.right;
  if (point_eq(p1, p2) || point_eq(p1, q2) || point_eq(q1, p2) ||
      point_eq(q1, q2))
    return false;

  int o1 = _orientation(p1, q1, p2);
  int o2 = _orientation(p1, q1, q2);
  int o3 = _orientation(p2, q2, p1);
  int o4 = _orientation(p2, q2, q1);

  if (o1 != o2 && o3 != o4)
    return true;

  if (o1 == 0 && _on_segment(p1, p2, q1))
    return true;
  if (o2 == 0 && _on_segment(p1, q2, q1))
    return true;
  if (o3 == 0 && _on_segment(p2, p1, q2))
    return true;
  if (o4 == 0 && _on_segment(p2, q1, q2))
    return true;

  return false;
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

std::unordered_set<std::pair<int, int>, PairHash>
_find_intersections(const std::vector<Segment> &segments) {
  std::unordered_set<std::pair<int, int>, PairHash> intersections;
  std::vector<Event> events;
  std::set<Event> active;
  events.reserve(2 * segments.size());
  for (int i = 0; i < segments.size(); i++) {
    events.push_back(Event(segments[i].left, i, true));
    events.push_back(Event(segments[i].right, i, false));
  }
  std::sort(events.begin(), events.end(), cmp_event);

  for (auto event = events.begin(); event != events.end(); event++) {
    if (event->is_left) {
      auto next = active.lower_bound(*event);
      auto prev = pred(active, next);
      if (next != active.end() &&
          _do_intersect(segments[event->index], segments[next->index])) {
        if (event->index < next->index) {
          intersections.emplace(event->index, next->index);
        } else {
          intersections.emplace(next->index, event->index);
        }
      }
      if (prev != active.end() &&
          _do_intersect(segments[event->index], segments[prev->index])) {
        if (event->index < prev->index) {
          intersections.emplace(event->index, prev->index);
        } else {
          intersections.emplace(prev->index, event->index);
        }
      }
      active.insert(*event);
    } else {
      auto it =
          active.find(Event(segments[event->index].left, event->index, true));
      auto next = succ(active, it);
      auto prev = pred(active, it);
      if (next != active.end() && prev != active.end() &&
          _do_intersect(segments[next->index], segments[prev->index])) {
        if (next->index < prev->index) {
          intersections.emplace(next->index, prev->index);
        } else {
          intersections.emplace(prev->index, next->index);
        }
      }
      active.erase(it);
    }
  }
  return intersections;
}

Point _find_intersection(const Segment &s1, const Segment &s2) {
  float a1, b1, c1, a2, b2, c2, det, x, y;
  a1 = s1.right.y - s1.left.y;
  b1 = s1.left.x - s1.right.x;
  c1 = a1 * s1.left.x + b1 * s1.left.y;
  a2 = s2.right.y - s2.left.y;
  b2 = s2.left.x - s2.right.x;
  c2 = a2 * s2.left.x + b2 * s2.left.y;
  det = a1 * b2 - a2 * b1;
  if (det == 0)
    return Point(0, 0);
  x = (b2 * c1 - b1 * c2) / det;
  y = (a1 * c2 - a2 * c1) / det;
  return Point(x, y);
}

bool _is_convex(const std::vector<Point> &polygon) {
  int orientation = 0;
  int triangle_orientation;
  for (int i = 0; i < polygon.size() - 2; i++) {
    triangle_orientation =
        _orientation(polygon[i], polygon[i + 1], polygon[i + 2]);
    if (triangle_orientation == 0)
      continue;
    if (orientation == 0)
      orientation = triangle_orientation;
    else if (orientation != triangle_orientation)
      return false;
  }
  triangle_orientation = _orientation(polygon[polygon.size() - 2],
                                      polygon[polygon.size() - 1], polygon[0]);
  if (triangle_orientation != 0 && triangle_orientation != orientation)
    return false;
  triangle_orientation =
      _orientation(polygon[polygon.size() - 1], polygon[0], polygon[1]);
  if (triangle_orientation != 0 && triangle_orientation != orientation)
    return false;
  return true;
}

std::vector<Triangle>
_triangle_convex_polygon(const std::vector<Point> &polygon) {
  std::vector<Triangle> result;
  for (int i = 1; i < polygon.size() - 1; i++) {
    if (_orientation(polygon[0], polygon[i], polygon[i + 1]) != 0) {
      result.push_back(Triangle(0, i, i + 1));
    }
  }
  return result;
}

std::vector<Segment> calc_edges(const std::vector<Point>& polygon){
  std::vector<Segment> edges;
  edges.reserve(polygon.size());
  for (int i = 0; i < polygon.size() - 1; i++) {
    edges.push_back(Segment(polygon[i], polygon[i + 1]));
  }
  edges.push_back(Segment(polygon[polygon.size() - 1], polygon[0]));
  return edges;
}

std::vector<Point> find_intersection_points(const std::vector<Point>& polygon) {
  /* find all edge intersedions and add mid points for all such intersection
   * place*/
  auto edges = calc_edges(polygon);

  auto intersections = _find_intersections(edges);
  if (intersections.size() == 0)
    return polygon;
  std::unordered_map<int, std::vector<Point>> intersections_points;
  for (auto p = intersections.begin(); p != intersections.end(); p++) {
    auto inter_point = _find_intersection(edges[p->first], edges[p->second]);
    intersections_points[p->first].push_back(inter_point);
    intersections_points[p->second].push_back(inter_point);
  }
  int points_count = polygon.size();
  for (auto iter_inter = intersections_points.begin();
       iter_inter != intersections_points.end(); iter_inter++) {
    points_count += iter_inter->second.size() - 1;
    iter_inter->second.push_back(edges[iter_inter->first].right);
    iter_inter->second.push_back(edges[iter_inter->first].left);
    std::sort(iter_inter->second.begin(), iter_inter->second.end(), cmp_point);
  }

  std::vector<Point> new_polygon;
  new_polygon.reserve(points_count);
  for (int i = 0; i < polygon.size(); i++) {
    auto point = polygon[i];
    new_polygon.push_back(point);
    if (intersections_points.count(i)) {
      auto new_points = intersections_points[i];
      if (new_points[0] == point) {
        for (int j = 1; j < new_points.size() - 1; j++) {
          new_polygon.push_back(new_points[j]);
        }
      } else {
        for (int j = new_points.size() - 2; j > 0; j++) {
          new_polygon.push_back(new_points[j]);
        }
      }
    }
  }
  return new_polygon;
}

bool is_normal_point(
  Point p,
  PointToEdges point_to_edges
  ){
  if (point_to_edges.at(p).size() != 2)
    return false;
  auto edges = point_to_edges.at(p);
  return (edges[0].second < p && p < edges[1].second) || (p < edges[0].second && p < edges[0].second);
}

bool is_merge_point(
  Point p,
  PointToEdges point_to_edges
  ){
  if (point_to_edges.at(p).size() != 2)
    return false;
  auto edges = point_to_edges.at(p);
  return (edges[0].second < p && edges[1].second < p);
}

bool is_split_point(
  Point p,
  PointToEdges point_to_edges
  ){
  if (point_to_edges.at(p).size() != 2)
    return false;
  auto edges = point_to_edges.at(p);
  return (p < edges[0].second && p < edges[1].second);
}


std::pair<std::vector<Triangle>, std::vector<Point>>
_triangulate_polygon(const std::vector<Point>& polygon) {
  if (_is_convex(polygon))
    return std::make_pair(_triangle_convex_polygon(polygon), polygon);

  // Implement sweeping line algorithm for triangulation
  // described on this lecture:
  // https://www.youtube.com/playlist?list=PLtTatrCwXHzEqzJMaTUFgqoCNllgwk4DH
  // 
  auto new_polygon = find_intersection_points(polygon);
  auto edges = calc_edges(new_polygon);
  PointToEdges point_to_edges; 
  for (int i=0; i<edges.size(); i++){
    point_to_edges[edges[i].left].push_back(std::make_pair(i, edges[i].right));
    point_to_edges[edges[i].right].push_back(std::make_pair(i, edges[i].left));
  }
  std::vector<Point> sorted_points = new_polygon;
  // copy to avoid modification of original vector
  std::sort(sorted_points.begin(), sorted_points.end(), cmp_point);
  for (auto point = sorted_points.begin(); point != sorted_points.end(); point++){
    is_normal_point(*point, point_to_edges);
  }
}