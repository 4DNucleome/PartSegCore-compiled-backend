#include <set>
#include <unordered_map>
#include <unordered_set>


struct Point{
    float x;
    float y;
    bool operator==(const Point& p) const {
        return x == p.x && y == p.y;
    }
    Point(float x, float y) : x(x), y(y) {}
    Point() {}

    bool operator<(const Point& p) const {
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

    Event(float x, float y, int index, bool is_left) : x(x), y(y), index(index), is_left(is_left) {}
    Event(const Point& p, int index, bool is_left): x(p.x), y(p.y), index(index), is_left(is_left) {}
    Event() {}

    bool operator<(const Event& e) const {
        if (x == e.x) {
            if(y==e.y){
                if (is_left == e.is_left){
                    return index < e.index;
                }
                return is_left > e.is_left;
            }
            return y < e.y;
        }
        return x < e.x;
    }
};

struct PairHash{
    std::size_t operator()(const std::pair<int, int>& p) const {
        return std::hash<int>()(p.first) * 31 + std::hash<int>()(p.second);
    }
};


struct Segment{
    Point left;
    Point right;
    Segment(Point p1, Point p2){
        if (p1 < p2){
            this->left = p1;
            this->right = p2;
        } else {
            this->left = p2;
            this->right = p1;
        }

    }
    Segment() {}
};


struct Triangle{
    int x;
    int y;
    int z;
    Triangle(int x, int y, int z) : x(x), y(y), z(z) {}
    Triangle() {}
};


bool point_eq(const Point& p, const Point& q){
    return p.x == q.x && p.y == q.y;
}

bool cmp_point(const Point& p, const Point& q){
    if(p.x == q.x){
        return p.y < q.y;
    }
    return p.x < q.x;
}

bool cmp_event(const Event& p, const Event& q){
    return p < q;
}


bool _on_segment(const Point& p, const Point& q, const Point& r){
    if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
        q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y))
        return true;
    return false;
}

int _orientation(const Point& p, const Point& q, const Point& r){
    float val = (q.y - p.y) * (r.x - q.x) -
                (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0;
    return (val > 0) ? 1 : 2;
}


bool _do_intersect(const Segment& s1, const Segment& s2){
    const Point& p1 = s1.left;
    const Point& q1 = s1.right;
    const Point& p2 = s2.left;
    const Point& q2 = s2.right;
    if (point_eq(p1, p2) || point_eq(p1, q2) || point_eq(q1, p2) || point_eq(q1, q2))
        return false;

    int o1 = _orientation(p1, q1, p2);
    int o2 = _orientation(p1, q1, q2);
    int o3 = _orientation(p2, q2, p1);
    int o4 = _orientation(p2, q2, q1);

    if (o1 != o2 && o3 != o4)
        return true;

    if (o1 == 0 && _on_segment(p1, p2, q1)) return true;
    if (o2 == 0 && _on_segment(p1, q2, q1)) return true;
    if (o3 == 0 && _on_segment(p2, p1, q2)) return true;
    if (o4 == 0 && _on_segment(p2, q1, q2)) return true;

    return false;
}

std::set<Event>::iterator pred(std::set<Event>& s, std::set<Event>::iterator it){
    if(it == s.begin()){
        return s.end();
    }
    return --it;
}

std::set<Event>::iterator succ(std::set<Event>& s, std::set<Event>::iterator it){
    return ++it;
}

std::unordered_set<std::pair<int, int>, PairHash>
_find_intersections(const std::vector<Segment>& segments){
    std::unordered_set<std::pair<int, int>, PairHash> intersections;
    std::vector<Event> events;
    std::set<Event> active;
    events.reserve(2 * segments.size());
    for (int i=0; i<segments.size(); i++){
        events.push_back(Event(segments[i].left, i, true));
        events.push_back(Event(segments[i].right, i, false));
    }
    std::sort(events.begin(), events.end(), cmp_event);

    for (auto event = events.begin(); event != events.end(); event++){
        if (event->is_left){
            auto next = active.lower_bound(*event);
            auto prev = pred(active, next);
            if (next != active.end() && _do_intersect(segments[event->index], segments[next->index])){
                if (event->index < next->index){
                    intersections.emplace(event->index, next->index);
                } else {
                    intersections.emplace(next->index, event->index);
                }
            }
            if (prev != active.end() && _do_intersect(segments[event->index], segments[prev->index])){
                if (event->index < prev->index){
                    intersections.emplace(event->index, prev->index);
                } else {
                    intersections.emplace(prev->index, event->index);
                }
            }
            active.insert(*event);
        } else {
            auto it = active.find(Event(segments[event->index].left, event->index, true));
            auto next = succ(active, it);
            auto prev = pred(active, it);
            if (next != active.end() && prev != active.end() && _do_intersect(segments[next->index], segments[prev->index])){
                if (next->index < prev->index){
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

Point _find_intersection(const Segment& s1, const Segment& s2){
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

bool _is_convex(const std::vector<Point>& polygon){
    int orientation = 0;
    int triangle_orientation;
    for (int i=0; i < polygon.size()-2; i++){
        triangle_orientation = _orientation(polygon[i], polygon[i+1], polygon[i+2]);
        if (triangle_orientation == 0)
            continue;
        if (orientation == 0)
            orientation = triangle_orientation;
        else if (orientation != triangle_orientation)
            return false;
    }
    triangle_orientation = _orientation(
        polygon[polygon.size() -2],
        polygon[polygon.size() -1],
        polygon[0]
    );
    if (triangle_orientation != 0 && triangle_orientation != orientation)
        return false;
    triangle_orientation = _orientation(
        polygon[polygon.size() - 1],
        polygon[0],
        polygon[1]
    );
    if (triangle_orientation != 0 && triangle_orientation != orientation)
        return false;
    return true;
}