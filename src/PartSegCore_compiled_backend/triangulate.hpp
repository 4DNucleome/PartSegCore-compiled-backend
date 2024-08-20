#include <set>
#include <unordered_map>
#include <unordered_set>

struct Event {
    float x;
    float y;
    int index;
    bool is_left;

    Event(float x, float y, int index, bool is_left) : x(x), y(y), index(index), is_left(is_left) {}
    Event() {}

    bool operator<(const Event& e) const {
        if(y==e.y){
            if(x==e.x){
                if (is_left == e.is_left){
                    return index < e.index;
                }
                return is_left < e.is_left;
            }
            return x < e.x;
        }
        return y < e.y;
    }
};

struct PairHash{
    std::size_t operator()(const std::pair<int, int>& p) const {
        return std::hash<int>()(p.first) * 31 + std::hash<int>()(p.second);
    }
};

struct Point{
    float x;
    float y;
    bool operator==(const Point& p) const {
        return x == p.x && y == p.y;
    }
    Point(float x, float y) : x(x), y(y) {}
    Point() {}

};

struct Segment{
    Point left;
    Point right;
    Segment(Point left, Point right) : left(left), right(right) {}
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

std::unordered_map<std::pair<int, int>, int, PairHash> _find_intersections(const std::vector<Segment>& segments){
    std::unordered_map<std::pair<int, int>, int, PairHash> intersections;
    std::vector<Event> events;
    events.reserve(2 * segments.size());
}
