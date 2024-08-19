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
