static const signed char neighbourhood[26][3] = {
    {0, -1, 0},  {0, 0, -1},  {0, 1, 0},   {0, 0, 1},
    {-1, 0, 0},  {1, 0, 0},

    {-1, -1, 0}, {1, -1, 0},  {-1, 1, 0},  {1, 1, 0},
    {-1, 0, -1}, {1, 0, -1},  {-1, 0, 1},  {1, 0, 1},
    {0, -1, -1}, {0, 1, -1},  {0, -1, 1},  {0, 1, 1},

    {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1},
    {1, 1, -1},  {1, -1, 1},  {-1, 1, 1},  {1, 1, 1}};

static const signed char neighbourhood2d[8][3] = {
    {0, -1, 0},  {0, 0, -1}, {0, 1, 0},  {0, 0, 1},
    {0, -1, -1}, {0, 1, -1}, {0, -1, 1}, {0, 1, 1},
};

static const char neigh_level[] = {4, 6, 18, 26, 26};

static const float distance[26] = {
    1,     1,     1,     1,     1,     1,

    1.414, 1.414, 1.414, 1.414, 1.414, 1.414, 1.414, 1.414,
    1.414, 1.414, 1.414, 1.414,

    1.732, 1.732, 1.732, 1.732, 1.732, 1.732, 1.732, 1.732,
};

static const signed char neighbourhood_connection[13][3] = {
    {1, 0, 0},  {0, 1, 0},  {0, 0, 1},  {1, 1, 0},  {1, 0, 1},
    {0, 1, 1},  {0, 1, -1}, {1, 0, -1}, {1, -1, 0}, {1, 1, 1},
    {1, -1, 1}, {1, 1, -1}, {1, -1, -1}};

static const char neigh_connection_level[] = {3, 9, 13};
