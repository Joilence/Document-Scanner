#include <vector>

/**
 * @brief 2D Point
 * @details Stores coordinates and its votes
 */
struct Point {
  double x, y, value;
  Point(double _x, double _y, double _value) : x(_x), y(_y), value(_value) {}
  Point(double _x, double _y) : x(_x), y(_y) { value = 0; }
  Point() {};
};
// bool comp(Point* a, Point* b) { return (*a).value > (*b).value; }
bool comp_ds(Point a, Point b) { return a.value > b.value; }
bool comp_as(Point a, Point b) { return a.value < b.value; }

/**
 * @brief 2D Line
 * @details Stores k and b, and two points of it to reduce time complexity
 */
struct Line {
  double k, b;
  double x0, y0, x1, y1;
  std::vector<Point> intersections;
  std::vector<Point> borderPoints;
  Line(double _k, double _b, int _x0, int _y0, int _x1, int _y1)
      : k(_k), b(_b), x0(_x0), y0(_y0), x1(_x1), y1(_y1) {}
};