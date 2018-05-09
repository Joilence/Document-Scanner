#include <iostream>
#include "PerspectiveTransform.cpp"

using namespace std;

int main() {
  int img_width = 445;
  int img_height = 600;
  vector<Point> corners(4);
  corners[0] = Point(0, 0);
  corners[1] = Point(img_width - 1, 0);
  corners[2] = Point(0, img_height - 1);
  corners[3] = Point(img_width - 1, img_height - 1);
  vector<Point> corners_trans(4);
  corners_trans[0] = Point(92, 77);
  corners_trans[1] = Point(407, 85);
  corners_trans[2] = Point(97, 526);
  corners_trans[3] = Point(396, 511);

  MatrixXd transform = getPerspectiveTransform(corners, corners_trans);
  cout << "transform = \n";
  cout << transform << endl;

  // vector<Point> points, points_trans;
  // for (int i = 0; i < img_width; ++i) {
  //   for (int j = 0; j < img_height; ++i) {
  //     points.push_back(Point(i, j));
  //   }
  // }


  // // Transform Points
  // for (int i = 0; i < points.size(); ++i) {
  //   double x = points[i].x;
  //   double y = points[i].y;
  //   // double denominator = a13 * x + a23 * y + a33;
  //   double denominator = transform(0, 2) * x + transform(1, 2) * y + transform(2, 2);
  //   // points[i].x = (a11 * x + a21 * y + a31) / denominator;
  //   points[i].x = (transform(0, 0) * x + transform(1, 0) * y + transform(2, 0)) / denominator;
  //   // points[i].y = (a12 * x + a22 * y + a32) / denominator;
  //   points[i].y = (transform(0, 1) * x + transform(1, 1) * y + transform(2, 1)) / denominator;
  // }

  // for (int i = 0; i < img_height; i++) {
  //   uchar* t = img_trans.ptr<uchar>(i);
  //   for (int j = 0; j < img_width; j++) {
  //     int tmp = i * img_width + j;
  //     int x = ponits[tmp * 2];
  //     int y = ponits[tmp * 2 + 1];
  //     if (x < 0 || x > (img_width - 1) || y < 0 || y > (img_height - 1))
  //       continue;
  //     uchar* p = img.ptr<uchar>(y);
  //     t[j * 3] = p[x * 3];
  //     t[j * 3 + 1] = p[x * 3 + 1];
  //     t[j * 3 + 2] = p[x * 3 + 2];
  //   }
  // }
  // perspectiveTransform(points, points_trans, transform);

  return 0;
}