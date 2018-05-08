//#include <iostream>
#include "PerspectiveTransform.cpp"

using namespace std;

int main() {
  int img_width = 1500;
  int img_height = 2000;
  vector<Point> corners(4);
  corners[0] = Point(0, 0);
  corners[1] = Point(img_width - 1, 0);
  corners[2] = Point(0, img_height - 1);
  corners[3] = Point(img_width - 1, img_height - 1);
  vector<Point> corners_trans(4);
  corners_trans[0] = Point(150, 250);
  corners_trans[1] = Point(771, 0);
  corners_trans[2] = Point(0, img_height - 1);
  corners_trans[3] = Point(650, img_height - 1);
  MatrixXd transform = getPerspectiveTransform(corners, corners_trans);
  cout << transform;
  return 0;
}