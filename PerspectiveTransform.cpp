// #include "OptImgProc.cpp"
#include </usr/local/include/eigen3/Eigen/Dense>
#include "datatype.cpp"
#include <vector>
#include <iostream>

using namespace Eigen;
using namespace std;

MatrixXd getPerspectiveTransform(vector<Point> src, vector<Point> dst) {
  MatrixXd A(8, 8);
  MatrixXd B(8, 1);

  double a[8][8], b[8];

  for (int i = 0; i < 4; ++i) {
    a[i][0] = a[i + 4][3] = src[i].x;
    a[i][1] = a[i + 4][4] = src[i].y;
    a[i][2] = a[i + 4][5] = 1;
    a[i][3] = a[i][4] = a[i][5] = a[i + 4][0] = a[i + 4][1] = a[i + 4][2] = 0;
    a[i][6] = -src[i].x * dst[i].x;
    a[i][7] = -src[i].y * dst[i].x;
    a[i + 4][6] = -src[i].x * dst[i].y;
    a[i + 4][7] = -src[i].y * dst[i].y;
    b[i] = dst[i].x;
    b[i + 4] = dst[i].y;
  }

  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      A(i, j) = a[i][j];
    }
  }
  for (int i = 0; i < 8; ++i) {
    B(i, 0) = b[i];
  }

  cout << A << endl;
  cout << B << endl;

  MatrixXd X = A.bdcSvd(ComputeThinU | ComputeThinV).solve(B);
  MatrixXd M(9, 1);
  for (int i = 0; i < 8; ++i) {
    M(i, 0) = X(i, 0);
  }
  M(8, 0) = 1;
  return M;
}

