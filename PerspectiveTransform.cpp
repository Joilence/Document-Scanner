// #include "OptImgProc.cpp"
#include </usr/local/include/eigen3/Eigen/Dense>
#include "datatype.cpp"
#include <vector>
#include <iostream>

using namespace Eigen;
using namespace std;


/* Calculates coefficients of perspective transformation
 * which maps (xi,yi) to (ui,vi), (i=1,2,3,4):
 *
 *      c00*xi + c01*yi + c02
 * ui = ---------------------
 *      c20*xi + c21*yi + c22
 *
 *      c10*xi + c11*yi + c12
 * vi = ---------------------
 *      c20*xi + c21*yi + c22
 *
 * Coefficients are calculated by solving linear system:
 * / x0 y0  1  0  0  0 -x0*u0 -y0*u0 \ /c00\ /u0\
 * | x1 y1  1  0  0  0 -x1*u1 -y1*u1 | |c01| |u1|
 * | x2 y2  1  0  0  0 -x2*u2 -y2*u2 | |c02| |u2|
 * | x3 y3  1  0  0  0 -x3*u3 -y3*u3 |.|c10|=|u3|,
 * |  0  0  0 x0 y0  1 -x0*v0 -y0*v0 | |c11| |v0|
 * |  0  0  0 x1 y1  1 -x1*v1 -y1*v1 | |c12| |v1|
 * |  0  0  0 x2 y2  1 -x2*v2 -y2*v2 | |c20| |v2|
 * \  0  0  0 x3 y3  1 -x3*v3 -y3*v3 / \c21/ \v3/
 *
 * where:
 *   cij - matrix coefficients, c22 = 1
 */
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
  M = Map<MatrixXd>(M.data(), 3, 3);
  return M;
}

void perspectiveTransform(vector<Point> & src, vector<Point> & dst, MatrixXd mtx) {

}

