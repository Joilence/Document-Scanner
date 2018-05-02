#include <CImg.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>

#define GRADLIMIT 20
#define THRESHOLD 550
#define MAX_DIS 20
#define OVERLAP_DIS 20
#define MAX_PEAK 8
#define MIN_DOUBLE 0.0000001

using std::cout;
using std::endl;
using std::string;
using std::vector;
using namespace cimg_library;
using cimg::PI;


/**
 * Color
 */
const unsigned char white[3] = {255, 255, 255}, black[3] = {0, 0, 0}, red[] = {220, 20, 60};;

/**
 * @brief 2D Point
 * @details Stores coordinates and its votes
 */
struct Point {
  double x, y, value;
  Point(double _x, double _y, double _value) : x(_x), y(_y), value(_value) {}
};
bool comp(Point* a, Point* b) { return (*a).value > (*b).value; }

/**
 * @brief 2D Line
 * @details Stores k and b, and two points of it to reduce time complexity
 */
struct Line {
  double k, b;
  double x0, y0, x1, y1;
  Line(double _k, double _b, int _x0, int _y0, int _x1, int _y1)
      : k(_k), b(_b), x0(_x0), y0(_y0), x1(_x1), y1(_y1) {}
};

/**
 * @brief Calculate Euclidean distance
 */
double distance(double x, double y) { return sqrt(x * x + y * y); }

/**
 * @brief Image Processor
 * @details 
 */
class OptImgProc {
 public:
  CImg<double> src;
  CImg<double> edged;
  CImg<double> hough;
  CImg<double> result;
  CImg<double> cimg_edged;

  std::vector<Point*> peaks;
  std::vector<Line*> lines;
  std::vector<Point*> intersections;

  OptImgProc(string path) {
    src = CImg<double>(path.c_str());
    src.resize(600, 800);
    edged = CImg<double>(src.width(), src.height(), 1, 1, 0);
    int maxDistance = distance(src.width(), src.height());
    hough = CImg<double>(360, maxDistance, 1, 1, 0);
    result = CImg<double>(src);
  }

  void findPeaks() {
    cimg_forXY(hough, x, y) {
      if (peaks.size() > MAX_PEAK) break;
      double votes = hough(x, y);
      if (votes > THRESHOLD) {
        bool hasNeighbor = false;
        for (int i = 0; i < peaks.size(); ++i) {
          double dis = distance(peaks[i]->x - x, peaks[i]->y - y);
          if (dis < MAX_DIS) {
            hasNeighbor = true;
            if (peaks[i]->value < votes) peaks[i] = new Point(x, y, votes);
            break;
          }
        }
        if (!hasNeighbor) peaks.push_back(new Point(x, y, votes));
      }
    }
    std::sort(peaks.begin(), peaks.end(), comp);
    print_peaks();
  }

  CImg<double> computeEdged() {
    double thresL = 13.5, thresH = 13.6;

    const unsigned char black[3] = {255, 255, 255};
    edged = CImg<double>(src.width(), src.height());
    edged.fill(0);

    CImgList<> grad = src.blur(3).normalize(0, 255).get_gradient();
    cimg_forXY(src, s, t) {
      double
        Gs = grad[0](s, t),                   //
        Gt = grad[1](s, t),                   //  The actual pixel is (s,t)
        Gst = cimg::abs(Gs) + cimg::abs(Gt),  //
        // For efficient computation we observe that 
        // |Gs|+ |Gt| ~=~ sqrt(Gs^2 + Gt^2)
        Gr, Gur, Gu, Gul, Gl, Gdl, Gd, Gdr;
        // right, up right, up, up left, left, down left, down, down right.
      if (Gst >= thresH) {
        edged.draw_point(s, t, black);
      } else if (thresL <= Gst && Gst < thresH) {
        // Neighbourhood of the actual pixel:
        Gr = cimg::abs(grad[0](s + 1, t)) + cimg::abs(grad[1](s + 1, t));
        Gl = cimg::abs(grad[0](s - 1, t)) + cimg::abs(grad[1](s - 1, t));
        Gur =
            cimg::abs(grad[0](s + 1, t + 1)) + cimg::abs(grad[1](s + 1, t + 1));
        Gdl =
            cimg::abs(grad[0](s - 1, t - 1)) + cimg::abs(grad[1](s - 1, t - 1));
        Gu = cimg::abs(grad[0](s, t + 1)) + cimg::abs(grad[1](s, t + 1));
        Gd = cimg::abs(grad[0](s, t - 1)) + cimg::abs(grad[1](s, t - 1));
        Gul =
            cimg::abs(grad[0](s - 1, t + 1)) + cimg::abs(grad[1](s - 1, t + 1));
        Gdr =
            cimg::abs(grad[0](s + 1, t - 1)) + cimg::abs(grad[1](s + 1, t - 1));
        if (Gr >= thresH || Gur >= thresH || Gu >= thresH || Gul >= thresH ||
            Gl >= thresH || Gdl >= thresH || Gu >= thresH || Gdr >= thresH) {
          edged.draw_point(s, t, black);
        }
      };
    }
    return edged;
  }

  CImg<double> computeHough() {
    CImg<double> img =
        src.get_norm().normalize(0, 255).resize(-100, -100, 1, 2, 2);
    const double rhomax = std::sqrt((double)(img.width() * img.width() +
                                             img.height() * img.height())) /
                          2,
                 thetamax = 2 * cimg::PI;

    CImgList<> grad = img.get_gradient();
    cimglist_for(grad, l) grad[l].blur(1.5);
    CImg<double> vote(500, 400, 1, 1, 0);
    vote.fill(0);
    cimg_forXY(img, x, y) {
      const double X = (double)x - img.width() / 2,
                   Y = (double)y - img.height() / 2, gx = grad[0](x, y),
                   gy = grad[1](x, y);
      double theta = std::atan2(gy, gx),
             rho =
                 std::sqrt(X * X + Y * Y) * std::cos(std::atan2(Y, X) - theta);
      if (rho < 0) {
        rho = -rho;
        theta += cimg::PI;
      }
      theta = cimg::mod(theta, thetamax);
      vote((int)(theta * vote.width() / thetamax),
           (int)(rho * vote.height() / rhomax)) +=
          (float)std::sqrt(gx * gx + gy * gy);
    }
    vote.blur(0.5);
    CImg<> vote2(vote);
    cimg_forXY(vote2, x, y) vote2(x, y) = (float)std::log(1 + vote(x, y));
    hough = CImg<double>(vote2).normalize(0, 1000);
    // hough.display();
    return hough;
  }

  bool isInImage(double cx, double cy) {
    return cx <= src.width() && cx >= 0 && cy >= 0 && cy <= src.height();
  }

  bool isNeighbor(double k, double b, double ki, double bi) {
    double cu = ((double)(0 - b) / k), 
           cd = ((double)(src.height() - b) / k),
           cr = ((double)(k * src.width() + b)),
           cl = b;
    double cui = ((double)(0 - bi) / ki),
           cdi = ((double)(src.height() - bi) / ki),
           cri = ((double)(ki * src.width() + bi)), 
           cli = bi;

    if (isInImage(cu, cui) && abs(cui - cu) < OVERLAP_DIS) return true;
    if (isInImage(cd, cdi) && abs(cdi - cd) < OVERLAP_DIS) return true;
    if (isInImage(cr, cri) && abs(cri - cr) < OVERLAP_DIS) return true;
    if (isInImage(cl, cli) && abs(cli - cl) < OVERLAP_DIS) return true;
    return false;
  }

  void findLines() {
    double rhomax = std::sqrt((double)(src.width() * src.width() +
                                        src.height() * src.height())) /
                              2,
           thetamax = 2 * cimg::PI;
    for (int i = 0; i < peaks.size(); ++i) {
      double rho = peaks[i]->y * rhomax / hough.height(),
             theta = peaks[i]->x * thetamax / hough.width(),
             x = result.width() / 2 + rho * std::cos(theta),
             y = result.height() / 2 + rho * std::sin(theta);
      double x0 = (double)(x + 1000 * std::sin(theta)),
             y0 = (double)(y - 1000 * std::cos(theta)),
             x1 = (double)(x - 1000 * std::sin(theta)),
             y1 = (double)(y + 1000 * std::cos(theta));
      // cout << x0 << " " << y0 << endl;
      // cout << x1 << " " << y1 << endl;
      double k = (double)((y1 - y0) / (x1 - x0 + MIN_DOUBLE));
      double b = (double)(y0 - k * x0);

      // sieve of lines
      cout << "---\n";
      cout << "Judging ... \n";
      cout << "y = " << k << " * x + " << b << endl;
      bool hasNeighbor = false;
      int parallel = 0;
      for (int i = 0; i < lines.size(); ++i) {
        double ki = lines[i]->k;
        double tan = abs(ki - k) / abs(1 + ki * k);
        if (tan < 1.5) {
          if (isNeighbor(k, b, lines[i]->k, lines[i]->b)) {
            hasNeighbor = true;
            cout << "Neighbourhood found.\n";
            cout << "y = " << lines[i]->k << " * x + " << lines[i]->b << endl;
            break;
          } else {
            parallel += 1;
            cout << "tan = " << tan << ", parallel found " << parallel << ".\n";
          }
        }
        if (parallel >= 2) {
          break;
        }
      }
      cout << "---\n";

      if (!hasNeighbor && parallel < 2)
        lines.push_back(new Line(k, b, x0, y0, x1, y1));
    }
    print_lines();
  }

  void computeIntersections() {
    for (int i = 0; i < lines.size(); ++i) {
      for (int j = i + 1; j < lines.size(); ++j) {
        double ki = lines[i]->k;
        double bi = lines[i]->b;
        double kj = lines[j]->k;
        double bj = lines[j]->b;

        double x = (bi - bj) / (kj - ki);
        double y = (ki * x + bi);
        intersections.push_back(new Point(x, y, 0));
      }
    }
  }

  void computeResult() {
    result = CImg<double>(src);
    findPeaks();
    findLines();
    computeIntersections();

    for (int i = 0; i < lines.size(); ++i) {
      double x0 = lines[i]->x0;
      double y0 = lines[i]->y0;
      double x1 = lines[i]->x1;
      double y1 = lines[i]->y1;
      // cout << x0 << " " << y0 << endl;
      // cout << x1 << " " << y1 << endl;
      result.draw_line(x0, y0, x1, y1, red, 1.0f, 0xF0F0F0F0);
          // .draw_line(x0, y0, x1, y1, black, 1.0f, 0x0F0F0F0F)
          // .draw_line(x0 + 1, y0, x1 + 1, y1, white, 1.0f, 0xF0F0F0F0)
          // .draw_line(x0 + 1, y0, x1 + 1, y1, black, 1.0f, 0x0F0F0F0F)
          // .draw_line(x0, y0 + 1, x1, y1 + 1, white, 1.0f, 0xF0F0F0F0)
          // .draw_line(x0, y0 + 1, x1, y1 + 1, black, 1.0f, 0x0F0F0F0F);
    }

    int rad = (src.height() + src.width()) / 200;
    for (int i = 0; i < intersections.size(); ++i) {
      result.draw_circle(intersections[i]->x, intersections[i]->y, rad, red);
    }

    // result.display();
  }

  void allProcess() {
    // computeEdged();
    computeHough();
    computeResult();
  }

  void print_peaks() {
    for (int i = 0; i < peaks.size(); ++i) {
      cout << "("
           << "【" << peaks[i]->value << "】" << peaks[i]->x << ", "
           << peaks[i]->y << ")"
           << " ";
    }
    cout << endl;
  }

  void print_lines() {
    cout << "Lines:\n";
    for (int i = 0; i < lines.size(); ++i) {
      cout << "y = " << lines[i]->k << " * x + " << lines[i]->b << endl;
    }
  }
};