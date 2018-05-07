#include <CImg.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>

#define GRADLIMIT 20
#define THRESHOLD 550         // baseline of peaks value
#define NEIGHBOR_PEAK_DIS 20  // decide whether two peaks are neighbor
#define NEIGHBOR_LINE_DIS 20  // decide whether two lines are neighbor
#define MAX_PEAK 8
#define MIN_DOUBLE 0.0000001  // make every denominator be non-zero

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
  vector<Point> intersections;
  vector<Point> borderPoints;
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

  std::vector<Point*> peaks;
  std::vector<Line*> lines;
  std::vector<Point*> intersections;

  OptImgProc(string path) {
    src = CImg<double>(path.c_str());
    edged = CImg<double>(src.width(), src.height(), 1, 1, 0);
    result = CImg<double>(src);
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
    grad[0].display();
    grad[1].display();
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
    hough.display();
    return hough;
  }

    void findPeaks() {
    cimg_forXY(hough, x, y) {
      double votes = hough(x, y);
      if (votes > THRESHOLD) {
        bool hasNeighbor = false;
        for (int i = 0; i < peaks.size(); ++i) {
          double dis = distance(peaks[i]->x - x, peaks[i]->y - y);
          if (dis < NEIGHBOR_PEAK_DIS) {
            hasNeighbor = true;
            if (peaks[i]->value < votes) peaks[i] = new Point(x, y, votes);
            break;
          }
        }

        // if has neighbor, then the point was either picked or passed
        // if not, continue
        if (!hasNeighbor) {
          // if peaks are not full, push back directly
          if (peaks.size() < MAX_PEAK) { //
            peaks.push_back(new Point(x, y, votes));
          // if peaks is full, compare with the last peak
          } else if (peaks.back()->value < votes) {
            peaks.back()->value = votes;
            peaks.back()->x = x;
            peaks.back()->y = y;
          }
        }
        std::sort(peaks.begin(), peaks.end(), comp);
      }
    }
    std::sort(peaks.begin(), peaks.end(), comp);
    print_peaks();
  }

  bool isInImage(double cx, double cy) {
    return cx <= src.width() && cx >= 0 && cy >= 0 && cy <= src.height();
  }

  bool isNeighbor(double k, double b, double ki, double bi) {
    // u = up, d = down, r = right, l = left
    double cu = ((double)(0 - b) / k), 
           cd = ((double)(src.height() - b) / k),
           cr = ((double)(k * src.width() + b)),
           cl = b;
    double cui = ((double)(0 - bi) / ki),
           cdi = ((double)(src.height() - bi) / ki),
           cri = ((double)(ki * src.width() + bi)), 
           cli = bi;

    if (isInImage(cu, 0) && isInImage(cui, 0) &&
        abs(cui - cu) < NEIGHBOR_LINE_DIS)
      return true;
    if (isInImage(cd, src.height()) && isInImage(cdi, src.height()) &&
        abs(cdi - cd) < NEIGHBOR_LINE_DIS)
      return true;
    if (isInImage(src.width(), cr) && isInImage(src.width(), cri) &&
        abs(cri - cr) < NEIGHBOR_LINE_DIS)
      return true;
    if (isInImage(0, cl) && isInImage(0, cli) &&
        abs(cli - cl) < NEIGHBOR_LINE_DIS)
      return true;
    return false;
  }

  void findLines() {
    double rhomax = std::sqrt((double)(src.width() * src.width() +
                                        src.height() * src.height())) /
                              2,
           thetamax = 2 * cimg::PI;

    // for every peak, compute a line
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
      cout << line_equation(k, b) << endl;
      bool hasNeighbor = false;
      int parallel = 0;
      for (int i = 0; i < lines.size(); ++i) {
        double ki = lines[i]->k;
        double tan = abs(ki - k) / abs(1 + ki * k);
        if (tan < 1.5) {
          if (isNeighbor(k, b, lines[i]->k, lines[i]->b)) {
            hasNeighbor = true;
            cout << "Neighbourhood found.\n";
            // cout << "y = " << lines[i]->k << " * x + " << lines[i]->b << endl;
            cout << line_equation(lines[i]->k, lines[i]->b) << endl;
            break;
          } else {
            parallel += 1;
            cout << "tan = " << tan << ", parallel found " << parallel << ".\n";
            cout << line_equation(lines[i]->k, lines[i]->b) << endl;
          }
        }
        if (parallel >= 2) {
          break;
        }
      }
      cout << "---\n";

      // if (!hasNeighbor && parallel < 2)
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
        if (isInImage(x, y)) {
          intersections.push_back(new Point(x, y, 0));
          lines[i]->intersections.push_back(Point(x, y, 0));
          lines[j]->intersections.push_back(Point(x, y, 0));
        }
      }
    }
  }

  void computeResult(bool debug = false) {
    result = CImg<double>(src);
    findPeaks();
    findLines();
    computeIntersections();

    drawResult(debug);

    int rad = (src.height() + src.width()) / 200;
    for (int i = 0; i < intersections.size(); ++i) {
      result.draw_circle(intersections[i]->x, intersections[i]->y, rad, red);
    }

    // result.display();
  }

  void drawResult (bool debug = false) {
    for (int i = 0; i < lines.size(); ++i) {
      double k = lines[i]->k;
      double b = lines[i]->b;

      vector<Point> linePoint;
      if (debug) {
        double cu = ((double)(0 - b) / k),
               cd = ((double)(src.height() - b) / k),
               cr = ((double)(k * src.width() + b)), cl = b;
        int count = 0;
        if (isInImage(cu, 0)) {
          linePoint.push_back(Point(cu, 0, 0));
          count++;
        }
        if (isInImage(cd, src.height())) {
          linePoint.push_back(Point(cd, src.height(), 0));
          count++;
        }
        if (count < 2 && isInImage(src.width(), cr)) {
          linePoint.push_back(Point(src.width(), cr, 0));
          count++;
        }
        if (count < 2 && isInImage(0, cl)) {
          linePoint.push_back(Point(0, cl, 0));
        }
      } else {
        linePoint.push_back(lines[i]->intersections[0]);
        linePoint.push_back(lines[i]->intersections[1]);
      }

      double x0 = linePoint[0].x;
      double y0 = linePoint[0].y;
      double x1 = linePoint[1].x;
      double y1 = linePoint[1].y;

      // cout << x0 << " " << y0 << endl;
      // cout << x1 << " " << y1 << endl;
      result.draw_line(x0, y0, x1, y1, white, 1.0f, 0xF0F0F0F0)
          .draw_line(x0, y0, x1, y1, black, 1.0f, 0x0F0F0F0F);

      // Widen the line
      int width = 7;
      for (int i = 1; i <= width; ++i) {
        result.draw_line(x0 + i, y0, x1 + i, y1, white, 1.0f, 0xF0F0F0F0)
          .draw_line(x0 + i, y0, x1 + i, y1, black, 1.0f, 0x0F0F0F0F)
          .draw_line(x0, y0 + i, x1, y1 + i, white, 1.0f, 0xF0F0F0F0)
          .draw_line(x0, y0 + i, x1, y1 + i, black, 1.0f, 0x0F0F0F0F);
      }

      // compute mid-point
      double mx = (x0 + x1) / 2;
      double my = (y0 + y1) / 2;
      // add annotation of line
      int font_height = 60;
      string line_info = line_equation(k, b);
      result.draw_text(mx, my, line_info.c_str(), white, black, 1, font_height);
    }
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
      // cout << "y = " << lines[i]->k << " * x + " << lines[i]->b << endl;
      cout << line_equation(lines[i]->k, lines[i]->b) << endl;
    }
  }

  string line_equation(double k, double b) {
    return "y = " + std::to_string(k) + " * x + " + std::to_string(b);
  }
};