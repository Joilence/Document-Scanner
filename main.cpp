// @Author: Jonathan Yang
// @Date: 2018.4.22

#include <CImg.h>
#include <stdlib.h>
#include <iostream>
#include "OptImgProc.cpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using namespace cimg_library;

int main(int argc, char const *argv[]) {
  vector<int> test_img;
  // test_img.push_back(0);
  // test_img.push_back(1);
  test_img.push_back(2);
  // test_img.push_back(3);
  // test_img.push_back(4);
  // test_img.push_back(5);
  // test_img.push_back(6);

  vector<OptImgProc> imgs;
  // for (int i = 0; i < TESTSNUM; ++i) {  
  for (int i = 0; i < test_img.size(); i++) {
    imgs.push_back(OptImgProc("Dataset/" + std::to_string(test_img[i]) + ".jpg"));
    // imgs[i].computeEdged();
    imgs[i].computeHough();
    imgs[i].computeResult();

    string save_path;
    // save_path = "result/" + std::to_string(test_img[i]) + "-edge.jpg";
    // imgs[i].edged.save(save_path.c_str());
    save_path = "result/" + std::to_string(test_img[i]) + "-hough.jpg";
    imgs[i].hough.normalize(0,255).save(save_path.c_str());
    save_path = "result/" + std::to_string(test_img[i]) + "-res.jpg";
    imgs[i].result.save(save_path.c_str());
    imgs[i].sort_intersections();
    imgs[i].computeTransMatrix();
    imgs[i].drawA4();
    save_path = "result/" + std::to_string(test_img[i]) + "-a4.jpg";
    imgs[i].A4.save(save_path.c_str());
  }

  return 0;
}