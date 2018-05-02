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

#define TESTSNUM 6

int main(int argc, char const *argv[]) {
  vector<OptImgProc> imgs;
  for (int i = 0; i < TESTSNUM; ++i) {  
    imgs.push_back(OptImgProc("Dataset/" + std::to_string(i) + ".jpg"));
    // imgs[i].grayScale();
    // imgs[i].blurize(3);  
    // imgs[i].computeEdged();
    imgs[i].computeHough();
    // imgs[i].cimgComputeHough();
    imgs[i].computeResult();
    // imgs[i].computeResult();

    string save_path;
    // save_path = "result/" + std::to_string(i) + "-edge.jpg";
    // imgs[i].edged.save(save_path.c_str());
    save_path = "result/" + std::to_string(i) + "-hough.jpg";
    imgs[i].hough.normalize(0,255).save(save_path.c_str());
    save_path = "result/" + std::to_string(i) + "-res.jpg";
    imgs[i].result.save(save_path.c_str());
  }

  return 0;
}