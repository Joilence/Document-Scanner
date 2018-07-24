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

  struct dirent *ptr;
  DIR *dir;
  string PATH = "./image/";
  dir = opendir(PATH.c_str());
  vector<string> test_img;
  while ((ptr = readdir(dir)) != NULL) {
    // skip '.' and '...'
    if (ptr->d_name[0] == '.') continue;
    // cout << ptr->d_name << endl;
    string filename = ptr->d_name;
    if (filename.find(".") >= filename.size()) {
      cout << filename << " ignored!\n";
      continue;
    }
    filename = filename.substr(0, filename.find("."));
    test_img.push_back(filename);
  }

  cout << test_img.size() << "images found." << endl;
  for(int i = 0; i < test_img.size(); ++i) {
    cout << test_img[i] << " ";
  }

  vector<OptImgProc> imgs;
  for (int i = 0; i < test_img.size(); i++) {
    // if (i == 3) break;
    cout << "####################\n";
    cout << "##### now on " << test_img[i] << " #####" << endl;
    cout << "####################\n";
    imgs.push_back(OptImgProc("image/" + test_img[i] + ".jpg"));
    // imgs[i].computeEdged();
    imgs[i].computeHough();
    imgs[i].computeResult();

    string save_path;
    save_path = "result/" + test_img[i] + "-edge.jpg";
    // imgs[i].edged.save(save_path.c_str());
    save_path = "result/" + test_img[i] + "-hough.jpg";
    // imgs[i].hough.normalize(0,255).save(save_path.c_str());
    save_path = "result/" + test_img[i] + "-res.jpg";
    // imgs[i].result.save(save_path.c_str());
    imgs[i].sort_intersections();
    imgs[i].computeTransMatrix();
    imgs[i].drawA4();
    save_path = "result/" + test_img[i] + "-a4.jpg";
    // imgs[i].A4.save(save_path.c_str());

    imgs[i].ostuSeg();
    // save_path = "seg/" + test_img[i] + "-ostu.jpg";
    // imgs[i].A4_ostuseg.save(save_path.c_str());

    imgs[i].iterSeg();
    save_path = "seg/" + test_img[i] + ".jpg";
    imgs[i].A4_iterseg.save(save_path.c_str());
  }

  return 0;
}