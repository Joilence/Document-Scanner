#include <fstream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/ml.hpp>
#include "opencv2/imgcodecs.hpp"

using namespace cv;
using namespace cv::ml;
using namespace std;

//大端转小端
int reverseInt(int i);

int main() {
  //读取测试样本集
  ifstream if_testImags("../../data/t10k-images-idx3-ubyte", ios::binary);
  //读取失败
  if (true == if_testImags.fail()) {
    cout << "Please check the path of file t10k-images-idx3-ubyte" << endl;
    return 1;
  }
  int magic_num, testImgsNum, nrows, ncols;
  //读取magic number
  if_testImags.read((char*)&magic_num, sizeof(magic_num));
  magic_num = reverseInt(magic_num);
  cout << "测试图像数据库t10k-images-idx3-ubyte的magic number为：" << magic_num
       << endl;
  //读取测试图像总数
  if_testImags.read((char*)&testImgsNum, sizeof(testImgsNum));
  testImgsNum = reverseInt(testImgsNum);
  cout << "测试图像数据库t10k-images-idx3-ubyte的图像总数为：" << testImgsNum
       << endl;
  //读取图像的行大小
  if_testImags.read((char*)&nrows, sizeof(nrows));
  nrows = reverseInt(nrows);
  cout << "测试图像数据库t10k-images-idx3-ubyte的图像维度row为：" << nrows
       << endl;
  //读取图像的列大小
  if_testImags.read((char*)&ncols, sizeof(ncols));
  ncols = reverseInt(ncols);
  cout << "测试图像数据库t10k-images-idx3-ubyte的图像维度col为：" << ncols
       << endl;

  //读取测试图像
  int imgVectorLen = nrows * ncols;
  Mat testFeatures = Mat::zeros(testImgsNum, imgVectorLen, CV_32FC1);
  Mat temp = Mat::zeros(nrows, ncols, CV_8UC1);
  for (int i = 0; i < testImgsNum; i++) {
    if_testImags.read((char*)temp.data, imgVectorLen);
    Mat tempFloat;
    //由于SVM需要的测试数据格式是CV_32FC1，在这里进行转换
    temp.convertTo(tempFloat, CV_32FC1);
    memcpy(testFeatures.data + i * imgVectorLen * sizeof(float), tempFloat.data,
           imgVectorLen * sizeof(float));
  }
  //归一化
  testFeatures = testFeatures / 255;
  //读取测试图像对应的分类标签
  ifstream if_testLabels("t10k-labels-idx1-ubyte", ios::binary);
  //读取失败
  if (true == if_testLabels.fail()) {
    cout << "Please check the path of file t10k-labels-idx1-ubyte" << endl;
    return 1;
  }
  int magic_num_2, testLblsNum;
  //读取magic number
  if_testLabels.read((char*)&magic_num_2, sizeof(magic_num_2));
  magic_num_2 = reverseInt(magic_num_2);
  cout << "测试图像标签数据库t10k-labels-idx1-ubyte的magic number为："
       << magic_num_2 << endl;
  //读取测试图像的分类标签的数量
  if_testLabels.read((char*)&testLblsNum, sizeof(testLblsNum));
  testLblsNum = reverseInt(testLblsNum);
  cout << "测试图像标签数据库t10k-labels-idx1-ubyte的标签总数为："
       << testLblsNum << endl;

  //由于SVM需要输入的标签类型是CV_32SC1，在这里进行转换
  Mat testLabels = Mat::zeros(testLblsNum, 1, CV_32SC1);
  Mat readLabels = Mat::zeros(testLblsNum, 1, CV_8UC1);
  if_testLabels.read((char*)readLabels.data, testLblsNum * sizeof(char));
  readLabels.convertTo(testLabels, CV_32SC1);

  //载入训练好的SVM模型
  Ptr<SVM> svm = SVM::load("mnist.xml");
  int sum = 0;
  //对每一个测试图像进行SVM分类预测
  for (int i = 0; i < testLblsNum; i++) {
    Mat predict_mat = Mat::zeros(1, imgVectorLen, CV_32FC1);
    memcpy(predict_mat.data,
           testFeatures.data + i * imgVectorLen * sizeof(float),
           imgVectorLen * sizeof(float));
    //预测
    float predict_label = svm->predict(predict_mat);
    //真实的样本标签
    float truth_label = testLabels.at<int>(i);
    //比较判定是否预测正确
    if ((int)predict_label == (int)truth_label) {
      sum++;
    }
  }

  cout << "预测准确率为：" << (double)sum / (double)testLblsNum << endl;

  //随机测试某一个图像看效果，输入为-2时退出，输入-1时则测试本地图片“2.jpg”，注意路径要放到源代码同级目录
  while (1) {
    int index;
    cout << "请输入要查看的测试图像下标" << endl;
    cin >> index;
    if (-1 == index) {
      Mat imgRead = imread("2.jpg", 0);
      Mat imgReadScal = Mat::zeros(nrows, ncols, CV_8UC1);
      Mat show_mat = Mat::zeros(nrows, ncols, CV_32FC1);

      resize(imgRead, imgReadScal, imgReadScal.size());

      imgReadScal.convertTo(show_mat, CV_32FC1);

      show_mat = show_mat / 255;

      //
      Mat predict_mat = Mat::zeros(1, imgVectorLen, CV_32FC1);
      // memcpy(show_mat.data, testFeatures.data + index*imgVectorLen *
      // sizeof(float), imgVectorLen * sizeof(float));
      memcpy(predict_mat.data, show_mat.data, imgVectorLen * sizeof(float));
      float response = svm->predict(predict_mat);

      imshow("test", show_mat);
      cout << "标签值为" << response << endl;

      waitKey(0);

    } else if (-2 == index) {
      break;
    } else {
      Mat show_mat = Mat::zeros(nrows, ncols, CV_32FC1);
      Mat predict_mat = Mat::zeros(1, imgVectorLen, CV_32FC1);
      memcpy(show_mat.data,
             testFeatures.data + index * imgVectorLen * sizeof(float),
             imgVectorLen * sizeof(float));
      memcpy(predict_mat.data,
             testFeatures.data + index * imgVectorLen * sizeof(float),
             imgVectorLen * sizeof(float));
      float response = svm->predict(predict_mat);

      imshow("test", show_mat);
      cout << "标签值为" << response << endl;

      waitKey(0);
    }
  }
  return 0;
}

//大端转小端
int reverseInt(int i)
{
    unsigned char c1, c2, c3, c4;

    c1 = i & 255;
    c2 = (i >> 8) & 255;
    c3 = (i >> 16) & 255;
    c4 = (i >> 24) & 255;

    return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;
}