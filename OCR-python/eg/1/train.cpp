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
  //读取训练样本集
  ifstream if_trainImags("../../data/train-images-idx3-ubyte", ios::binary);
  //读取失败
  if (true == if_trainImags.fail()) {
    cout << "Please check the path of file train-images-idx3-ubyte" << endl;
    return 1;
  }
  int magic_num, trainImgsNum, nrows, ncols;
  //读取magic number
  if_trainImags.read((char*)&magic_num, sizeof(magic_num));
  magic_num = reverseInt(magic_num);
  cout << "训练图像数据库train-images-idx3-ubyte的magic number为：" << magic_num
       << endl;
  //读取训练图像总数
  if_trainImags.read((char*)&trainImgsNum, sizeof(trainImgsNum));
  trainImgsNum = reverseInt(trainImgsNum);
  cout << "训练图像数据库train-images-idx3-ubyte的图像总数为：" << trainImgsNum
       << endl;
  //读取图像的行大小
  if_trainImags.read((char*)&nrows, sizeof(nrows));
  nrows = reverseInt(nrows);
  cout << "训练图像数据库train-images-idx3-ubyte的图像维度row为：" << nrows
       << endl;
  //读取图像的列大小
  if_trainImags.read((char*)&ncols, sizeof(ncols));
  ncols = reverseInt(ncols);
  cout << "训练图像数据库train-images-idx3-ubyte的图像维度col为：" << ncols
       << endl;

  //读取训练图像
  int imgVectorLen = nrows * ncols;
  Mat trainFeatures = Mat::zeros(trainImgsNum, imgVectorLen, CV_32FC1);
  Mat temp = Mat::zeros(nrows, ncols, CV_8UC1);
  for (int i = 0; i < trainImgsNum; i++) {
    if_trainImags.read((char*)temp.data, imgVectorLen);
    Mat tempFloat;
    //由于SVM需要的训练数据格式是CV_32FC1，在这里进行转换
    temp.convertTo(tempFloat, CV_32FC1);
    memcpy(trainFeatures.data + i * imgVectorLen * sizeof(float),
           tempFloat.data, imgVectorLen * sizeof(float));
  }
  //归一化
  trainFeatures = trainFeatures / 255;
  //读取训练图像对应的分类标签
  ifstream if_trainLabels("train-labels-idx1-ubyte", ios::binary);
  //读取失败
  if (true == if_trainLabels.fail()) {
    cout << "Please check the path of file train-labels-idx1-ubyte" << endl;
    return 1;
  }
  int magic_num_2, trainLblsNum;
  //读取magic number
  if_trainLabels.read((char*)&magic_num_2, sizeof(magic_num_2));
  magic_num_2 = reverseInt(magic_num_2);
  cout << "训练图像标签数据库train-labels-idx1-ubyte的magic number为："
       << magic_num_2 << endl;
  //读取训练图像的分类标签的数量
  if_trainLabels.read((char*)&trainLblsNum, sizeof(trainLblsNum));
  trainLblsNum = reverseInt(trainLblsNum);
  cout << "训练图像标签数据库train-labels-idx1-ubyte的标签总数为："
       << trainLblsNum << endl;

  //由于SVM需要输入的标签类型是CV_32SC1，在这里进行转换
  Mat trainLabels = Mat::zeros(trainLblsNum, 1, CV_32SC1);
  Mat readLabels = Mat::zeros(trainLblsNum, 1, CV_8UC1);
  if_trainLabels.read((char*)readLabels.data, trainLblsNum * sizeof(char));
  readLabels.convertTo(trainLabels, CV_32SC1);

  /*
  //Add some random test code
  while (1)
  {
  int index;
  cout << "请输入要查看的训练图像下标" << endl;
  cin >> index;
  if (-1 == index)
  {
  break;
  }
  Mat show_mat = Mat::zeros(nrows, ncols, CV_32FC1);
  memcpy(show_mat.data, trainFeatures.data + index*imgVectorLen * sizeof(float),
  imgVectorLen * sizeof(float)); Mat show_char; show_mat.convertTo(show_char,
  CV_8UC1); imshow("test", show_mat); cout << "标签值为" <<
  trainLabels.at<int>(index);
  }
  waitKey(0);
  */

  // 训练SVM分类器
  //初始化
  Ptr<SVM> svm = SVM::create();
  //多分类
  svm->setType(SVM::C_SVC);
  // kernal选用RBF
  svm->setKernel(SVM::RBF);
  //设置经验值
  svm->setGamma(0.01);
  svm->setC(10.0);
  //设置终止条件，在这里选择迭代200次
  svm->setTermCriteria(TermCriteria(TermCriteria::MAX_ITER, 200, FLT_EPSILON));
  //训练开始
  svm->train(trainFeatures, ROW_SAMPLE, trainLabels);

  cout << "训练结束，正写入xml:" << endl;
  //保存模型
  svm->save("mnist.xml");



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