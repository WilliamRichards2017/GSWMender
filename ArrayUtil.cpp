#include "ArrayUtil.h"
#include <iostream>
#include <vector>

using namespace std;

vector<vector<int> > ArrayUtil::buildArray2D(unsigned height, unsigned width){
  vector<vector<int> > array(height, vector<int>(width, 0));
  return array;
}

void ArrayUtil::printArray2D(vector<vector<int> > vec){
  for (int i = 0; i < vec.size(); i++){
    for (int j = 0; j < vec[i].size(); j++){
      cout << vec[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}


