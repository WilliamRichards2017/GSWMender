#include "Traceback.h"

vector<vector<int> > buildArray2D(unsigned height, unsigned width){
  vector<vector<int> > array(height, std::vector<int>(width, 0));
  return array;
}

// get the x,y coordinates of the maximum value in Score matrix
pair<int, int>  Traceback::getMaxCoords(){
  int maxV = -1;
  int x = -1;
  int y = -1;

  for (int i = 0; i < MVM_.size(); i++){
    for (int j = 0; j < MVM_[i].size(); j++){
      if(MVM_[i][j] > maxV){
	y = i;
	x = j;
	maxV = MVM_[i][j];
      }
    }
  }
  pair<int, int> coords = std::make_pair(x,y);
  return coords;
}

vector<vector<vector<int> > > Traceback::buildTBMs(){
  map<Node *, vector< vector< vector<int> > >, less<Node *> > GS = ga_->getScoreMatrix();
  int l2 = ga_->getQueryLength();
  for (auto it = std::begin(subjectNodes_); it != std::end(subjectNodes_); ++it){
    int l1 = (*it)->getSequence().length();
    MVM_ = buildArray2D(l2+1,l1+1);
    vector<vector<int> > TBM = buildArray2D(l2+1,l1+1);
    vector< vector< vector<int> > > S = GS[*it];
    for (int i2=0; i2<=l2; i2++) {
      for (int i1=0; i1<=l1; i1++) {
	MVM_[i2][i1] = max(max(S[i1][i2][1],S[i1][i2][2]),S[i1][i2][0]);
      }
    }
    std::pair<int, int> coords = getMaxCoords();
    int x = coords.first;
    int y = coords.second;
    //start at max value coords
    while(x > -1 && y > -1){
      TBM[y][x] = 1;
      //Move diagonal
      if(x == 0 || y == 0){
	break;
      }
      else if(MVM_[y-1][x-1] >= MVM_[y-1][x] && MVM_[y-1][x-1] >= MVM_[y][x-1]){
	y--;
	x--;
      }
      //Move up
      else if(MVM_[y-1][x] > MVM_[y-1][x-1] && MVM_[y-1][x] > MVM_[y][x-1]){
	y--;
      }
      //move left
      else{
	x--;
      }
    } // end of while
    TBMs_.push_back(TBM);
  }// end of node loop
}

