#include "gsw.h"
#include "Variant.h"
#include "Traceback.h"
#include "PileUp.h"
#include "Class-GraphAlignment.h"
#include "ArrayUtil.h"
#include <vector>
#include <string>
#include <assert.h>

vector<Variant> buildAllVariants(){

  string query1 = "CTTTCTTTCTTTCCTTCCTT";
  string query2 = "CTTTCTTTCTTTTCCTTT";
  string query3 = "TTTCTTTCTTTCTTTTCCTTT";


  std::pair<string,string> sv = std::make_pair("TT","TTCCT");

  int pos1 = 10;
  int pos2 = 10;
  int pos3 = 13;

  Variant v1 = {query1, pos1};
  Variant v2 = {query2, pos2};
  Variant v3 = {query3, pos3};
  // Variant v4 = {query4, pos4};


  vector<Variant> variants;
  variants.push_back(v1);
  variants.push_back(v2);
  variants.push_back(v3);
  //variants.push_back(v4);

  return variants;
}


//check that all the dimmensions of each correcsponding node in our traceback matches
void checkDimsMatch(vector<Traceback> tracebacks) {
  int n = 0;
  for(auto it = std::begin(tracebacks); it != std::end(tracebacks); ++it){
    assert(it->TBMs_.size() == it->MVMs_.size());
    assert((&it->TBMs_)[0].size() == (&it->MVMs_)[0].size());
    assert((&(&it->TBMs_)[0])[0].size() == (&(&it->MVMs_)[0])[0].size());
    n++;
  }
  cout << "Passed " << n << "/" << n << " dimension matching check tests\n";
}

void checkPos(vector<Variant> variants, std::pair<string,string> sv){
  int n = 0;
  for(auto it = std::begin(variants); it != std::end(variants); ++it){
    assert(it->ref[it->pos-1] == sv.first[0]);
    n++;
  }
  cout << "Passed " << n << "/" << n << " breakpoint position tests\n";
}

void checkTracebackVectorSize(vector<vector<vector<int> > > tbs){
  assert(tbs.size()==4);
  cout << "Passed 1/1 correct TB vector length\n";
}

void testMaxValue(vector<vector<int> > tbm, int m){
  for(int i = 0; i < tbm.size(); i++){
    for(int j = 0; j < tbm[i].size(); j++){
      assert(tbm[i][j] <= m && tbm[i][j] >= 0);
    }
  }
}

void checkTBMaxValue(vector<vector<vector<int> > > tbs, vector<Variant> variants){
  int maxV = variants.size();
  int c = 0;
  for(auto it = std::begin(tbs); it != std::end(tbs); ++it){
    testMaxValue((*it), maxV);
    c++;
  }
  cout << "Passed " << c << "/" << c << " traceback value bound tests\n";
}

void checkGraphEquality(vector<Traceback> tbs){
  assert(tbs.size() > 0);
  if(tbs.size() == 1){
    cout << "Passed all graph equality checks\n";
    return;
  }
  else{
    Traceback t1 = tbs.back();
    tbs.pop_back();
    Traceback t2 = tbs.back();
    
    vector<Node *> sn1 = t1.subjectNodes_;
    vector<Node *> sn2 = t2.subjectNodes_;

    assert(sn1.size() == 4);
    assert(sn2.size() == 4);

    for(unsigned i = 0; i < 4; i++){
      assert(sn1[i]->getSequence() == sn2[i]->getSequence());
    }
    checkGraphEquality(tbs);
  }
  return;
}


void printTracebacks(vector<vector<vector<int> > > tbs){
  int z = 0;
  cout << "\n";
  for(auto it = std::begin(tbs); it != std::end(tbs); ++it){
    cout << "Print out node: " << z << std::endl;
    ArrayUtil::printArray2D(tbs[z]);
    z++;
  }
}


void runAllTests(){
  static int D = 0;
  static int H = 1;
  static int V = 2;
  vector<Variant> variants = buildAllVariants();
  std::pair<string,string> sv = std::make_pair("TT","TTCCT");
  checkPos(variants, sv);
  Variant ref = {"CTTTCTTTCTTTCCTTCCTT", 10};
  PileUp *p = new PileUp(variants, sv, ref);
  checkGraphEquality(p->tbs_);
  checkTracebackVectorSize(p->sumMatrix_);
  checkTBMaxValue(p->sumMatrix_, variants);
  checkDimsMatch(p->tbs_);
  printTracebacks(p->sumMatrix_);
}

int main(){
  runAllTests();
}
