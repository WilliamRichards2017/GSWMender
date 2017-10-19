#include "gsw.h"
#include "Variant.h"
#include "Traceback.h"
#include "PileUp.h"
#include "Class-GraphAlignment.h"
#include "ArrayUtil.h"

#include <vector>
#include <string>
#include <assert.h>

struct Variants{
  vector<Variant> variants_;
};

vector<Traceback> buildTracebackVector(vector<Variant> variants){
  vector<Traceback> tracebackVec;

  for(auto it = std::begin(variants); it != std::end(variants); ++it){
    PileUp * p = new PileUp();
    vector<Node *> subjectNodes = p->buildDiamondGraph(p->getNodes(*it));
    GraphAlignment * ga = new GraphAlignment::GraphAlignment(subjectNodes, it->ref, 2, -2, -3, -2, false);
    Traceback t = {subjectNodes, ga};
    tracebackVec.push_back(t);
  }
  return tracebackVec;
}

vector<Variant> buildAllVariants(){

  string query1 = "TAAAGCTGTGTTGTGCT";
  string query2 = "TAAAGCTGTTGTGCT";
  string query3 = "TAAAGCGTGTGTGCT";
  string query4 = "TAAAGCGTGTTGTGCT";

  std::pair<string,string> sv1 = std::make_pair("CGATTGTTT","TGTGT");
  std::pair<string,string> sv2 = std::make_pair("CGATTGTTT", "TGT");
  std::pair<string,string> sv3 = std::make_pair("CGATTGTTT", "GTG");
  std::pair<string,string> sv4 = std::make_pair("CGATTGTTT", "GTGT");

  int pos = 6;

  Variant v1 = {query1, sv1, pos};
  Variant v2 = {query2, sv2, pos};
  Variant v3 = {query3, sv3, pos};
  Variant v4 = {query4, sv4, pos};


  vector<Variant> variants;
  variants.push_back(v1);
  variants.push_back(v2);
  variants.push_back(v3);
  variants.push_back(v4);

  return variants;
}

vector<vector<vector<int> > >sumTBs(vector<Variant> variants){
  vector<Traceback> tracebacks = buildTracebackVector(variants);
  PileUp p = {tracebacks};
  return p.sumTracebacks();

}

void check_bounds() {

}

void checkPos(vector<Variant> variants){
  for(auto it = std::begin(variants); it != std::end(variants); ++it){
    assert(it->ref[it->pos-1] == it->sv.first[0]);
  }
  cout << "Passed all checks in checkPos\n";
}

void checkTracebackVectorSize(vector<vector<vector<int> > > tbs){
  assert(tbs.size()==4);
  cout << "Passing correct TB vector length\n";
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
    int z = 0;
    testMaxValue((*it), maxV);
  }
  cout << "passed TBM minimum and maximum bounds of values\n";
}


void runAllTests(){
  static int D = 0;
  static int H = 1;
  static int V = 2;

  vector<Variant> variants = buildAllVariants();
  checkPos(variants);
  vector<vector<vector<int> > > tbs = sumTBs(variants);
  checkTracebackVectorSize(tbs);
  int z = 0;
  for(auto it = std::begin(tbs); it != std::end(tbs); ++it){
    ArrayUtil::printArray2D(tbs[z]);
    z++;
  }
  checkTBMaxValue(tbs, variants);
}

int main(){
  runAllTests();
}
