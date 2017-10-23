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


//check that all the dimmensions of each correcsponding node in our traceback matches
void checkDimsMatch(vector<Traceback> tracebacks) {
  int n = 0;
  for(auto it = std::begin(tracebacks); it != std::end(tracebacks); ++it){
    assert(it->TBMs_.size() == it->MVMs_.size());
    assert((&it->TBMs_)[0].size() == (&it->MVMs_)[0].size());
    assert((&(&it->TBMs_)[0])[0].size() == (&(&it->MVMs_)[0])[0].size());
    n++;
  }
  cout << "Passed " << n << "/" << n << "dimension matching check tests\n";
}

void checkPos(vector<Variant> variants){
  int n = 0;
  for(auto it = std::begin(variants); it != std::end(variants); ++it){
    assert(it->ref[it->pos-1] == it->sv.first[0]);
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
  checkPos(variants);
  PileUp *p = new PileUp(variants);
  //vector<vector<vector<int> > > tbs = sumTBs(variants);
  //checkTracebackVectorSize(tbs);
  //checkTBMaxValue(tbs, variants);
  //printTracebacks(tbs);
}

int main(){
  runAllTests();
}
