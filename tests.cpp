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
  int refPos = 184258400;
  //int refPos= 184258390;
  string query1 = "GGATGGCACAGACATCCACTTTAGCAAAGTGGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCTTGCAGCAGAGTTGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCCTTCACCATCTAAGGAA"; //d
  string query2 = "GATGGCACAGACATCCACTTTAGCAAAGTGGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCTTGCAGCAGAGTTGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCCTTCACCATCTAAGGAAAGG"; //d
  string query3 = "GGCACAGACATCCACTTTAGCAAAGTGGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCTTGCAGCAGAGTTGCAAGAGGTCTTGGGACCTGTGGTCCTAATGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGAC"; //d
  string query4 = "CACAGACATCCACTTTAGCATTGTGGTGTGCCTACTTGGATCAGTAGATCACACTGAGTTCAAGTCTCCTTGCAGCTGAGTTGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCCTTCACCATCTATGGTTATGTGTCC";  //d
  string query5 = "CAGACATCCACTTTAGCAAAGTGGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCGTGCAGCAGAGTTGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCCCTCACCATCTAAGGAAAGGTGTCCCC"; //d
  string query6 = "CATCCACTTTAGCAAAGTGGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCTTGCAGCAGAGTTGCAAGAGGTCTTGGGACCTGTGGTCCTAATGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCC";
  string query7 = "CCACTTTAGCAAAGTGGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCTTGCAGCAGAGGTGCAAGAGGTCTTGGGACCTGTGGTCCTAATGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCCTTC";
  string query8 = "ACTTTAGCAAAGTGGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCTTGCAGCAGAGTTGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCCTTCACCATCTAAGGAAAGGTGTCCCCCCATTCCCA";
  string query9 = "CTTTAGCAAAGTGGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCTTGCAGCAGAGTTGCAAGAGGGCTTGGGACCTGTGGTCCTAATGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCCTTCACC";
  string query10 = "TTAGCAAAGTGGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCGTGCAGCAGAGTTGCAAGAGGTCTTGGGACCTGTGGTCCTAATGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCCTTCACCA";
  string query11 = "AGCAAAGTGGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCTTGCAGCAGAGTTGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCCTTCACCATCTAAGGAAAGGTGTCCCCCCATTCCCAATGGT";
  string query12 = "GGAGTGCCTACTTGGAGCAGCAGAGCACACTGAGTACAAGTCCCCTTGCAGCAGAGTTGCAAGAGGTCTTGGGACCTGTGGTCCTAATGCAAGATAAGGCCACGGGGCCTGAGGCACCCCTAGACCCTGAGCCTTCACCATCTAAGGAA";

  int pos1 = refPos-184258313;
  int pos2 = refPos-184258314;
  int pos3 = refPos-184258317;
  int pos4 = refPos-184258319;
  int pos5 = refPos-184258321;
  int pos6 = refPos-184258325;
  int pos7 = refPos-184258328;
  int pos8 = refPos-184258330;
  int pos9 = refPos-184258331;
  int pos10 = refPos-184258333;
  int pos11 = refPos-184258335;
  int pos12 = refPos-184258343;

  Variant v1 = {query1, pos1};
  Variant v2 = {query2, pos2};
  Variant v3 = {query3, pos3};
  Variant v4 = {query4, pos4};
  Variant v5 = {query5, pos5};
  Variant v6 = {query6, pos6};
  Variant v7 = {query7, pos7};
  Variant v8 = {query8, pos8};
  Variant v9 = {query9, pos9};
  Variant v10 = {query10, pos10};
  Variant v11 = {query11, pos11};
  Variant v12 = {query12, pos12};
  
  vector<Variant> variants;
  variants.push_back(v1);
  variants.push_back(v2);
   variants.push_back(v3);
  variants.push_back(v4);
  variants.push_back(v5);
  variants.push_back(v6);
  //TODO: figure out how to deal with position check on varaints right before SV
  // variants.push_back(v7);
  variants.push_back(v8);
  variants.push_back(v9);
  variants.push_back(v10);
  variants.push_back(v11);
  variants.push_back(v12);
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

void checkGraphEquality(vector<Traceback> tbs, int i){
  assert(tbs.size() > 0);
  if(tbs.size() == 1){
    cout << "Passed " << i << "/" << i << " graph equality checks\n";
    return;
  }
  else{
    Traceback t1 = tbs.back();
    tbs.pop_back();
    Traceback t2 = tbs.back();
    vector<Node *> sn1 = t1.subjectNodes_;
    vector<Node *> sn2 = t2.subjectNodes_;
    assert(sn1.size() == 4 && sn2.size() == 4);
    for(unsigned i = 0; i < 4; i++){
      assert(sn1[i]->getSequence() == sn2[i]->getSequence());
    }
    checkGraphEquality(tbs, ++i);
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

void checkTrimmedVariants(vector<Variant> variants){
  int s = variants.back().ref.size();
  int c = 0;
  for(auto it = std::begin(variants); it != std::end(variants); ++it){
    assert(it->ref.size() ==s);
    c++;
  }
  cout << "Passed " << c << "/" << c << " trimmed variant tests\n";
}


void runAllTests(){
  vector<Variant> variants = buildAllVariants();
  std::pair<string,string> sv = std::make_pair("TTGCAAGAGGTCTTGGGACCTGTGGTCCTAA","T");
  //std::pair<string,string> sv = std::make_pair("TGCAGCAGAGTTGCAAGAGGTCTTGGGACCT","T");
  checkPos(variants, sv);
  Variant ref = {"TGAGTACAAGTCCCTTGCAGCAGAGTTGCAAGAGGTCTTGGGACCTGTGGTCCTAATGCAAGATAAGGCCA", 27};
  //Variant ref = {"TGAGTACAAGTCCCTTGCAGCAGAGTTGCAAGAGGTCTTGGGACCTGTGGTCCTAATGCAAGATAAGGCCA", 17};

  PileUp *p = new PileUp(variants, sv, ref);
  checkGraphEquality(p->tbs_,1);
  checkTracebackVectorSize(p->sumMatrix_);
  checkTBMaxValue(p->sumMatrix_, variants);
  checkDimsMatch(p->tbs_);
  checkTrimmedVariants(p->variants_);
  printTracebacks(p->sumMatrix_);
}

int main(){
  runAllTests();
}
