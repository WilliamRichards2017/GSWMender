#ifndef __PILEUP_H_INCLUDED__
#define __PILEUP_H_INCLUDED__

#include "Traceback.h"
#include "Variant.h"

class Traceback;
class Node;

class PileUp{

 public:
  PileUp(vector<Variant>,std::pair<string,string>,Variant);
  ~PileUp();
  Variant ref_;
  std::pair<string,string> sv_;
  vector<Traceback> tbs_;
  vector<Variant> variants_;
  vector<Node *> subjectNodes_;
  vector<vector<vector<int> > > sumMatrix_;
  
 private:
  vector<string> getNodes();
  vector<Node *> buildDiamondGraph(vector<string>);
  vector<vector<vector<int> > > sumTracebacks();
  void trimVariants();
  void buildTracebackVector(vector<Variant>);
  void deleteGraph();
};



#endif // PILEUP_H_INCLUDED
