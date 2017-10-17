#ifndef __PILEUP_H_INCLUDED__
#define __PILEUP_H_INCLUDED__

#include gsw.cpp

class Traceback;
class Variant;
class Node;

class PileUp{

 public:
  vector<Traceback *> tbs_;
  vector<Variant *> variants_;
  vector<Node *> subjectNodes_;
 private:
  vector<string> getNodes(Variant v);
  vector<vector<string> > getAllNodes();
  vector<Node *>buildDiamondGraph(vector);
  vector<vector<Node *> > buildAllGraphs(vector<vector<string> >);
  vector<vector<vector<int> > > sumTracebacks();
  void deleteGraph(vector<Node *>);
  vector<vector<vector<int> > > sumTracebacks();
}



#endif __PILEUP_H_INCLUDED__
