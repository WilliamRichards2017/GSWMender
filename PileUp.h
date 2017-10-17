#ifndef __PILEUP_H_INCLUDED__
#define __PILEUP_H_INCLUDED__

#include "Traceback.h"
#include "Variant.h"

class Traceback;
class Node;

class PileUp{

 public:
  vector<Traceback> tbs_;
  vector<Variant> variants_;
  vector<Node *> subjectNodes_;
  vector<vector<vector<int> > > sumTracebacks();
 private:
  vector<string> getNodes(Variant v);
  vector<vector<string> > getAllNodes();
  vector<Node *>buildDiamondGraph(vector<string>);
  vector<vector<Node *> > buildAllGraphs(vector<vector<string> >);
  void deleteGraph();
};



#endif // PILEUP_H_INCLUDED
