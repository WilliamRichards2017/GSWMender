#ifndef __PILEUP_H_INCLUDED__
#define __PILEUP_H_INCLUDED__

#include gsw.cpp

//forward declared depedancies
class Traceback;
class Variant;

class PileUp{

 public:
  vector<Traceback> tbs;
  vector<Variant> variants;
 private:
  vector<vector<string> > getAllNodes();
  vector<vector<Node *> > buildAllGraphs(vector<vector<string> >);
  vector<vector<vector<int> > > sumTracebacks();
}



#endif __PILEUP_H_INCLUDED__
