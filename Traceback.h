#ifndef __TRACEBACK_H_INCLUDED__
#define __TRACEBACK_H_INCLUDED__


#include "Class-GraphAlignment.h"

class Traceback {

 public:
  vector<Node *> subjectNodes_;
  GraphAlignment * ga_;
  vector<vector<int> > MVM_;
  vector<vector<vector<int> > > TBMs_;

 private:
  pair<int,int> getMaxCoords();
  void buildTBMs();

};

#endif // __TRACEBACK_H_INCLUDED__
