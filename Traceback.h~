#ifndef __TRACEBACK_H_INCLUDED__
#define __TRACEBACK_H_INCLUDED__


#include "Class-GraphAlignment.h"

class Traceback {

 public:
  vector<Node *> subjectNodes_;
  GraphAlignment * ga_;
  vector<vector<int> > MVM_;
  vector<vector<vector<int> > > TBMs_;
  vector<vector<vector<int> > > buildTBMs();

 private:
  pair<int,int> getMaxCoords();


};

#endif // __TRACEBACK_H_INCLUDED__
