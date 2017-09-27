#include "Traceback.h"
#include "Class-Node.h"

Traceback::Traceback(vector<Node*> subjectNodes, GraphAlignment *ga){

  vector< vector < vector< int> > > scoreMatrix = ga->GS;

  for (vector<Node *>::const_iterator iter = contributorNodes.begin(); iter != contributorNodes.end(); iter++){
    Node * contributorNode = * iter;
  }   
  
}
