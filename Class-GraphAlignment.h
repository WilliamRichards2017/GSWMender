//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Class-GraphAlignment
// Class definition for GraphAlignment object
// Copyright 2015 Gabor T. Marth, University of Utah
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <iterator>
#include <algorithm>

#include "Class-Node.h"

using std::ios;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::istream;
using std::fstream;
using std::cin;
using std::cout;
using std::clog;
using std::endl;
using std::string;
using std::vector;
using std::deque;
using std::map;

using namespace std;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// utility routines
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// type definitions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Match
//------------------------------------------------------------------------------
struct Match {
  string type;
  int length;
};


//------------------------------------------------------------------------------
// NodeMap
//------------------------------------------------------------------------------
struct NodeMap {
  string cigar;
  int offset;
  int score;
  int maxScore;
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// class definition
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


class GraphAlignment {

public:

  //----------------------------------------------------------------------------
  // public functions
  //----------------------------------------------------------------------------

  // constructor
  GraphAlignment(
		 vector<Node *>, // subject nodes representing Variant Graph
		 string, // query
		 int, // match score
		 int, // mismatch score
		 int, // gap initiation score
		 int, // gap extension score
		 bool // debug flag
		 );
  
  // getAssBlockOffset
  int getScore();

  // getQueryLength
  int getQueryLength();

  // printMatrix
  void printMatrix(Node *, ostream&);

  // print traceback matrix
  void printTBMatrix(Node *, ostream&);

  // getGlobalAlignment
  string getGlobalAlignment();

  // getGlobalCigar()
  string getGlobalCigar();

  // getMatchedNodes()
  vector<Node *> getMatchedNodes();

  // getNodeCigar
  string getNodeCigar(Node *);

  // getNodeoffset
  int getNodeOffset(Node *);

  int *getTraceBackPointer(Node *);

  map<Node *, vector< vector< vector<int> > >, less<Node *> > getScoreMatrix(); // Score matrices                                                                                                                                                           
  map<Node *, vector< vector< vector< vector<bool> > > >, less<Node *> > getTBMatrix(); // traceback matrices  

  map<Node *, vector< vector< vector<int> > >, less<Node *> > GS; // Score matrices                                                                                                                                                        
  map<Node *, vector< vector< vector< vector<bool> > > >, less<Node *> > GT; // traceback matrices  
  

  //----------------------------------------------------------------------------
  // public variables
  //----------------------------------------------------------------------------
  
 private:

  //----------------------------------------------------------------------------
  // private variables
  //----------------------------------------------------------------------------

  // trace-back direction indices
  int D;
  int H;
  int V;

  // subject and query related variables
  string s2; // query sequence
  int l2; // query length


  // optimum score related variables
  int SM; // best alignment score
  Node * NodeM; // subject node with the best alignment score
  int I1M; // subject index in score matrix corresponding to best alignment score
  int I2M; // query index in score matrix corresponding to best alignment score
  int StateM; // Match/horizontal gap/vertical gap state corresponding to best alignment score
  string globalCigar;
  string globalAlignment;
  vector<Node *> matchedNodes;
  map<Node *, string, less<Node *> > nodeCigarMap;
  map<Node *, int, less<Node *> > nodeOffsetMap;
};
