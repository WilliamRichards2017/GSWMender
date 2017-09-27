//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Class-Node
// Class definition for Variang Graph Node object
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
//------------------------------------------------------------------------------
// class definition
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


class Node {

public:

  //----------------------------------------------------------------------------
  // public functions
  //----------------------------------------------------------------------------

  // constructor
  Node(
       string, // id
       string, // sequence
       vector<Node *>, // list of contributor nodes
       int // node type, 0 = same, 1 = ref, 2 = alt  
       );
  
  // getId
  string getId();

  // getSequence
  string getSequence();

  bool isSame();

  bool isAlt();

  bool isRef();

  void pushLast(Node *, Node *);
  
  void pullFirst(Node *, Node *);

  void undoPush(Node *, Node *);

  void undoPull(Node *, Node *);

  void setSequence(string);

  // getContributorNodes
  vector<Node *> getContributorNodes();

  // printNode
  void printNode(ostream&);

  //----------------------------------------------------------------------------
  // public variables
  //----------------------------------------------------------------------------

private:

  //----------------------------------------------------------------------------
  // private variables
  //----------------------------------------------------------------------------
  string id;
  string sequence;
  vector<Node *> contributors;
  int type;
};

