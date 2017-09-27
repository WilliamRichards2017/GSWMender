//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Class-Alignment
// Class definition for Alignment object
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


class Alignment {

public:

  //----------------------------------------------------------------------------
  // public functions
  //----------------------------------------------------------------------------

  // constructor
  Alignment(
	    string, // subject
	    string, // query
	    int, // match score
	    int, // mismatch score
	    int, // gap initiation score
	    int // gap extension score
	   );
  
  // getAssBlockOffset
  int getScore();

  // getSubjectLength
  int getSubjectLength();

  // getQueryLength
  int getQueryLength();

  // printMatrix
  void printMatrix(ostream&);

  // getBaseAlignment
  string getBaseAlignment();

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
  string s1; // subject sequence
  string s2; // query sequence
  int l1; // subject length
  int l2; // query length

  // score and traceback matrices
  vector< vector< vector<int> > > S; // Score matrices
  vector< vector< vector< vector<bool> > > > T; // traceback matrices

  // optimum score related variables
  int SM; // best alignment score
  int I1M; // subject index in score matrix corresponding to best alignment score
  int I2M; // query index in score matrix corresponding to best alignment score
};

