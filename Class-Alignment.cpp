//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Class-Alignment
// Class code for Alignment object
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
#include <cmath>

// "boost" regular expression library
#include <boost/regex.hpp>
#include <boost/format.hpp>


#include "Class-Alignment.h"

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
using std::min;
using std::max;

using namespace std;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// utility routines
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// constants
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// static class-wide variable initializations
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// class methods
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// constructor
//------------------------------------------------------------------------------
Alignment::Alignment(
		     string subject, // subject sequence
		     string query, // query sequence
		     int M, // match score
		     int X, // mismatch score
		     int GI, // gap initiation score
		     int GE // gap extension score
		     ) {
  //----------------------------------------------------------------------------
  // derived quantities
  //----------------------------------------------------------------------------
  s1 = subject;
  s2 = query;
  l1 = s1.length();
  l2 = s2.length();

  //----------------------------------------------------------------------------
  // initialize constants
  //----------------------------------------------------------------------------
  D = 0;
  H = 1;
  V = 2;

  //----------------------------------------------------------------------------
  // initialize scores and related quantities
  //----------------------------------------------------------------------------
  I1M = 0;
  I2M = 0;
  SM = 0;

  //----------------------------------------------------------------------------
  // initialize score matrix and traceback matrices
  //----------------------------------------------------------------------------

  // initialize score matrices with all 0s
  for (int i1=0; i1<=l1; i1++) {
    vector< vector<int> > zeroVectorVector;
    for (int i2=0; i2<=l2; i2++) {
      vector<int> zeroVector;
      for (int i3=0; i3<=2; i3++) {
	zeroVector.push_back(0);
      }
      zeroVectorVector.push_back(zeroVector);
    }
    S.push_back(zeroVectorVector);
  }

  // initialize gap score matrices with negative numbers
  for (int i1=0; i1<=l1; i1++) {
    S.at(i1).at(0).at(H) = l1 * (X + GI + GE);
  }
  for (int i2=0; i2<=l2; i2++) {
    S.at(0).at(i2).at(V) = l2 * (X + GI + GE);
  }

  // initialize affine gap penalty trace-back
  for (int i1=0; i1<=l1; i1++) {
    vector< vector< vector<bool> > > falseVectorVectorVector;
    for (int i2=0; i2<=l2; i2++) {
      vector< vector<bool> > falseVectorVector;
      for (int i3=0; i3<=2; i3++) {
	vector<bool> falseVector;
	for (int i4=0; i4<=2; i4++) {
	  falseVector.push_back(false);
	}
	falseVectorVector.push_back(falseVector);
      }
      falseVectorVectorVector.push_back(falseVectorVector);
    }
    T.push_back(falseVectorVectorVector);
  }

  //----------------------------------------------------------------------------
  // fill score matrices
  //----------------------------------------------------------------------------
  for (int i1=1; i1<=l1; i1++) {
    for (int i2=1; i2<=l2; i2++) {
 
      // calculate diagonal score (match or mismatch)
      string b1 = s1.substr(i1-1, 1); // subject base
      string b2 = s2.substr(i2-1, 1); // query base	
      int mx = 0;
      if (b1 == b2) {mx = M;} else {mx = X;} // update with match or mismatch score

      // calculate and register Match optimal score
      int mdd = S.at(i1-1).at(i2-1).at(D) + mx;
      int mdh = S.at(i1-1).at(i2-1).at(H) + mx;
      int mdv = S.at(i1-1).at(i2-1).at(V) + mx;
      int md = max(max(max(mdd, mdh), mdv), 0);
      S.at(i1).at(i2).at(D) = md;

      // determine and register Match tracebacks
      if (md > 0) {
	//if (mdd == md) {T[i1).at(i2).at(D).at(D] = true;}
	if (mdd == md) {T.at(i1).at(i2).at(D).at(D) = true;}
	if (mdh == md) {T.at(i1).at(i2).at(D).at(H) = true;}
	if (mdv == md) {T.at(i1).at(i2).at(D).at(V) = true;}
      }

      // calculate and register horizontal gap optimal score
      int ihd = S.at(i1-1).at(i2).at(D) + GI;
      int ihh = S.at(i1-1).at(i2).at(H) + GE;
      int ih =  max(ihd, ihh);
      S.at(i1).at(i2).at(H) = ih;

      // determine and register Horizontal Gap tracebacks
      if (ih > 0) {
	if (ihd == ih) {T.at(i1).at(i2).at(H).at(D) = true;}
	if (ihh == ih) {T.at(i1).at(i2).at(H).at(H) = true;}
      }
    
      // calculate and register vertical gap optimal score
      int ivd = S.at(i1).at(i2-1).at(D) + GI;
      int ivv = S.at(i1).at(i2-1).at(V) + GE;
      int iv = max(ivd, ivv);
      S.at(i1).at(i2).at(V) = iv;

      // determine and register Horizontal Gap tracebacks
      if (iv > 0) {
	if (ivd == iv) {T.at(i1).at(i2).at(V).at(D) = true;}
	if (ivv == iv) {T.at(i1).at(i2).at(V).at(V) = true;}
      }
    
      // update largest affine score, and the indices or the corresponding cell
      if (md > SM) {
	SM = md;
	I1M = i1;
	I2M = i2;
      }
    }
  }
}

//------------------------------------------------------------------------------
// getScore()
//------------------------------------------------------------------------------
int Alignment::getScore() {
  return SM;
}

//------------------------------------------------------------------------------
// getSubjectLength()
//------------------------------------------------------------------------------
int Alignment::getSubjectLength() {
  return l1;
}

//------------------------------------------------------------------------------
// getQueryLength()
//------------------------------------------------------------------------------
int Alignment::getQueryLength() {
  return l2;
}

//------------------------------------------------------------------------------
// getBaseAlignment()
//------------------------------------------------------------------------------
string Alignment::getBaseAlignment() {

  // initialize data
  int i1 = I1M;
  int i2 = I2M;
  int state = D;

  int s = S.at(i1).at(i2).at(D);
  string a1 = "";
  string a2 = "";

  // construct alignment strings according to traceback
  while (s > 0 && i2 > 0 && i1 > 0) {
    // diagonal trace-back
    if (state == D) {
      string b1 = s1.substr(i1-1, 1); // subject base
      string b2 = s2.substr(i2-1, 1); // query base
      a1 = b1 + a1; // add subject base to alignment
      a2 = b2 + a2; // add quary base to alignment
      if (T.at(i1).at(i2).at(D).at(D)) {state = D; s = S.at(i1-1).at(i2-1).at(D);}
      else if (T.at(i1).at(i2).at(D).at(H)) {state = H; s = S.at(i1-1).at(i2-1).at(H);}
      else if (T.at(i1).at(i2).at(D).at(V)) {state = V; s = S.at(i1-1).at(i2-1).at(V);}
      i1--; // decrement subject index
      i2--; // decrement query index
    }
    // horizontal trace-back
    else if (state == H) {
      string b1 = s1.substr(i1-1, 1); // subject base
      string b2 = "-"; // query base is a gap
      a1 = b1 + a1; // add subject base to alignment
      a2 = b2 + a2; // add quary base to alignment
      if (T.at(i1).at(i2).at(H).at(D)) {state = D; s = S.at(i1-1).at(i2).at(D);}
      else if (T.at(i1).at(i2).at(H).at(H)) {state = H; s = S.at(i1-1).at(i2).at(H);}
      i1--; // decrement subject index only
    }
    // vertical trace-back
    else if (state == V) {
      string b1 = "-"; // subject base is a gap
      string b2 = s2.substr(i2-1, 1); // query base
      a1 = b1 + a1; // add subject base to alignment
      a2 = b2 + a2; // add quary base to alignment
      if (T.at(i1).at(i2).at(V).at(D)) {state = D; s = S.at(i1).at(i2-1).at(D);}
      else if (T.at(i1).at(i2).at(V).at(V)) {state = V; s = S.at(i1).at(i2-1).at(V);}
      i2--; // decrement query index only
    }
  }

  // return base aligment string
  return a1 + "\n" + a2;
}

//------------------------------------------------------------------------------
// printMatrix()
//------------------------------------------------------------------------------
void Alignment::printMatrix(ostream &out) {
  out << "test" << endl;
  out << "Score matrix" << endl;
  
  // first, print the table header
  out << boost::format("%1$2s | %2$27s | ") % "" % "-"; // print blank entry
  for (int i1=1; i1 <= l1; i1++) {
    out << boost::format("%1$27s | ") % s1.substr(i1-1, 1); // print subject base
  }
  out << endl; 
  out << boost::format("%1$s") % "----";
  for (int i1=0; i1 <= l1; i1++) {
    out << boost::format("%1$30s") % "------------------------------"; // print 
  }
  out << endl;
    
  // then print the table contents
  for (int i2=0; i2<=l2; i2++) {
    if (i2 == 0) {
      out << boost::format("%1$2s | ") % "-"; // print blank entry
    }
    else {
      out << boost::format("%1$2s | ") % s2.substr(i2-1, 1); // print query base
    }
    for (int i1=0; i1<=l1; i1++) {
      int sd = S.at(i1).at(i2).at(D); int tdd = T.at(i1).at(i2).at(D).at(D); int tdh = T.at(i1).at(i2).at(D).at(H); int tdv = T.at(i1).at(i2).at(D).at(V);
      int sh = S.at(i1).at(i2).at(H); int thd = T.at(i1).at(i2).at(H).at(D); int thh = T.at(i1).at(i2).at(H).at(H);
      int sv =  S.at(i1).at(i2).at(V); int tvd = T.at(i1).at(i2).at(V).at(D); int tvv = T.at(i1).at(i2).at(V).at(V);
      out << boost::format("%1$4d(%2$1d%3$1d%4$1d) %5$4d(%6$1d%7$1d) %8$4d(%9$1d%10$1d) | ") % sd % tdd % tdh % tdv % sh % thd % thh % sv % tvd % tvv;
    }
    out << endl;
    out << boost::format("%1$4s") % "----";
    for (int i1=0; i1 <= l1; i1++) {
      out << boost::format("%1$-30s") % "------------------------------"; // print 
    }
    out << endl;
  }
}


