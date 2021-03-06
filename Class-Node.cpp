//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Class-Node
// Class code for Varioant Graph Node object
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
Node::Node(
	   string idIn, // ID
	   string sequenceIn, // node sequence
	   vector<Node *> contributorsIn, // contributor nodes
	   int nodeType // 0 = same, 1 = ref, 2 = alt
	   ) {
  //----------------------------------------------------------------------------
  // put input into object fields
  //----------------------------------------------------------------------------
  id = idIn;
  sequence = sequenceIn;
  contributors = contributorsIn;
  type = nodeType;
}

//------------------------------------------------------------------------------
// getId()
//------------------------------------------------------------------------------
string Node::getId() {
  return id;
}

bool Node::isSame() {
  return type == 0;
}

bool Node::isRef() {
  return type == 1;
}

bool Node::isAlt() {
  return type == 2;
}

//------------------------------------------------------------------------------
// getSubjectLength()
//------------------------------------------------------------------------------
string Node::getSequence() {
  return sequence;
}

void Node::setSequence(string seq) {
  sequence = seq;
}

void Node::pushLast(Node *node1, Node *node2){
  string c = string(1, node1->getSequence().back());
  node2->setSequence(c.append(node2->getSequence()));
  node1->setSequence(node1->getSequence().substr(0,node1->getSequence().length()-1));
}

void Node::pullFirst(Node *node1, Node *node2){
  string c = string(1, node2->getSequence().at(0));
  node1->setSequence(node1->getSequence().append(c));
  node2->setSequence(node2->getSequence().substr(1, node2->getSequence().length()));
}

void Node::undoPush(Node *node1, Node *node2){
  pullFirst(node1, node2);
}

void Node::undoPull(Node *node1, Node *node2){
  pushLast(node1, node2);
}

//------------------------------------------------------------------------------
// getContributors()
//------------------------------------------------------------------------------
vector<Node *> Node::getContributorNodes() {
  return contributors;
}

//------------------------------------------------------------------------------
// printNode()
//------------------------------------------------------------------------------
void Node::printNode(ostream &out) {
  out << "Node ID=" << id << endl;
  out << "Node sequence=" << sequence << endl;
  out << "Contributor nodes=(";
  for (vector< Node *>::const_iterator iter = contributors.begin(); iter != contributors.end(); iter++) {
    Node * node = * iter;
    out << node->getId() << ",";
  }
  out << ")" << endl;
}


