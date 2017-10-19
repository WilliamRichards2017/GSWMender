//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Gsw
// Copyright 2015 Gabor T. Marth, University of Utah
// Modified by Will Richards 2017
// All rights reserved.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// standard includes
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>

// "tclap" commandline parsing library
#include <tclap/CmdLine.h>

// "boost" regular expression library
#include <boost/regex.hpp>
#include <boost/format.hpp>

// include classes
#include "Class-GraphAlignment.h"
#include "Traceback.h"
#include "PileUp.h"

// uses
using namespace std;
using namespace TCLAP;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Begining of GSWMender helper functions and container structs
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

struct Graph{
  GraphAlignment *alignment;
  vector<Node *> nodes;
};


struct ArgStruct {
  string shortId;
  string longId;
  string description;
  bool required;
  string defaultValueString;
  int defaultValueInt;
  double defaultValueDouble;
  bool defaultValueBool;
  string type;
  bool multi;
  vector<string> constraint;
};

struct AlleleData {
  string base;
  short qual;
};

//Generate the nodes of a diamond graph from a Variant (query, variant, position)
vector<string > getNodes(Variant v){
  vector<string> strings;
  string s1 = v.ref.substr(0, v.pos);
  string s2 = v.sv.first;
  string s3 = v.sv.second;
  string s4 = v.ref.substr(v.pos + v.sv.first.length(), v.ref.length());
  strings.push_back(s1);
  strings.push_back(s2);
  strings.push_back(s3);
  strings.push_back(s4);

  cout << " ref path is: " << s1 + s2 + s4 << std::endl;
  cout << "alt path is:  " << s1 + s3 + s4 << std::endl << std::endl;
  return strings;
}



vector<Node *> buildDiamondGraph(vector<string> strings){

  vector<Node *> subjectNodes;
  vector<Node * > contributors1;
  Node * node1 = new Node(
			  "node1",
			  strings[0],
			  contributors1,
			    0
			  );
  subjectNodes.push_back(node1);
  vector<Node *> contributors2;
  contributors2.push_back(node1);
  Node * node2 = new Node(
			  "node2",
			  strings[1],
			  contributors2,
			    1
			  );
  subjectNodes.push_back(node2);
  vector<Node *> contributors3;
  contributors3.push_back(node1);
  Node * node3 = new Node(
			  "node3",
			  strings[2],
			  contributors3,
			    2
			  );
  subjectNodes.push_back(node3);
  vector<Node *> contributors4;
  contributors4.push_back(node2);
  contributors4.push_back(node3);
  Node * node4 =  new Node(
			   "node4",
			   strings[3],
			   contributors4,
			     0
			   );
  subjectNodes.push_back(node4);
  return subjectNodes;
}

void deleteGraph(vector<Node *> nodes){
  for(auto it = std::begin(nodes); it != std::end(nodes); ++it){
    delete * it;
  }
}

GraphAlignment *updateGA(GraphAlignment *ga, vector<Node *> subjectNodes, string query, int M, int X, int GI, int GE, bool debug){
  delete ga;
  ga = new GraphAlignment(subjectNodes, query, M, X, GI, GE, debug);
  return ga;
}


int getLongestString(vector<string> strings){
  int i = -1;
  for(auto it = std::begin(strings); it != std::end(strings); ++it){
    if(it->length() > i){
      i = it->length();
    }
  }
  return i;
}

std::pair<int,int> getLongestVariants(vector<std::pair<int,int> >){

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Begining of GSWMender core functions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
Graph pullFromSame(vector<Node *> subjectNodes, GraphAlignment *ga, string query, int M, int X, int GI, int GE, bool debug){
  for(vector<Node *>::const_iterator iter = subjectNodes.begin(); iter != subjectNodes.end(); iter++){
    Node * node = * iter;
    vector<Node *> contributorNodes = node->getContributorNodes();
    int prevScore = ga->getScore();
    int nextScore = prevScore;

    if(node->getContributorNodes().size() > 1){
      while(nextScore >= prevScore){
        Node *ref = node->getContributorNodes()[0];
        if(node->getSequence().size() > 0){
          ref->pullFirst(ref, node);
        }
        else {
          break;
        }
        prevScore = ga->getScore();
        ga = updateGA(ga, subjectNodes, query, M, X, GI, GE, debug);
        nextScore = ga->getScore();
        if(prevScore > nextScore){
          ref->undoPull(ref, node);
          ga = updateGA(ga, subjectNodes, query, M, X, GI, GE, debug);
          break;
        }
      }
    }
  }
  Graph g = {ga, subjectNodes};
  return g;
}

Graph pushToRef(vector<Node * > subjectNodes, GraphAlignment *ga, string query, int M, int X, int GI, int GE, bool debug){
  for(vector<Node *>::const_iterator iter = subjectNodes.begin(); iter != subjectNodes.end(); iter++){
    Node * node = * iter;
    vector<Node *> contributorNodes = node->getContributorNodes();
    ga = updateGA(ga, subjectNodes, query, M, X, GI, GE, debug);
    int prevScore = ga->getScore();
    int nextScore = prevScore;
    if(node->isRef()){
      while(nextScore >= prevScore){
        if (contributorNodes[0]->getSequence().length() > 0){
          node->pushLast(contributorNodes[0], node);
        }
        else {break;}
	prevScore = ga->getScore();
        ga = updateGA(ga, subjectNodes, query, M, X, GI, GE, debug);
        nextScore = ga->getScore();
        if(prevScore > nextScore){
          node->undoPush(contributorNodes[0], node);
          ga = updateGA(ga, subjectNodes, query, M, X, GI, GE, debug);
          break;
        }
      }
    }
  }
  Graph g = {ga, subjectNodes};
  return g;
}

Graph refit(vector<Node *> subjectNodes, GraphAlignment *ga, string query, int M, int X, int GI, int GE, bool debug){
  ga = updateGA(ga, subjectNodes, query, M, X, GI, GE, debug);
  cout << "optimal alignment score before refit: " << ga->getScore() << std::endl;
  cout << "Global Alignment:" << endl << ga->getGlobalAlignment() << endl;
  Graph g = pushToRef(subjectNodes, ga, query, M, X, GI, GE, debug);
  ga = g.alignment;
  subjectNodes = g.nodes;
  g = pullFromSame(subjectNodes, ga, query, M, X, GI, GE, debug);
  return g;
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Util  helper functions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

vector<Traceback> buildTracebackVector(vector<Variant> variants){
  vector<Traceback> tracebackVec;

  for(auto it = std::begin(variants); it != std::end(variants); ++it){
    vector<Node *> subjectNodes = buildDiamondGraph(getNodes(*it));
    GraphAlignment * ga = new GraphAlignment(subjectNodes, it->ref, 2, -2, -3, -2, &pause);
    Traceback t = {subjectNodes, ga};
    tracebackVec.push_back(t);
  }
  return tracebackVec;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// static variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

static string ProgramName("gsw");
static string ProgramDescription("Smith-Waterman optimal alignment algorithm");
static string ProgramVersion("0.4.1");
static string ProgramDate("2015-04-18");
static string ProgramDeveloper("Gabor T. Marth");
static string ProgramInstitution("University of Utah");
static string ProgramCopyrightDates("2015");

static vector<ArgStruct> ArgList;

class MyOutput : public StdOutput {
public:

  virtual void failure(CmdLineInterface& c, ArgException& e)
  {
    cerr << "################################################################################" << endl;
    cerr << "### Program: " << ProgramName << " Version: " <<  ProgramVersion << " Release date: " << ProgramDate << endl;
    cerr << "### " << ProgramDescription << endl;
    cerr << "### " << "Copyright: " << ProgramCopyrightDates << " " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cerr << "### All rights reserved." << endl;
    cerr << "###" << endl;
    cerr << "### Command line error:" << e.what() << endl;
    cerr << "### For usage please type: " << c.getProgramName() << " --help" << endl;
    cerr << "################################################################################" << endl;
  }

  virtual void usage(CmdLineInterface& c)
  {

    cout << "################################################################################" << endl;
    cout << "### Program: " << ProgramName << " Version: " <<  ProgramVersion << " Release date: " << ProgramDate << endl;
    cout << "### " << ProgramDescription << endl;
    cout << "### " << "Copyright: " << ProgramCopyrightDates << " " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cout << "### All rights reserved." << endl;
    cout << "###" << endl;
    cout << "### Usage: " << c.getProgramName() << " [arguments], where:" << endl;
    for(vector<ArgStruct>::const_iterator it = ArgList.begin();
	it != ArgList.end(); it++) {
      ArgStruct arg = *it;

      string idString = "";
      if (arg.longId != "") {idString += "--" + arg.longId;}
      if (arg.shortId != "") {idString += ", -" + arg.shortId;}

      string multiString = "single-valued";
      if (arg.multi) {multiString = "multi-valued";}

      if (arg.required) {
	cout << "### " << idString << " [" << arg.type << ", required, no default, " << multiString << "]" << endl;
      }
      else {
	cout << "### " << idString << " [" << arg.type << ", optional, default=" << arg.defaultValueString << ", " << multiString << "]" << endl;
      }
      if (arg.constraint.size() > 0) {
	cout << "###     Permitted values: (";
	bool first = true;
	for (vector<string>::const_iterator iter = arg.constraint.begin();
	     iter != arg.constraint.end(); iter++) {
	  string value = *iter;
	  if (! first) {
	    cout << "|";
	  }
	  first = false;
	  cout << value;
	}
	cout << ")" << endl;
      }
      cout << "###     Description: " << arg.description << endl;
    }
    cout << "################################################################################" << endl;
  }

  virtual void version(CmdLineInterface& c)
  {
    cerr << "################################################################################" << endl;
    cerr << "### Program: " << ProgramName << " Version: " <<  ProgramVersion << " Release date: " << ProgramDate << endl;
    cout << "### " << ProgramDescription << endl;
    cout << "### " << "Copyright: " << ProgramCopyrightDates << " " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cout << "### All rights reserved." << endl;
    cout << "###" << endl;
    cerr << "################################################################################" << endl;
  }
}; // end of myOutput class


int main (int argc, char *argv[]) {

  // constants
  static int D = 0;
  static int H = 1;
  static int V = 2;

  // Create new CmdLine object
  CmdLine cmd("", ' ', ProgramVersion);

  // add program-specific command line arguments

  // initialize arg
  ArgStruct arg;

  // subject
  ArgStruct argSubject;
  arg = argSubject;
  arg.shortId = "s";
  arg.longId = "subject";
  arg.description = "Subject sequence";
  arg.required = false;
  //arg.defaultValueString = "ACGT";
  arg.defaultValueString = "CCCCACCCACCCC";
  arg.type = "string";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_subject(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // query
  ArgStruct argQuery;
  arg = argQuery;
  arg.shortId = "q";
  arg.longId = "query";
  arg.description = "Query sequence";
  arg.required = false;
  //arg.defaultValueString = "CTATTTTAGTAGGTTGTTA";
  arg.defaultValueString = "TAAAGCCGATTGTTTTGTGCT";
  //arg.defaultValueString = "ACGT";
  arg.type = "string";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_query(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // match score
  ArgStruct argM;
  arg = argM;
  arg.shortId = "";
  arg.longId = "M";
  arg.description = "Match score";
  arg.required = false;
  arg.defaultValueString = "2";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_M(arg.shortId, arg.longId, arg.description, arg.required, 2, arg.type, cmd);

  // mismatch score
  ArgStruct argX;
  arg = argX;
  arg.shortId = "";
  arg.longId = "X";
  arg.description = "Mismatch score";
  arg.required = false;
  arg.defaultValueString = "-2";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_X(arg.shortId, arg.longId, arg.description, arg.required, -2, arg.type, cmd);

  // gap initiation score
  ArgStruct argGI;
  arg = argX;
  arg.shortId = "";
  arg.longId = "GI";
  arg.description = "Gap initiation score";
  arg.required = false;
  arg.defaultValueString = "-3";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_GI(arg.shortId, arg.longId, arg.description, arg.required, -3, arg.type, cmd);

  // gap extension score
  ArgStruct argGE;
  arg = argX;
  arg.shortId = "";
  arg.longId = "GE";
  arg.description = "Gap extension score";
  arg.required = false;
  arg.defaultValueString = "-1";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_GE(arg.shortId, arg.longId, arg.description, arg.required, -1, arg.type, cmd);

  // print matrix?
  ArgStruct argMatrix;
  arg = argMatrix;
  arg.shortId = "";
  arg.longId = "matrix";
  arg.description = "Print score matrix?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_matrix(arg.shortId, arg.longId, arg.description, cmd, false);

  // debug
  ArgStruct argDebug;
  arg = argDebug;
  arg.shortId = "";
  arg.longId = "debug";
  arg.description = "Print debugging messages?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_debug(arg.shortId, arg.longId, arg.description, cmd, false);

  //----------------------------------------------------------------------------
  // register (but not add to cmd) special arguments that are automatically
  // added to cmd
  //----------------------------------------------------------------------------

  // help
  ArgStruct argHelp;
  arg = argHelp;
  arg.shortId = "h";
  arg.longId = "help";
  arg.description = "Print usage statement?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);

  // version
  ArgStruct argVersion;
  arg = argVersion;
  arg.shortId = "";
  arg.longId = "version";
  arg.description = "Print program version?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);

  // add custom output handler
  MyOutput my;
  cmd.setOutput(&my);

  // parse command line and catch possible errors
  try {
    cmd.parse(argc,argv);
  }
  catch ( ArgException& e ) {
    cerr << "ERROR: " << e.error() << " " << e.argId() << endl;
  }

  // assign command line parameters

  string subject = cmd_subject.getValue();
  string query = cmd_query.getValue();
  int M = cmd_M.getValue();
  int X = cmd_X.getValue();
  int GI = cmd_GI.getValue();
  int GE = cmd_GE.getValue();
  bool matrix = cmd_matrix.getValue();
  bool debug = cmd_debug.getValue();

  // report command line and parameters used
  map<bool, string, less<bool> > bool2String;
  bool2String[false] = "false";
  bool2String[true] = "true";

  if (debug) {
    cerr << "Command line:";
    for (int i=0; i<argc; i++) {
      cerr << " " << argv[i];
    }
    cerr << endl;
    cerr << endl;
    cerr << "Complete list of parameter values:" << endl;
    cerr << "  --subject = " << subject << endl;
    cerr << "  --query = " << query << endl;
    cerr << "  --M = " << M << endl;
    cerr << "  --X = " << X << endl;
    cerr << "  --GI = " << GI << endl;
    cerr << "  --GE = " << GE << endl;
    cerr << "  --matrix = " <<  bool2String[matrix] << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // main code starts here
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  /*string query1 = "ATCCAGTAATCCGGGATCCAT";
  string query2 = "ATCCAGTATCCGGGATCCAT";
  string query3 = "TATCCAGTATCCGGGTT"; */
  
  /* string query1 = "CCCCCCCCCCCC";
  string query2 = "CCCCCCCCCCCC";
  string query3 = "CCCCCCCCCCCC";
  string query4 = "CCCCACCCACCC";
  string query5 = "CCCCACCCACCC"; 
  */

  string query1 = "TAAAGCTGTGTTGTGCT";
  string query2 = "TAAAGCTGTTGTGCT";
  string query3 = "TAAAGCGTGTGTGCT";
  string query4 = "TAAAGCGTGTTGTGCT";
  //string query5 = "TAAAGCTGTGTGTGTGCT";

  vector<string> queries;
  queries.push_back(query1);
  queries.push_back(query2);
  queries.push_back(query3);
  queries.push_back(query4);
  //queries.push_back(query5);

  //cout << "\nLongest string is: " << getLongestString(queries) << std::endl;

  /*pair<string, string> sv1 = std::make_pair("AATCC", "A");
  pair<string, string> sv2 = std::make_pair("ATCC", "A");
  pair<string, string> sv3 = std::make_pair("A", "A");*/

  /*std::pair<string,string> sv1 = std::make_pair("CCCCCC","CCCCCC");
  std::pair<string,string> sv2 = std::make_pair("CCCCCC","CCCCCC");
  std::pair<string,string> sv3 = std::make_pair("CCCCCC","CCCCCC");
  std::pair<string,string> sv4 = std::make_pair("CACCCA","CCCCCC");*/

  std::pair<string,string> sv1 = std::make_pair("CGATTGTTT","TGTGT");
  std::pair<string,string> sv2 = std::make_pair("CGATTGTTT", "TGT");
  std::pair<string,string> sv3 = std::make_pair("CGATTGTTT", "GTG");
  std::pair<string,string> sv4 = std::make_pair("CGATTGTTT", "GTGT");

  int pos = 6;

  vector<int> positions;
  positions.push_back(pos);
  positions.push_back(pos);
  positions.push_back(pos);
  positions.push_back(pos);

  vector<Variant> variants;
  Variant v1 = {query1, sv1, pos};
  Variant v2 = {query2, sv2, pos};
  Variant v3 = {query3, sv3, pos};
  Variant v4 = {query4, sv4, pos};

  variants.push_back(v1);
  variants.push_back(v2);
  variants.push_back(v3);
  variants.push_back(v4);

  //uncomment to debug node values
  /*vector<string> s = getNodes(v1);
  for(auto it = std::begin(s); it != std::end(s); ++it){
    cout << *it << std::endl;
    }*/

  //Create a graph alignment 
  vector<Node *> subjectNodes = buildDiamondGraph(getNodes(v1));
  GraphAlignment *ga = new GraphAlignment(subjectNodes, query1, M, X, GI, GE, debug);

  vector<Traceback> tracebacks = buildTracebackVector(variants);

  PileUp p = {tracebacks};
  vector<vector<vector<int> > > pileup = p.sumTracebacks();


  
  /*Graph g = refit(subjectNodes, ga, query, M, X, GI, GE, debug);
  ga = g.alignment;
  subjectNodes = g.nodes;
  
  Traceback TBAfter = {subjectNodes, ga};
  vector<vector<vector<int> > > after = TBAfter.buildTB();
  */

 

  cout << "Optimal score of GSW for first variant: " << ga->getScore() << endl;
  cout << "Global Cigar:" << ga->getGlobalCigar() << endl;
  cout << "Global Alignment:" << endl << ga->getGlobalAlignment() << endl;

  vector<Node *> matchedNodes = ga->getMatchedNodes();
  cout << "Graph node alignments:" << endl;
  for (vector<Node *>::const_iterator iter = subjectNodes.begin(); iter != subjectNodes.end(); iter++) {
    Node * node = * iter;
    //string cigar = ga->getNodeCigar(node);
    //int offset = ga->getNodeOffset(node);
    //cout << "  Node=" << node->getId() << " CIGAR=" << cigar << " offset=" << offset << endl;
    ga->printMatrix(node, cout);
  }
}

