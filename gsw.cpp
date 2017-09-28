//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// gsw
// Copyright 2015 Gabor T. Marth, University of Utah
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

// uses
using namespace std; 
using namespace TCLAP; 

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// templates
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// subroutines
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// typedefs
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

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

struct Graph{
  GraphAlignment *alignment;
  vector<Node *> nodes;
};

struct Variant{
  string ref;
  std::pair<string, string> sv;
  int pos;

};

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
  cout << "alt path is:  " << s1 + s3 + s4 << std::endl;
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

struct Traceback {

  const vector<Node *> _subjectNodes;
  GraphAlignment * ga;

  int** buildArray2D(unsigned height, unsigned width){
    int** array = 0;
    array = new int*[height];
    
    for(int h = 0; h < height; h++){
      array[h] = new int[width];

      for (int w = 0; w < width; w++){
	array[h][w] = 0;
      }
    }
    return array;
   
  }


  pair<int, int>  getMaxCoords(int **MVM, int h, int w){
    int maxV = -1;
    int x = -1;
    int y = -1;

    for(int i = h; i > 0; i--){
      for(int j = w; j > 0; j--){
	if(MVM[i][j] > maxV){
          y = i;
          x = j;
          maxV = MVM[i][j];
        }
      }
    }
    pair<int, int> coords = std::make_pair(x,y);
    return coords;
    
  }

  
  vector<std::pair<int, int> > buildMatrixSizeVector() {
    int h = ga->getQueryLength()+1;
    vector<std::pair<int,int> > dimsVec;
  
    for(auto it = std::begin(_subjectNodes); it != std::end(_subjectNodes); ++it){
      Node * node = * it;
      //int w = node->getSequence().length()+1;
      //std::pair<int, int> dims = std::make_pair(h,w);
      //dimsVec.push_back(dims);

      dimsVec.push_back(std::make_pair(h,node->getSequence().length()+1));
    }
    return dimsVec;
  }
  
  vector<int **> buildTB(){

    map<Node *, vector< vector< vector<int> > >, less<Node *> > GS = ga->getScoreMatrix();


    int l2 = ga->getQueryLength();
    vector<int** > TBMs;

    for (vector<Node *>::const_iterator iter = _subjectNodes.begin(); iter != _subjectNodes.end(); iter++) {
      Node * node = * iter;

      int l1 = node->getSequence().length();
      int** MVM = buildArray2D(l2+1, l1+1);
      int** TBM = buildArray2D(l2+1, l1+1);


      //initialize Maximum value matrix

      vector< vector< vector<int> > > S = GS[node];

      for (int i1=1; i1<=l1; i1++) {
        for (int i2=1; i2<=l2; i2++) {
          MVM[i2][i1] = max(max(S[i1][i2][1],S[i1][i2][2]),S[i1][i2][0]);
        }
      }

      std::pair<int, int> coords = getMaxCoords(MVM, l2, l1);
      int x = coords.first;
      int y = coords.second;
      
      cout << "coords are: " << x << ", " << y << std::endl; 
      //start at max value coords
      while(x > 0 && y > 0){
	TBM[y][x] = 1;
	
	//Move up 
	if(MVM[y-1][x] > MVM[y-1][x-1] && MVM[y-1][x] > MVM[y][x-1]){
	  y--;
	}
	//Move diagonal
	else if(MVM[y-1][x-1] > MVM[y-1][x] && MVM[y-1][x-1] > MVM[y][x-1]){
	  y--;
	  x--;
	}
	//move left
	else{
	  x--;
	}
      } // end of while

      for (int i = 0; i < l2+1; i++){
	for(int j = 0; j < l1+1; j++){
	  std::cout << TBM[i][j] << ' ';
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
      TBMs.push_back(TBM);
    }// end of node loop
    return TBMs;
  }// end of buildTB()
};


struct PileUp{
  //traceback matrices
  vector<Traceback> tbs;
  vector<Variant> variants;

  vector<vector<string> > getAllNodes(){

    vector<vector<string> > allNodes;
    for(auto it = std::begin(variants); it != std::end(variants); ++it){
      Variant v = *it;
      vector<string> nodes = getNodes(v);
      allNodes.push_back(nodes);
    }
    return allNodes;
  }

  
  vector<vector<Node *> >buildAllGraphs(vector<vector<string> > allStrings){

    vector<vector<Node *> > allGraphs;
    
    for(auto it = std::begin(allStrings); it != std::end(allStrings); ++it){
      vector<string> strings = *it;
      vector<Node *> graph = buildDiamondGraph(strings);
      allGraphs.push_back(graph);
    }
    return allGraphs;
  }

  vector <int**> sumTracebacks() {
    
    //NOTE: Must free this memory by hand, DONT FORGET!!!!!
    //int** sumMatrix = (int **)calloc(w*h, sizeof(int*));

    vector<int**>  sumMatrix;
    
    vector<vector<string> > strings = getAllNodes();
    vector<vector<Node *> > nodes = buildAllGraphs(strings);

    int count = 0;
    for(auto it = std::begin(tbs); it != std::end(tbs); ++it){
      Traceback tb = *it;
      vector<int**> matrices = tb.buildTB();
      vector<std::pair<int, int> > dims = tb.buildMatrixSizeVector();
      vector<Node *> subjectNodes = tb._subjectNodes;
      
      
      
    }

    /*int sumMatrix[h][w];
    memset(sumMatrix, 0, sizeof sumMatrix);

    for(auto it = std::begin(tbs); it != std::end(tbs); ++it){
      int** tb = *it;
      for(unsigned i = 0; i < h; i++){
	for(unsigned j = 0; j < w; j++) {
	  sumMatrix[i][j] += tb[i][j];
	}
      }
    } // end of auto

    std::cout << "---------SumMatrix----------\n";
    for (int i = 0; i < h; i++){
      for(int j = 0; j < w; j++){
	std::cout << sumMatrix[i][j] << ' ';
      }
      std::cout << std::endl;
    }
    std::cout << "---------SumMatrix----------\n";
    */

    return sumMatrix;
  }
};


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
};

void print_any_array(int array[], size_t elements) {
  for(size_t i = 0; i < elements; i ++) {
    std::cout << array[i] << ' ';
  }
}

GraphAlignment *updateGA(vector<Node *> subjectNodes, string query, int M, int X, int GI, int GE, bool debug){
  GraphAlignment * ga;
  ga = new GraphAlignment(subjectNodes, query, M, X, GI, GE, debug);
  return ga;
}

Graph pushToRef(vector<Node * > subjectNodes, GraphAlignment *ga, string query, int M, int X, int GI, int GE, bool debug){
  for(vector<Node *>::const_iterator iter = subjectNodes.begin(); iter != subjectNodes.end(); iter++){
    Node * node = * iter;
    vector<Node *> contributorNodes = node->getContributorNodes();
    ga = updateGA(subjectNodes, query, M, X, GI, GE, debug);
    int prevScore = ga->getScore();
    int nextScore = prevScore;
    
    if(node->isRef()){
      while(nextScore >= prevScore){
	if (contributorNodes[0]->getSequence().length() > 0){
	  node->pushLast(contributorNodes[0], node);
       	}
	else{
	  break;
	}
	prevScore = ga->getScore();
	ga = updateGA(subjectNodes, query, M, X, GI, GE, debug);
	nextScore = ga->getScore();

	if(prevScore > nextScore){
	  node->undoPush(contributorNodes[0], node);
	  ga = updateGA(subjectNodes, query, M, X, GI, GE, debug);
	  break;
	}
      }
    }
  }
  Graph g = {ga, subjectNodes};
  return g;
}

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
	ga = updateGA(subjectNodes, query, M, X, GI, GE, debug);
	nextScore = ga->getScore();
	if(prevScore > nextScore){
	  ref->undoPull(ref, node);
	  ga = updateGA(subjectNodes, query, M, X, GI, GE, debug);
	  break;
	}
      }
    }
  }
  Graph g = {ga, subjectNodes};
  return g;
}

vector<Graph * > buildAllGraphs(vector<Variant *> variants) {
  for(vector<Variant *>::const_iterator iter = variants.begin(); iter!= variants.end(); iter++){
    Variant * variant = * iter;
  }

} 

Graph refit(vector<Node *> subjectNodes, GraphAlignment *ga, string query, int M, int X, int GI, int GE, bool debug){
  cout << "optimal alignment score before refit: " << ga->getScore() << std::endl;

  cout << "Global Alignment:" << endl << ga->getGlobalAlignment() << endl;

  Graph g = pushToRef(subjectNodes, ga, query, M, X, GI, GE, debug);
  ga = g.alignment;
  subjectNodes = g.nodes;
  g = pullFromSame(subjectNodes, ga, query, M, X, GI, GE, debug);
  return g;
}



vector<vector< int> > *build_traceback_matrix(GraphAlignment *ga){

}

int main (int argc, char *argv[]) {

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  static int D = 0;
  static int H = 1;
  static int V = 2;


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // command line
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Create new CmdLine object
  //----------------------------------------------------------------------------
  CmdLine cmd("", ' ', ProgramVersion);
    
  //----------------------------------------------------------------------------
  // add program-specific command line arguments
  //----------------------------------------------------------------------------

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
  arg.defaultValueString = "CTATTTTAGTAGTTGTTGTTA";
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
  arg.defaultValueString = "ATCGAAGATCCATGT";
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

  //----------------------------------------------------------------------------
  // add custom output handler
  //----------------------------------------------------------------------------
  MyOutput my;
  cmd.setOutput(&my);

  //----------------------------------------------------------------------------
  // parse command line and catch possible errors
  //----------------------------------------------------------------------------
  try {
    cmd.parse(argc,argv);
  } 
  catch ( ArgException& e ) { 
    cerr << "ERROR: " << e.error() << " " << e.argId() << endl; 
  }
  
  //----------------------------------------------------------------------------
  // assign command line parameters
  //----------------------------------------------------------------------------

  string subject = cmd_subject.getValue();
  string query = cmd_query.getValue();
  int M = cmd_M.getValue();
  int X = cmd_X.getValue();
  int GI = cmd_GI.getValue();
  int GE = cmd_GE.getValue();
  bool matrix = cmd_matrix.getValue();
  bool debug = cmd_debug.getValue();

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // check and fix command line options
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // report command line and parameters used
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
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

  ///  string ref = "CTATTTTAGTAGTTGTTGTTA";
  /// pair<string, string> sv = std::make_pair("GTT", "G");
  /// int pos = 11;

  string query1 = "ATCCAGTAATCCGGGATCCAT";
  string query2 = "ATCCAGTATCCGGGATCCAT";

  pair<string, string> sv1 = std::make_pair("AATCC", "A");
  pair<string, string> sv2 = std::make_pair("ATCC", "A");

  int pos = 7;


  Variant v1 = {query1, sv1, pos};
  Variant v2 = {query2, sv2, pos};

  vector<string> strings1 = getNodes(v1);
  vector<string> strings2 = getNodes(v2);

    //cout << "\n" << strings[0] << ", " << strings[1] << ", " << strings[2] << ", " << strings[3] << std::endl;
  
  vector<Node *> subjectNodes = buildDiamondGraph(strings1);
  vector<Node *> subjectNodes2 = buildDiamondGraph(strings2);


  GraphAlignment * ga;
  GraphAlignment * ga2;
  ga = new GraphAlignment(subjectNodes, query1, M, X, GI, GE, debug);
  ga2 = new GraphAlignment(subjectNodes2, query2, M, X, GI, GE, debug);
  

  for (vector<Node *>::const_iterator iter = subjectNodes.begin(); iter!= subjectNodes.end(); iter++){
    //ga->printMatrix(* iter, cout);
    //ga->printTBMatrix(* iter, cout);
  }
  
  Traceback t1 = {subjectNodes, ga};
  // vector<int**> TB1 = t1.buildTB();

  Traceback t2 = {subjectNodes2, ga2};
  //vector<int**> TB2 = t2.buildTB();

  vector<Traceback> tracebacks;
  tracebacks.push_back(t1);
  tracebacks.push_back(t2);
  

  PileUp p = {tracebacks};

  p.sumTracebacks();

  Graph g = refit(subjectNodes, ga, query, M, X, GI, GE, debug);
  ga = g.alignment;
  subjectNodes = g.nodes;

  ga = updateGA(subjectNodes, query, M, X, GI, GE, debug);
  
  Traceback TBAfter = {subjectNodes, ga};
  vector<int**> after = TBAfter.buildTB();

  cout << "Optimal score of GSW: " << ga->getScore() << endl;
  cout << "Global Cigar:" << ga->getGlobalCigar() << endl;
  cout << "Global Alignment:" << endl << ga->getGlobalAlignment() << endl;
  
  vector<Node *> matchedNodes = ga->getMatchedNodes();
  cout << "Graph node alignments:" << endl;
  for (vector<Node *>::const_iterator iter = matchedNodes.begin(); iter != matchedNodes.end(); iter++) {
    Node * node = * iter;
    string cigar = ga->getNodeCigar(node);
    int offset = ga->getNodeOffset(node);
    //cout << "  Node=" << node->getId() << " CIGAR=" << cigar << " offset=" << offset << endl;
    //ga->printMatrix(node, cout);
  }
}

