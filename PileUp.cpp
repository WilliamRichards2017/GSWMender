#include "PileUp.h"


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
  return strings;
}

vector<vector<string> > getAllNodes(){
  vector<vector<string> > allNodes;
  for(auto it = std::begin(variants_); it != std::end(variants_); ++it){
    Variant v = *it;
    vector<string> nodes = getNodes(v);
    allNodes.push_back(nodes);
  }
  return allNodes;
}

vector<Node *> buildDiamondGraph(vector<string> strings){
  vector<Node * > contributors1;
  Node * node1 = new Node(
			  "node1",
			  strings[0],
			  contributors1,
			  0
			  );
  subjectNodes_.push_back(node1);
  vector<Node *> contributors2;
  contributors2.push_back(node1);
  Node * node2 = new Node(
			  "node2",
			  strings[1],
			  contributors2,
			  1
			  );
  subjectNodes_.push_back(node2);
  vector<Node *> contributors3;
  contributors3.push_back(node1);
  Node * node3 = new Node(
			  "node3",
			  strings[2],
			  contributors3,
			      2
			  );
  subjectNodes_.push_back(node3);
  vector<Node *> contributors4;
  contributors4.push_back(node2);
  contributors4.push_back(node3);
  Node * node4 =  new Node(
			   "node4",
			   strings[3],
			   contributors4,
			        0
			   );
  subjectNodes_.push_back(node4);
  return subjectNodes;
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
	
void deleteGraph(){
  for(auto it = std::begin(subjectNodes_); it != std::end(subjectNodes_); ++it){
    delete * it;
  }
}

vector<vector<vector<int> > > sumTracebacks() {
  vector<vector<vector<int> > >  sumMatrix;
  vector<vector<string> > strings = getAllNodes();
  vector<vector<Node *> > nodes = buildAllGraphs(strings);
  int count = 0;
  for(auto it = std::begin(tbs); it != std::end(tbs); ++it){
    Traceback tb = *it;
    vector<vector<vector<int> > > matrices = tb.buildTB();
    vector<Node *> subjectNodes = tb._subjectNodes;
    unsigned c = 0;
    //iterate through dimensions vector to build up empty 2Ds                                                                        
    for(auto it = std::begin(matrices); it != std::end(matrices); ++it){
      vector<vector<int> > m = *it;
      vector<vector<int> > matrix = buildArray2D(m.size(), m[0].size());
      sumMatrix.push_back(matrix);
      cout << "built arrayy \n";
      for (unsigned i = 0; i < m.size(); i++){
	for(unsigned j = 0; j < m[0].size(); j++){
	  sumMatrix[c][i][j] += matrices[c][i][j];
	}
      }

      //cout << "printing out node " << c << std::endl;                                                                              
      printArray2D(sumMatrix[c]);
      c++;
    } // end of dims loop                                                                                                            
  } // end of traceback loop;                                                                                                        
  cout << "leacing sumTracebacks\n";
  return sumMatrix;
}
