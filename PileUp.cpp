#include "PileUp.h"
#include "Traceback.h"
#include "ArrayUtil.h"


PileUp::PileUp(vector<Variant> variants) : variants_(variants){
  buildTracebackVector(variants_);
  sumMatrix_  = sumTracebacks();                                                       

}

vector<string> PileUp::getNodes(Variant v){
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

vector<vector<string> > PileUp::getAllNodes(){
  vector<vector<string> > allNodes;
  for(auto it = std::begin(variants_); it != std::end(variants_); ++it){
    Variant v = *it;
    vector<string> nodes = getNodes(v);
    allNodes.push_back(nodes);
  }
  return allNodes;
}

vector<Node *> PileUp::buildDiamondGraph(vector<string> strings){
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

vector<vector<Node *> > PileUp::buildAllGraphs(vector<vector<string> > allStrings){
  vector<vector<Node *> > allGraphs;
  for(auto it = std::begin(allStrings); it != std::end(allStrings); ++it){
    vector<string> strings = *it;
    vector<Node *> graph = buildDiamondGraph(strings);
    allGraphs.push_back(graph);
  }
  return allGraphs;
 }

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


 void PileUp::buildTracebackVector(vector<Variant> variants){
   cout << "size of variant vector is: " <<  variants.size() << std::endl;
   for(auto it = std::begin(variants); it != std::end(variants); ++it){
     vector<Node *> subjectNodes = buildDiamondGraph(getNodes(*it));
     cout << "size of subject nodes is " << subjectNodes.size() << std::endl;
     GraphAlignment * ga = new GraphAlignment(subjectNodes, it->ref, 2, -2, -3, -2, false);
     Traceback t = {subjectNodes, ga};
     tbs_.push_back(t);
   }
 }

 void PileUp::deleteGraph(){
   for(auto it = std::begin(subjectNodes_); it != std::end(subjectNodes_); ++it){
     delete * it;
   }
 }

vector<vector<vector<int> > > PileUp::sumTracebacks() {
  cout << "inside sumTracebacks\n";
  vector<vector<vector<int> > >  sumMatrix;
  vector<vector<string> > strings = getAllNodes();
  vector<vector<Node *> > nodes = buildAllGraphs(strings);
  int count = 0;
  for(auto it = std::begin(tbs_); it != std::end(tbs_); ++it){
    Traceback tb = *it;
    cout << "iterating through tracebacks\n";
    vector<vector<vector<int> > > matrices = tb.buildTBMs();
    cout << "size of tbms is " << matrices.size();
    vector<Node *> subjectNodes = tb.subjectNodes_;
    unsigned c = 0;
    //iterate through dimensions vector to build up empty 2Ds                                                                        
    for(auto it = std::begin(matrices); it != std::end(matrices); ++it){
      vector<vector<int> > m = *it;
      vector<vector<int> > matrix = ArrayUtil::buildArray2D(m.size(), m[0].size());
      if(count == 0){
	sumMatrix.push_back(matrix);
      }
      
      for (unsigned i = 0; i < m.size(); i++){
	for(unsigned j = 0; j < m[0].size(); j++){
	  sumMatrix[c][i][j] += matrices[c][i][j];
	  //cout <<  "c, i j, count = " << c << ", " << i << ", " << j << ", " << count << std::endl;
	}
      }
      
      cout << "printing out node " << c << std::endl;
      ArrayUtil::printArray2D(sumMatrix[c]);
      c++;
    } // end of dims loop                                                                                                            
    count++;
  } // end of traceback loop;                                                                                                        
  return sumMatrix;
}
