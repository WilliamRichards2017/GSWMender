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

vector<vector<Node *> >buildAllGraphs(vector<vector<string> > allStrings){
  vector<vector<Node *> > allGraphs;
  for(auto it = std::begin(allStrings); it != std::end(allStrings); ++it){
    vector<string> strings = *it;
    vector<Node *> graph = buildDiamondGraph(strings);
    allGraphs.push_back(graph);
  }
  return allGraphs;
}
