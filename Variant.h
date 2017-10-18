#ifndef __Variant_H_INCLUDED__
#define __Variant_H_INCLUDED__

#include <string>

using namespace std;

struct Variant{
  string ref;
  std::pair<string, string> sv;
  int pos;

};

#endif //  __Variant_H_INCLUDED__
