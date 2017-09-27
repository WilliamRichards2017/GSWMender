#include "Class-GraphAlignment.h"

class Traceback {

 public:
  Traceback(
	    int* // pointer to array
	    );
  
  GraphAlignment ga;

  int *build_traceback_matrix(int* arrPointer);
}
