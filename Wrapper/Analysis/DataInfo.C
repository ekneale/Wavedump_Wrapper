#define DataInfo_C
#include "DataInfo.h"

int DataInfo::GetNSamples(char test,
			  char digitiser
			  ){
  
  switch(digitiser){     
  case ('V'):
    if (test=='A')
      return 5100;
    else
      return 110;
  case ('D'):
    return 1024;
  default:
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }
  
}
