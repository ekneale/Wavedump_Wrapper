#ifndef DataInfo_h
#define DataInfo_h

#include <iostream>

using namespace std;

class DataInfo {
 public :
  
  DataInfo(){}
  ~DataInfo(){}
  int   GetNSamples(char test = 'G',
		    char digitiser = 'V');
  int   GetNVDCBins(char digitiser = 'V'); 
  float GetVoltageRange(char digitiser = 'V');
  float GetnsPerSample(char digitiser = 'V');
  float GetmVPerBin(char digitiser = 'V');
  int   GetSampleRateInMHz(char digitiser = 'V');
  
};

#endif


#ifdef DataInfo_C 

int  DataInfo::GetNVDCBins(char digitiser){

  if     ( digitiser == 'V' )
    return 16384;
  else if( digitiser == 'D' )
    return 4096;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }
  
}

float DataInfo::GetVoltageRange(char digitiser){
  
  if     ( digitiser == 'V' )
    return 2.0;
  else if( digitiser == 'D' )
    return 1.0;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }

}

float DataInfo::GetnsPerSample(char digitiser){
  
  return (1.0e3 / (float)GetSampleRateInMHz(digitiser));

}

int DataInfo::GetSampleRateInMHz(char digitiser){

  if     ( digitiser == 'V' )
    return 500;
  else if( digitiser == 'D' )
    return 5000;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }

}

float DataInfo::GetmVPerBin(char digitiser){
  return ( 1.0e3 * GetVoltageRange(digitiser) / 
	   GetNVDCBins(digitiser) );
}

#endif 
