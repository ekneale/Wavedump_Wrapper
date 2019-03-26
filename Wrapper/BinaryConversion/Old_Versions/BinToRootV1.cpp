/***************************************************
 * A program to process wavedump output files
 *
 * Author 
 *  Gary Smith
 *  https://github.com/gsmith23
 *  28 11 18
 *
 * Adapted from
 * SPE_Gen.cpp 
 *  Tomi Akindele
 *  https://github.com/akindele
 *  04 10 2018
 *
 * Purpose
 *  This program reads the binary 
 *  file and writes a TTree
 *  and histograms to a root file. 
 *
 * How to build
 *  $ make BinToRoot
 *
 * How to run
 *  $ ./BinToRoot  
 *
 * Dependencies
 *  root.cern
 *  Makefile
 *
 */ 

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "TH2.h"

using namespace std;

void printDAQConfig(digitiser){
  
  cout << endl;
  cout << "   " 
       << getNVDCBins(digitiser)        
       << " VDC Bins "    << endl;
  cout << "   " 
       << getVoltageRange(digitiser)    
       << " Volts Range "  << endl;
  cout << "   " 
	 << getSampleRateInMHz(digitiser) 
       << " MHz sample rate " << endl;
  cout << endl;
}

int  getNSamples(char digitiser){

  if     ( digitiser == 'V' )
    return 110;
  else if( digitiser == 'D' )
    return 1024;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }
  
}

int  getNVDCBins(char digitiser){

  if     ( digitiser == 'V' )
    return 16384;
  else if( digitiser == 'D' )
    return 4096;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }
  
}

float getVoltageRange(char digitiser){
  
  if     ( digitiser == 'V' )
    return 2.0;
  else if( digitiser == 'D' )
    return 1.0;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }


}

int getSampleRateInMHz(char digitiser){

  if     ( digitiser == 'V' )
    return 500;
  else if( digitiser == 'D' )
    return 5000;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }

}


float getnsPerSample(char digitiser){
  
  return (1.0e3 / (float)getSampleRateInMHz(digitiser));

}

float getDelay(char digitiser){
  
  if      ( digitiser == 'V' )
    return 54.;
  else if ( digitiser == 'D' )
    return 105.;
  else 
    return 0.;
}

float getmVPerBin(char digitiser){
  return ( 1.0e3 * getVoltageRange(digitiser) / 
	   getNVDCBins(digitiser) );
}

float getCharge(int intVDC, char digitiser = 'V'){
  
  float nsPerSample = getnsPerSample(digitiser);
  
  float mVolts = (float)intVDC * getmVPerBin(digitiser);
  
  float charge = -1. * nsPerSample * mVolts;
  
  if( digitiser == 'V' || digitiser == 'D' )
    return charge; 
  else
    return 0.;

}


bool isCorrectDigitiser(int header, char digitiser){
  
  if( ( header == 244  &&  digitiser == 'V' ) || 
      ( header == 4120 &&  digitiser == 'D' ) )
    return true;
  else{
    cerr << " Error: digitiser choice does not match header info " << endl;
    return false;
  }
}

// integration windows fixed wrt trigger
// pedestal window before signal window only
float Accumulate_V1(short VDC, float time){

  // Integrate pedestal using 
  // 40 ns window before signal
  // And signal in 40 ns window
  if      ( time >= -40 && 
	    time <    0 )
    return((int)-VDC);
  else if ( time >= 0 &&
	    time <  40 ){
    return((int)VDC);
  }
  else
    return 0.;
  
}

// integration window fixed wrt trigger
// pedestal windows before and after signal window
float Accumulate_V2(short VDC, float time){

  // Integrate pedestal using 
  // 20 ns window before signal
  // and 20 ns window after
  // And signal in 40 ns window
  if      ( time >= -20 && 
	    time <    0 )
    return((float)-VDC);
  else if ( time >= 0 &&
	    time < 40  )
    return((float)VDC);
  else if ( time >= 40 &&
	    time <  60 )
    return((float)-VDC);
  else
    return 0.;
   
}

// integration windows wrt pulse peak
// -10 ns before to 20 ns after
// (timeRel is time relative to minT)
// pedestal window before
float Accumulate_V3(short VDC, float timeRel){
  
  if      ( timeRel >= -50 && 
	    timeRel <  -10 )
    return((float)-VDC);
  else if ( timeRel >= -10 &&
	    timeRel <   30 ){
    return((float)VDC);
  }
  else
    return 0.;

}

// integration window wrt pulse peak
// pedestal windows before and after signal window
// ( double counting bug as in SPE_Gen )
float Accumulate_V4(short VDC, 
		    short sample,
		    short minT,
		    char  digitiser){
  int w = 0;
  
  if     ( digitiser == 'V' )
    w = 10;
  else if( digitiser == 'D')
    w = 100;
    
  float result = 0;
  
  if ( sample >= (minT - 3*w) && 
       sample <= (minT - 2*w) )
    result += ((float)-1*VDC)*3./4.;
  
  if ( sample >= (minT - 2*w) && 
       sample <= (minT - 1*w) )
    result += ((float)-1*VDC)*3./4.;
  
  if ( sample >= (minT - 1*w) &&
       sample <= (minT      ) )
    result += (float)VDC;
  
  if ( sample >= (minT      ) &&
       sample <= (minT + 1*w) )
    result += (float)VDC;
  
  if ( sample >= (minT + 1*w) &&
       sample <= (minT + 2*w) )
    result += (float)VDC;
  
  if ( sample >= (minT + 2*w) &&
       sample <= (minT + 3*w) )
    result += ((float)-1*VDC)*3./4.;
  
  if ( sample >= (minT + 3*w) && 
       sample <= (minT + 4*w) )
    result += ((float)-1*VDC)*3./4.;
  
  return result;
  
}

// Read in binary file and 
// write to root file
int ProcessBinaryFile(string filePath,
		      int    verbosity = 0,
		      char   digitiser = 'V'
		      ){

  cout << " Processing  " << filePath << endl;

  bool  testMode  = false;
  bool  keepGoing = true;
  int   maxEvents = 50000;
  
  if     ( verbosity == 1 )
    maxEvents = 10;
  else if( verbosity == 2 )
    maxEvents = 1;

  // Read from here
  ifstream fileStream(filePath);

  if(!fileStream.good()){
    cerr << endl;
    cerr << " check filename " << endl;
    return -1;
  }
  
  printDAQConfig(digitiser);
    
  // Write to here
  TFile * outFile    = new TFile("outputFile.root",
				 "RECREATE",
				 "Wavedump Data");
  //---------------------
  // Event Level Data
  TTree * eventTree  = new TTree("eventTree",
				 "event-level variables");
  
  // accumulated and coverted data
  TH1F * hCharge1 = new TH1F("hCharge1","hCharge1;Charge (mV nS);Counts",
			     128,-560.,2000.);
  
  TH1F * hCharge2 = new TH1F("hCharge2","hCharge2;Charge (mV nS);Counts",
			     128,-560.,2000.);
  
  TH1F * hCharge3 = new TH1F("hCharge3","hCharge3;Charge (mV nS);Counts",
			     128,-560.,2000.);

  TH1F * hCharge4 = new TH1F("hCharge4","hCharge4;Charge (mV nS);Counts",
			     128,-560.,2000.);
  
  TH2F * hPulses = new TH2F("hPulses",
			    "hPulses;Sample;VDC",
			    getNSamples(digitiser),0,(getNSamples(digitiser)-1),
			    getNVDCBins(digitiser),0.,(getNVDCBins(digitiser)-1));
  
  TH2F * hTV = new TH2F("hTV",
			"hTV;Time (ns);Voltage (mV)",
			getNSamples(digitiser),
			0.,(getNSamples(digitiser)-1)*getnsPerSample(digitiser),
			getNVDCBins(digitiser),
			0.,getVoltageRange(digitiser)*1.0e3);

  int   event  = -1;
  
  // Note that pulses are negative polarity and
  // the sign is preserved here, therefore 
  // minima will be where the signal was largest
  short minVDC = 32767, maxVDC = -32768;  
  short minT   = 32767, maxT   = -32768;  
  int   intVDC1 = 0,  intVDC2 = 0,  intVDC3 = 0,  intVDC4 = 0;
  float charge1 = 0., charge2 = 0., charge3 = 0., charge4 = 0.;  

  // read in samples per event
  // vector may be better
  short pulse[getNSamples(digitiser)];
  
  eventTree->Branch("event",&event,"event/I");
  eventTree->Branch("minVDC",&minVDC,"minVDC/S");
  eventTree->Branch("maxVDC",&maxVDC,"maxVDC/S");
  eventTree->Branch("minT",&minT,"minT/S");
  eventTree->Branch("maxT",&maxT,"maxT/S");
  
  TString arrayString = "";
  arrayString.Form("pulse[%d]/S",getNSamples(digitiser));
  
  eventTree->Branch("pulse",pulse,arrayString);
  
  //---------------------
  // Sample Level Data

  short VDC = 0, sample = 0;
  
  int   fileHeader = 0;
  float floatVDC   = 0.;
  
  // Time (delay subtracted)
  float time = 0.;
  
  // read in data from streamer object
  // until the end of the file
  while ( fileStream.is_open() && 
	  fileStream.good()    && 
	  !fileStream.eof()    &&
	  keepGoing 
	  ){
    
    //-------------------
    // file-level data
    event++;
    
    //-------------------
    // event-level data
    intVDC1 = 0, intVDC2 = 0, intVDC3 = 0, intVDC4 = 0;
    charge1 = 0.,charge2 = 0.,charge3 = 0.,charge4 = 0.;
    
    // VDC range
    minVDC =  32767;
    maxVDC = -32768;  

    // time walk?
    minT   =  32767;
    // random?
    maxT   = -32768;

    //-------------------
    // pulse-level data
    sample = 0;
    VDC    = 0;    
    
    // read in header info which comes 
    // as six four bit sized chunks
    for (int iHeader = 0 ; iHeader < 6 ; iHeader++ ){
      
      fileHeader = 0;
      
      // read in 4 bit segment
      fileStream.read( (char*)&fileHeader,
		       sizeof(int)); // read 4 bits
      
      // check first header value matches expectations
      // NB other values may be acceptable so modify
      // isCorrectDigitiser() as appropriate
      if ( iHeader == 0 && 
	   event   == 0 &&
	   !isCorrectDigitiser(fileHeader,
			       digitiser)  ) {
	
	return -1;
      }
    } // end: for (int i = 0 ; i < intsPerHeader
    
    // read in pulse which comes 
    // in 2 (VME) or 4 (Desktop) bit chunks
    for (short iSample = 0; iSample < getNSamples(digitiser); iSample++){
      
      VDC = 0;
      
      // read 2 bits
      if     ( digitiser == 'V' ){
	fileStream.read((char*)&VDC,2); 
      }
      // read 4 bits
      else if( digitiser == 'D' ){
	fileStream.read((char*)&floatVDC,sizeof(float));
	VDC = (short)floatVDC;
      }
      else {
	
	cerr << " Error:  incorrect digitiser " << endl;
	return -1;
	
      }
      
      // for writing
      sample = iSample;
      pulse[sample] = VDC;
      
      if     ( VDC < minVDC ){
	minVDC = VDC;
	minT   = sample;
      }
      else if( VDC > maxVDC ) {
	maxVDC = VDC;
	maxT   = sample;
      }
      
      // use time with delay subtracted so that 
      // any runs can use same windows
      time = sample*getnsPerSample(digitiser) - getDelay(digitiser);
      
      // fixed window accumulations
      // add / subtract / skip sample VDC value
      intVDC1 += Accumulate_V1(VDC,time);

      intVDC2 += Accumulate_V2(VDC,time);
      
      //--------------------------------
      // Sample by sample data here
      
      if(verbosity > 1)
	cout << " VDC(" << iSample << ") = " << VDC << endl;
      
    } // end: for (short iSample = 0; iSa
    
    float minTime = (float)minT*getnsPerSample(digitiser);
    float timeMinusMinTime = 0;
    
    // Loop over pulse again
    // necessary for when event-level variables are used
    // such as pulse minimum
    for (short iSample = 0; iSample < getNSamples(digitiser); iSample++){
      
      time = (float)iSample*getnsPerSample(digitiser) - getDelay(digitiser);

      timeMinusMinTime = time - minTime + getDelay(digitiser);;
	
      if(verbosity > 0){
	cout << " minTime = " << minTime << endl;
	cout << " time    = " << time << endl;
      }
      
      if( minTime > (getDelay(digitiser) - 20.)  && 
	  minTime < (getDelay(digitiser) + 20.) ){
	
	// accumulations wrt pulse peak (minimum)
	intVDC3 += Accumulate_V3(pulse[iSample],
				 timeMinusMinTime);
      }
      
      if( minTime > 60. && 
	  minTime < 120. ){

	intVDC4 += Accumulate_V4(pulse[iSample],
				 iSample,
				 minT,
				 digitiser);
      }
      
      if( !(event%1000) && event < 100000 ){
	hPulses->Fill(iSample,pulse[iSample]);
	hTV->Fill(time + getDelay(digitiser),
		  pulse[iSample]*getmVPerBin(digitiser));
      }
      
    } // end: for (short iSample = 0; iS...
    
    if( minVDC < 0 && maxVDC > 0 )
      cout << " Warning: pulse is zero crossing " << endl;
    
    // convert accumulated VDC to charge
    charge1 = getCharge(intVDC1,digitiser);
    charge2 = getCharge(intVDC2,digitiser);
    
    if( verbosity > 0 ){
      cout << endl;
      cout << " charge1 = " << charge1 << endl;
      cout << " charge2 = " << charge2 << endl;
    }
    
    hCharge1->Fill(charge1);
    hCharge2->Fill(charge2);
    
    if( minTime > (getDelay(digitiser) - 20.) && 
	minTime < (getDelay(digitiser) + 20.) ){
      
      charge3 = getCharge(intVDC3,digitiser);
      
      if( verbosity > 0 )
	cout << " charge3 = " << charge3 << endl;
      
    }
    
    if( minTime > 60. && 
	minTime < 120. ){
      
      charge4 = getCharge(intVDC4,digitiser);
      
      if( verbosity > 0 )
	cout << " charge4 = " << charge4 << endl;
      
      hCharge3->Fill(charge3);
      hCharge4->Fill(charge4);      
    }
    
    //--------------------------------
    // Write event by event data here
    
    eventTree->Fill();
    
    if(verbosity > 0){
      
      cout << endl;
      cout << " Event " << event << endl;
      cout << " minVDC(" << minT << ") = " << minVDC << endl;
      cout << " maxVDC(" << maxT << ") = " << maxVDC << endl;
      cout << endl;
      
      if(verbosity > 1){
	
	cout << endl;
      }
    
    }
    
    if( testMode &&
	(event+1) == maxEvents )
      keepGoing = false;
    
  } // end: while loop

  // close wavedump file
  fileStream.close();	

  //--------------------------------
  // Write file info here

  eventTree->Write();
  eventTree->Delete();

  outFile->Write();
  outFile->Close();
  
  return (event+1);
}


string getFilePath(){
  return  "../../Data/wave_0.dat";
}


int main(int argc, char **argv)
{
  
  // ------------------
  // Optional variables
  // - see ProcessBinaryFile()
  
  // 0 - silence (default) 
  // 1 - event-by-event
  // 2 - sample-by-sample
  int  verbosity   = 0;
  
  // Printing
  // 'V' - VME (default)
  // 'D' - Desktop
  char   digitiser = 'V';
  digitiser = 'D';
  
  //-------------------
  
  string filePath = getFilePath();
  
  int nEvents = ProcessBinaryFile(filePath,
				  verbosity,
				  digitiser);
  
  cout << " Output file contains " << nEvents << " events " << endl; 
}
