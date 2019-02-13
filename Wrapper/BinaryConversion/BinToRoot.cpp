/***************************************************
 * A program to process wavedump output files
 *
 * Author 
 *  gary.smith@ed.ac.uk
 *  28 11 18
 *
 * Adapted from
 * SPE_Gen.cpp 
 *  Tomi Akindele
 *  https://github.com/akindele
 *  04 10 2018
 *
 * Purpose
 *  This program reads a wavedump binary 
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
 *
 * Modified 
 *  gary.smith@ed.ac.uk
 *  11 12 18
 *
 */ 

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TLine.h"

using namespace std;

int  getNSamples(char digitiser,
		 char test = 'G'){

  switch(digitiser){     
  case ('V'):
    if (test=='A')
      return 5100;
    else
      return 110;
  case ('D'):
    return 1024;
  case ('d'):
    return 1024;
  default:
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }
  
}

int  getNVDCBins(char digitiser){

  if     ( digitiser == 'V' )
    return 16384;
  else if( digitiser == 'D' ||
	   digitiser == 'd' )
    return 4096;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }
  
}

float getVoltageRange(char digitiser){
  
  if     ( digitiser == 'V' )
    return 2.0;
  else if( digitiser == 'D' ||
	   digitiser == 'd' )
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
  else if( digitiser == 'd' )
    return 1000;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }

}

void printDAQConfig(char digitiser){
  

  TString strDigitiser = "VME";
  
  if (digitiser=='d')
    strDigitiser = "Desktop";

  cout << endl;
  cout << "   " 
       << strDigitiser
       << " Digitiser "    << endl;
  cout << "   " 
       << getNVDCBins(digitiser)        
       << " VDC Bins "    << endl;
  cout << "   " 
       << getVoltageRange(digitiser)    
       << " Volts Range "  << endl;
  cout << "   " 
	 << getSampleRateInMHz(digitiser) 
       << " MHz sample rate " << endl;

}

float getnsPerSample(char digitiser){
  
  return (1.0e3 / (float)getSampleRateInMHz(digitiser));

}

float getPulseLength(char digitiser,
		     char test){

  return (float)getNSamples(digitiser,test)*getnsPerSample(digitiser);

}

float getDelay(int run = 0){
  
  if      ( run == 0 ) // 12th October
    return 60.; 
  else if ( run < 10 ) // 16th October
    return 50.;
  else if ( run < 20 ) // Underground
    return 90.;
  else                 // Surface
    return 60.;
}

float getmVPerBin(char digitiser){
  return ( 1.0e3 * getVoltageRange(digitiser) / 
	   getNVDCBins(digitiser) );
}

float getCharge(int intVDC, char digitiser = 'V'){
  
  float nsPerSample = getnsPerSample(digitiser);
  
  float mVolts = (float)intVDC * getmVPerBin(digitiser);
  
  float charge = -1. * nsPerSample * mVolts;
  
  if( digitiser == 'V' || 
      digitiser == 'D' )
    return charge; 
  else
    return 0.;

}


bool isCorrectDigitiser(int header,
			char digitiser,
			int test){
  
  if( digitiser == 'd' )
    digitiser = 'D';

  switch( test ){
  case ('A'): 
    if( header == 10224 && digitiser == 'V')
      return true;
    else{
      cerr << " Error: digitiser choice does not match header info " << endl;
      return false;
    }
  default:
    if( ( header == 244  &&  digitiser == 'V' ) || 
	( header == 4120 &&  digitiser == 'D' ) )
      return true;
    else{
      cerr << " Error: digitiser choice does not match header info " << endl;
      return false;
    }
  }
}

// Must leave room for pedestal window
bool isPeakInRange(float minTime, char digitiser,
		   char test){
				 
  if( minTime >= 50.   && 
      minTime <  (getPulseLength(digitiser,test) - 30.0)  )
    return true;
  else
    return false;
}

float gateWidth(){
  return 50;
}

// integration windows fixed wrt trigger
// baseline window before signal window only
float Accumulate_Fixed(short VDC, float time){

  // Integrate baseline using 
  // 50 ns window before signal
  // (max size given run 1 delay)
  // And signal in 50 ns window
  // NB using factors of two
  // given VME sampling rate.
  // Using greater than 40 to 
  // allow for trigger jitter
  if      ( time >= -gateWidth()  && 
	    time <    0 )
    return((int)-VDC);
  else if ( time >= 0 &&
	    time <  gateWidth() ){
    return((int)VDC);
  }
  else
    return 0.;
  
}

// Integration windows wrt pulse peak
// -10 ns before to 20 ns after
// timeRel is time relative to minT (in ns)
// baseline window before
float Accumulate_Peak(short VDC, float timeRel){
  
  if      ( timeRel >= -50 && 
	    timeRel <  -10 )
    return((int)-VDC);
  else if ( timeRel >= -10 &&
	    timeRel <   30 ){
    return((int)VDC);
  }
  else
    return 0.;

}


TString getRunFolderName(int run){

  if     ( run < 10  )
    return "RUN00000%d/";
  else if( run < 100 )
    return "RUN0000%d/";
  else if( run < 1000 )
    return "RUN000%d/";
  else if( run < 10000 )
    return "RUN00%d/";
  else if( run < 100000 )
    return "RUN0%d/";
  else 
    return "RUN%d/";
}

TString getPMTFolderName(int pmt){

  if     ( pmt < 10  )
    return "PMT000%d/";
  else if( pmt < 100 )
    return "PMT00%d/";
  else if( pmt < 1000 )
    return "PMT0%d/";
  else
    return "PMT%d/";
}

TString getTestFolderName(char test){
  
  switch(test){
  case('S'):
    return "SPEtest/";
  case('N'):
    return "Nominal/";
  case('D'):
    return "DarkRateTest/";
  case('A'):
    return "APTest/";
  case('G'):
    return "GainTest/";
  default:
    cerr << " Error: invalid test type " << endl;
    return "";
  }
}


TString getBinaryFilePath(TString filePath = "../../Data/",
			  int  run  = 1, 
			  int  pmt  = 1, 
			  int  loc  = 0,
			  char test = 'S',
			  int  hvStep = 0){
  
  TString rtnFilePath = "";
  
  // append file path
  filePath += getRunFolderName(run);
  filePath += getPMTFolderName(pmt);
  filePath += getTestFolderName(test);
  
  // append with filename
  switch(test){
  case('G'):
    filePath += "wave_%d_hv%d.dat";
    rtnFilePath.Form(filePath,run,pmt,loc,hvStep);
    break;
  default:
    filePath += "wave_%d.dat";
    rtnFilePath.Form(filePath,run,pmt,loc);
  }

  return rtnFilePath;
}

TString getRawRootFilePath(TString filePath = "./",
			   int  run  = 0, 
			   int  pmt  = 1, 
			   int  loc  = 0,
			   char test = 'S',
			   int  hvStep = 0){
  
  // default option
  if( run == 0 )
    return "outputFile.root";

  // Use same folder structure as
  // binary files ( folder creation
  // not implemented )
  // or write to single folder
  bool singleOutputFolder = true;

  TString rtnFilePath = "";

  if ( singleOutputFolder ){
    // append with filename
    // containing info
    switch(test){
    case('G'):
      filePath += "Run_%d_PMT_%d_Loc_%d_HV_%d";
      rtnFilePath.Form(filePath,run,pmt,loc,hvStep);
      break;
    default:
      filePath += "Run_%d_PMT_%d_Loc_%d_Test_%c";
      rtnFilePath.Form(filePath,run,pmt,loc,test);
    }
  }
  else {
    
    filePath += getRunFolderName(run);
    filePath += getPMTFolderName(pmt);
    filePath += getTestFolderName(test);
    
    // append with filename
    // containing info
    switch(test){
    case('G'):
      filePath += "Run_%d_PMT_%d_Loc_%d_HV_%d";
      rtnFilePath.Form(filePath,run,pmt,run,pmt,loc,hvStep);
      break;
    default:
      filePath += "Run_%d_PMT_%d_Loc_%d_Test_%c";
      rtnFilePath.Form(filePath,run,pmt,run,pmt,loc,test);
    }
  }

  return rtnFilePath;
  
}

// Read in binary file and 
// write to root file
int ProcessBinaryFile(TString inFilePath,
		      TString outFilePath,
		      int run = 0, int loc = 0,
		      int pmt = 1, int hvStep = 0,
		      char test = 'S',
		      int  verbosity = 0,
		      char digitiser = 'V'
		      ){

  inFilePath = getBinaryFilePath(inFilePath,
				 run, pmt, loc, test, hvStep);
  
  cout << endl;
  cout << " Processing  " << inFilePath << endl;

  //----------------------
  // Variables for testing
  bool  testMode  = false;
  bool  keepGoing = true;
  int   maxEvents = 100;
  
  if ( test == 'A' )
    maxEvents = 2;
  
  if     ( verbosity == 1 )
    maxEvents = 10;
  else if( verbosity == 2 )
    maxEvents = 1;
  //----------------------


  // Read from here
  ifstream fileStream(inFilePath);

  if(!fileStream.good()){
    cerr << endl;
    cerr << " Error: check filename " << endl;
    return -1;
  }
  
  if( verbosity > 0 ) 
    printDAQConfig(digitiser);
  
  outFilePath = getRawRootFilePath(outFilePath,
				   run, pmt, loc, test, hvStep);

  
  TString outFileName = outFilePath;
  outFilePath += ".root";
  
  cout << endl;
  cout << " Writing     " << outFilePath << endl;
  
  // Write to here
  TFile * outFile    = new TFile(outFilePath,
				 "RECREATE",
				 outFileName);

  TString eventTreeNameTemp = "", eventTreeName = "";;
  TString canvasNameTemp    = "", canvasName    = "";
  
  TString hQ_FixedNameTemp = "", hQ_FixedName = "";
  TString hQ_PeakNameTemp  = "", hQ_PeakName  = "";
  
  if(hvStep == 0){
    eventTreeNameTemp = "Events_Run_%d_PMT_%d_Loc_%d_Test_%c";
    eventTreeName.Form(eventTreeNameTemp,run,pmt,loc,test);
    canvasNameTemp = "Canvas_Run_%d_PMT_%d_Loc_%d_Test_%c";
    canvasName.Form(canvasNameTemp,run,pmt,loc,test);
    
    hQ_FixedNameTemp = "hQ_Fixed_Run_%d_PMT_%d_Loc_%d_Test_%c";
    hQ_FixedName.Form(hQ_FixedNameTemp,run,pmt,loc,test);
    
    hQ_PeakNameTemp = "hQ_Peak_Run_%d_PMT_%d_Loc_%d_Test_%c";
    hQ_PeakName.Form(hQ_PeakNameTemp,run,pmt,loc,test);
    
  }
  else{
    eventTreeNameTemp = "Events_Run_%d_PMT_%d_Loc_%d_HV_%d";
    eventTreeName.Form(eventTreeNameTemp,run,pmt,loc,hvStep);
    canvasNameTemp = "Canvas_Run_%d_PMT_%d_Loc_%d_HV_%d";
    canvasName.Form(canvasNameTemp,run,pmt,loc,hvStep);

    hQ_FixedNameTemp = "hQ_Fixed_Run_%d_PMT_%d_Loc_%d_HV_%d";
    hQ_FixedName.Form(hQ_FixedNameTemp,run,pmt,loc,hvStep);
    
    hQ_PeakNameTemp = "hQ_Peak_Run_%d_PMT_%d_Loc_%d_HV_%d";
    hQ_PeakName.Form(hQ_PeakNameTemp,run,pmt,loc,hvStep);

  }
  //---------------------
  // Event Level Data
  TTree * eventTree  = new TTree(eventTreeName,
				 eventTreeName);
  
  TCanvas * canvas = new TCanvas(canvasName,
				 canvasName);
  
  // increase binning by 4 times
  Int_t nBinsX = 512;
  Float_t rangeX[2] = {-500.,2000.};

  TH1F * hQ_Fixed = new TH1F(hQ_FixedName,
			     "Fixed gate;Charge (mV nS);Counts",
			     nBinsX,rangeX[0],rangeX[1]);
  
  TH1F * hQ_Peak = new TH1F(hQ_PeakName,
			    "Gate around peak;Charge (mV nS);Counts",
			    nBinsX,rangeX[0],rangeX[1]);
  
  TH2F * hPulses = new TH2F("hPulses",
			    "Subset of Raw Pulses;Sample;VDC",
			    getNSamples(digitiser,test),0,
			    (getNSamples(digitiser,test)-1),
			    getNVDCBins(digitiser),0.,
			    (getNVDCBins(digitiser)-1));
  
  TH2F * hTV = new TH2F("hTV",
			"Subset of calibrated pulses;Time (ns);Voltage (mV)",
			getNSamples(digitiser,test),
			0.,(getNSamples(digitiser,test)-1)*getnsPerSample(digitiser),
			getNVDCBins(digitiser),
			0.,getVoltageRange(digitiser)*1.0e3);
  int   event  = -1;

  // Note that pulses are negative polarity and
  // the sign is preserved here, therefore 
  // minima will be where the signal peaked
  short minVDC = 32767, maxVDC = -32768;  
  short minT   = 32767, maxT   = -32768;  
  int   intVDCfixed = 0, intVDCpeak = 0;
  float qFixed = 0., qPeak = 0.;

  // read in samples per event
  // vector may be better
  short pulse[getNSamples(digitiser,test)];
  
  eventTree->Branch("event",&event,"event/I");
  eventTree->Branch("minVDC",&minVDC,"minVDC/S");
  eventTree->Branch("maxVDC",&maxVDC,"maxVDC/S");
  eventTree->Branch("minT",&minT,"minT/S");
  eventTree->Branch("maxT",&maxT,"maxT/S");
  
  TString arrayString = "";
  arrayString.Form("pulse[%d]/S",getNSamples(digitiser,test));
  
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
    
    if( (event > 0) && 
	((event % 1000 == 0 && event < 5000) ||
	 (event % 1000000 == 0))              ){
      cout << endl;
      cout << " event count " << event << endl;;
    }

    //-------------------
    // event-level data
    intVDCfixed = 0, intVDCpeak = 0;
    qFixed = 0.,qPeak = 0.;
    
    // VDC range
    minVDC =  32767;
    maxVDC = -32768;  

    // time walk?
    minT   =  32767;
    // random?
    maxT   = -32768;

    //-------------------
    // sample-level data
    sample = 0;
    VDC    = 0;    
    
    if(event < 10 && 
       verbosity > 0){
      cout << endl;
      cout << " event      = " << event << endl;
    }
    
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
			       digitiser,test)  
	   ) {
	
	return -1;
      }
      
      if(event < 10 && 
	 verbosity > 0){
	cout << " fileHeader = " << fileHeader << endl;
      }
      
    } // end: for (int i = 0 ; i < intsPerHeader
    
    // read in pulse which comes 
    // in 2 (VME) or 4 (Desktop) bit chunks
    for (short iSample = 0; iSample < getNSamples(digitiser,test); iSample++){
      
      VDC = 0;
      
      // read 2 bits
      if     ( digitiser == 'V' ){
	fileStream.read((char*)&VDC,2); 
      }
      // read 4 bits
      else if( digitiser == 'D' ||
	       digitiser == 'd'){
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
      // any run can use same fixed windows
      time = sample * getnsPerSample(digitiser) - getDelay(run);
      
      // fixed window accumulations
      // add / subtract / skip sample VDC value
      intVDCfixed += Accumulate_Fixed(VDC,time);

      if(verbosity > 1)
	cout << " VDC(" << iSample << ") = " << VDC << endl;
      
    } // end: for (short iSample = 0; iSa
    
    // time of the signal peak
    float minTime = minT*getnsPerSample(digitiser);
    float timeRelPeak = 0;
    
    // Loop over pulse again
    // necessary for when event-level variables are used
    // such as pulse minimum
    for (short iSample = 0; iSample < getNSamples(digitiser,test); iSample++){
      
      time = iSample*getnsPerSample(digitiser) - getDelay(run);

      timeRelPeak = time - minTime + getDelay(run);
	
      // accumulations wrt pulse peak (minimum)      
      if( isPeakInRange(minTime,digitiser,test) )
	intVDCpeak += Accumulate_Peak(pulse[iSample],
				      timeRelPeak);
      
      // sample vs VDC and time vs voltage plots
      // for checking signals (delay etc) 
      // plot pulses for 1000 events
      // from first 100,000
      if( ( event < 100000    )  &&
	  ( event % 100  == 0 )
	  ){
	
	hPulses->Fill(iSample,pulse[iSample]);
	
	hTV->Fill(time + getDelay(run),
		  pulse[iSample]*getmVPerBin(digitiser));
	
      }
      
    } // end: for (short iSample = 0; iS...
    
    hQ_Fixed->Fill(getCharge(intVDCfixed,digitiser));
    
    if( isPeakInRange(minTime,digitiser,test) )
      hQ_Peak->Fill(getCharge(intVDCpeak,digitiser));
    
    
    //--------------------------------
    // Write event by event data here
    
    eventTree->Fill();
    
    if(verbosity > 0){
      cout << endl;
      cout << " Event " << event << endl;
      cout << " minVDC(" << minT << ") = " << minVDC << endl;
      cout << " maxVDC(" << maxT << ") = " << maxVDC << endl;
      cout << " qFixed = " << qFixed << endl;
      cout << " qPeak = " << qPeak << endl;
      cout << endl;
    }
    
    if( testMode &&
	(event+1) == maxEvents )
      keepGoing = false;
    
  } // end: while loop

  // close wavedump file
  fileStream.close();	

  canvas->Divide(2,2);

  canvas->cd(1);
  
  hQ_Fixed->SetAxisRange(-500., 2500.,"X");

  gPad->SetLogy(1);
  hQ_Fixed->Draw();
  
  canvas->cd(2);
  gPad->SetLogy(1);
  
  hQ_Peak->SetAxisRange(-500., 2500.,"X");
  hQ_Peak->Draw();
  
  canvas->cd(3);
  
  float minY = getVoltageRange(digitiser)*(16 - 2)/32*1.0e3;
  float maxY = getVoltageRange(digitiser)*(16 + 1)/32*1.0e3 ;

  hTV->SetAxisRange(minY,maxY,"Y");

  float minX = 0.;
  float maxX = 500.; // 220.
  
  if(test=='A'){
    hTV->SetAxisRange(minX,maxX,"X");
  }
  
  hTV->Draw("colz");
  
  float lineYMin = minY * (16 - 1)/(16 - 2); 
  float lineYMax = maxY * (16 + 0.25)/(16 + 1); 
  
  float lineXMin = getDelay(run) - gateWidth();
  float lineXMax = getDelay(run) + gateWidth();
  
  TLine *lPedMin = new TLine(lineXMin,lineYMin,
			     lineXMin,lineYMax); 
  lPedMin->SetLineColor(kRed);
  
  TLine *lSigMin = new TLine(getDelay(run),lineYMin,
			     getDelay(run),lineYMax); 
  lSigMin->SetLineColor(kBlue);
  
  TLine *lSigMax = new TLine(lineXMax,lineYMin,
			     lineXMax,lineYMax); 
  lSigMax->SetLineColor(kBlue);
  
  lPedMin->Draw("same");
  lSigMin->Draw("same");
  lSigMax->Draw("same");
  
  canvas->cd(4);
  
  minY = getNVDCBins(digitiser)*10./16 ;
  maxY = getNVDCBins(digitiser)*6./16 ;
  
  //hPulses->SetAxisRange(minY,maxY,"Y");
  hPulses->SetAxisRange(7500.,8500.,"Y");
    
  hPulses->Draw("colz");
  
  canvasName = "./Plots/";
  canvasName += canvas->GetName();
  canvasName += ".png";

  canvas->SaveAs(canvasName);
  
  //--------------------------------
  // Write file info here

  hPulses->Delete();
  hTV->Delete();

  eventTree->Write();
  eventTree->Delete();

  outFile->Write();
  outFile->Close();
  
  return (event+1);
}


string getDefaultFilePath(){
  return  "../../Data/wave_0.dat";
}

int main(int argc, char **argv)
{
  
  // ------------------
  // Optional variables
  // - see ProcessBinaryFile()

  // Printing  
  // 0 - silence (default) 
  // 1 - event-by-event
  // 2 - sample-by-sample
  int  verbosity   = 0;
  
  // 'V' - VME (default)
  // 'D' - Desktop
  char   digitiser = 'd';

  //-------------------
  int  run = 160000; 
  int  pmt = 16;
  int  loc = 0;
  char test = 'D';

  // Default - not gain test setting
  // Set automatically for test type
  // in loop below.
  int  hvStep = 0; 
  int  nSteps = 1;

  static const int nTests = 1;

  // 'S' SPE, 'G' Gain, 'D' Dark
  // 'A' After, 'N' Nominal, 
  char testList[nTests] = {'D'};
  //char testList[nTests] = {'S','N','G','D','A'};
  
  static const int nRuns = 1; // 5
  int  runList[nRuns] = {160001};
  
  // PMTS 130 131 132 133
  //int  runList[nRuns] = {2,3,10,22,23}; 
  //int  runList[nRuns] = {4,11,12,20,21};
  

  static const int nPMTsA = 1;//80;
  static const int nPMTsB = 0;//20;
  static const int nPMTs  = nPMTsA + nPMTsB;
  
  //int  pmtAList[nPMTsA] = {130,131,132,133};  
  //int  pmtAList[nPMTsA] = {90,159,166,171};  

  int  pmtAList[1] = {16};
  // 
  int  pmtBList[1] = {0}; 
  
  // //-------------
  // // RUN 1
  // // PMT 139 missing SPE data
  // int  pmtAList[nPMTsA] = {83 , 88,108,107,
  // 			   73 , 76, 84, 87,
  // 			   66 , 78, 82,103,
  // 			   104,106,112,141,
  // 			   61 , 65, 75,105,
  // 			   74 ,111,140,142,
  // 			   143,145,146,147,
  // 			   63 , 67,158,160,
  // 			   139,161,164,165,
  // 			   90 ,159,166,171,
  // 			   81 ,167,169,170,
  // 			   50 , 53,162,163,
  // 			   55 , 56, 92, 94,
  // 			   57 , 51, 54, 59,
  // 			   96 , 97, 98, 99,
  // 			   153,148,154,157,
  // 			   1  ,  3,  6,  7,
  // 			   34 , 37, 39, 42,
  // 			   26 , 27, 28, 29,
  // 			   130,131,132,133};
  
  
  // int  pmtBList[nPMTsB] = {102,149,150,152,
  // 			   9  , 10, 12, 14,
  // 			   43 , 47, 48, 49,
  // 			   30 , 31, 32, 33,
  // 			   134,135,136,138};			   
  
  //-------------
  
  int  locAList[4] = {0,1,2,3};
  int  locBList[4] = {4,5,6,7};
  
  TString inputDirectory  = "/scratch/Gary/Data/";
  
  TString outputDirectory = "/scratch/Gary/Data/";
  
  int nEvents = -2;
  
  cout << endl;
  cout << " ----------------------- "  << endl;
  cout << " Processing Binary Data " << endl;
  
  for(int iRun = 0 ; iRun < nRuns ; iRun++ ){

    run = runList[iRun];
    
    // iPMT is index, ie not PMT serial number
    for (int iPMT = 0 ; iPMT < nPMTs ; iPMT++){
      
      if( iPMT < nPMTsA ){
	pmt = pmtAList[iPMT];
	loc = locAList[iPMT%4];
      }
      else{
	pmt = pmtBList[iPMT-nPMTsA];
	loc = locBList[(iPMT-nPMTsA)%4];
      }
      
      if(pmt == 0) 
	continue;

      for ( int iTest = 0 ; iTest < nTests ; iTest++ ){
	
	test = testList[iTest];
	
	// Gain test settings for
	// Further loop over HV steps
	if( test=='G'){
	  nSteps = 6;
	  hvStep = 1;
	}
	else{
	  nSteps = 1;
	  hvStep = 0;
	}
	
	for ( ; hvStep < nSteps ; hvStep++ )
	  nEvents = ProcessBinaryFile(inputDirectory,
				      outputDirectory, 
				      run, loc, 
				      pmt, hvStep,
				      test,
				      verbosity, 
				      digitiser);
	
	cout << endl;
	cout << " Output file contains " << nEvents << " events " << endl; 
	
      } // end:  for ( int iTest =...
    } // end: for (int iPMT = 0 ...
  } // end: for(int iRun = 0 ;...
  
  cout << endl;
  cout << " Completed Processing    " << endl;
  cout << " ----------------------- " << endl;
  
}
