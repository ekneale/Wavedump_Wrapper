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
 *  Option 1: 
 *  $ ./BinToRoot 
 *  input data when prompted
 *
 * Dependencies
 *  root.cern
 *
 * Modified 
 *  gary.smith@ed.ac.uk
 *  15 02 19
 *
 *  Pulse/pulse -> Waveform/waveform - previous name was incorrect
 *  Note that this must be propagated to output file anaysis code.
 *
 *  removed unused variables qFixed and qPeak
 *
 *  New function ExecuteProcessing() acts as 
 *  intermediate step before ProcessBinaryFile()
 *  accommodating a cleaner main() function
 *
 *  New function ExecuteProcessingRun() acts
 *  as intermediate step before ExecuteProcessing()
 *  to accommodate processing entire runs
 *
 *  ProcessBinaryFile()
 *  added negPulsePol as optional argument 
 *  added samplingSetting  as optional argument 
 *
 *  Added new histograms to canvas 
 *  hPeakT_ns - time of pulse peak in ns
 *  hQFixed_PeakV - peak voltage (baseline subtracted)
 *  against charge accumulated around gate
 */ 

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TLine.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"

using namespace std;

int  GetNSamples(char digitiser,
		 char test){
  
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

int  GetNVDCBins(char digitiser){

  if     ( digitiser == 'V' )
    return 16384;
  else if( digitiser == 'D' )
    return 4096;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }
  
}

float GetVoltageRange(char digitiser){
  
  if     ( digitiser == 'V' )
    return 2.0;
  else if( digitiser == 'D' )
    return 1.0;
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }


}

int GetSampleRateInMHz(char digitiser,
		       char samplingSetting){

  if     ( digitiser == 'V' )
    return 500;
  else if( digitiser == 'D' ){
    if     (samplingSetting == 'S')
      return 5000;
    else if(samplingSetting == 'L') 
      return 1000;
    else{
      cerr << "Error: Unknown sampling setting " << endl;
      return 0;
    }
  }
  else{
    cerr << "Error: Unknown digitiser " << endl;
    return 0;
  }

}

void printDAQConfig(char digitiser,
		    char samplingSetting,
		    bool negPulsePol){
  

  TString strDigitiser = "VME";
  
  if (digitiser=='D')
    strDigitiser = "Desktop";

  cout << endl;
  cout << "   " 
       << strDigitiser
       << " Digitiser "    << endl;
  cout << "   " 
       << GetNVDCBins(digitiser)        
       << " VDC Bins "    << endl;
  cout << "   " 
       << GetVoltageRange(digitiser)    
       << " Volts Range "  << endl;
  cout << "   " 
       << GetSampleRateInMHz(digitiser, 
			     samplingSetting) 
       << " MHz sample rate " << endl;
  
  if ( negPulsePol )
    cout << "   Negative Pulse Polarity" 
	 << endl;
  else
    cout << "   Positive Pulse Polarity" 
	 << endl;
  
}

float GetnsPerSample(char digitiser,
		     char samplingSetting){
  
  return (1.0e3 / (float)GetSampleRateInMHz(digitiser,
					    samplingSetting));

}

float GetWaveformLength(char digitiser,
			char test,
			char samplingSetting){
  
  return (float)GetNSamples(digitiser,
			    test)*GetnsPerSample(digitiser,
						 samplingSetting);
  
}

float GetDelay(int run = 0){
  
  if      ( run == 0 ) // 12th October
    return 60.; 
  else if ( run < 10 ) // 16th October
    return 50.;
  else if ( run < 20 ) // Underground
    return 90.;
  else if ( run > 999) // Bine
    return 100.;
  else                 // Surface
    return 60.;
}

float GetmVPerBin(char digitiser){
  return ( 1.0e3 * GetVoltageRange(digitiser) / 
	   GetNVDCBins(digitiser) );
}

float GetCharge(int intVDC, char digitiser,
		char samplingSetting, bool negPulsePol){
  
  float nsPerSample = GetnsPerSample(digitiser,
				     samplingSetting);
  
  float mVolts = (float)intVDC * GetmVPerBin(digitiser);

  float charge = 0.;
  
  if(negPulsePol)
    charge = -1. * nsPerSample * mVolts;
  else
    charge = nsPerSample * mVolts;
  
  return charge; 
  
}

bool isCorrectDigitiser(int header,
			char digitiser,
			int test){
  
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

// Must leave room for baseline window
bool isPeakInRange(float peakT_ns, char digitiser,
		   char test, char samplingSetting){
  
  float waveformLength = GetWaveformLength(digitiser,test,samplingSetting);
  
  if( peakT_ns >= 50.   && 
      peakT_ns <  ( waveformLength - 30.0 ) )
    return true;
  else
    return false;
}

float gateWidth(){
  return 50;
}

// integration windows fixed wrt trigger
// baseline window before signal window only
int Accumulate_Fixed(short VDC, float time){

  // Integrate baseline using 
  // 50 ns window before signal
  // (max size given run 1 delay)
  // And signal in 50 ns window
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

int Accumulate_Baseline(short VDC, float time){

  // Integrate baseline using 
  // 50 ns window before signal
  if      ( time >= -gateWidth()  && 
	    time <    0 )
    return((int)VDC);
  else
    return 0;
}



// Integration windows wrt waveform peak
// -10 ns before to 30 ns after
// timeRel is time relative to minT (in ns)
// baseline window before
int Accumulate_Peak(short VDC, float timeRel){
  
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




TString GetRunFolderName(int run){

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

TString GetPMTFolderName(int pmt){

  if     ( pmt < 10  )
    return "PMT000%d/";
  else if( pmt < 100 )
    return "PMT00%d/";
  else if( pmt < 1000 )
    return "PMT0%d/";
  else
    return "PMT%d/";
}

TString GetTestFolderName(char test){
  
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


TString GetBinaryFilePath(TString filePath = "../../Data/",
			  int  run  = 1, 
			  int  pmt  = 1, 
			  int  loc  = 0,
			  char test = 'S',
			  int  hvStep = 0){
  
  TString rtnFilePath = "";
  
  // append file path
  filePath += GetRunFolderName(run);
  filePath += GetPMTFolderName(pmt);
  filePath += GetTestFolderName(test);
  
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

TString GetRawRootFilePath(TString filePath = "./",
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
    
    filePath += GetRunFolderName(run);
    filePath += GetPMTFolderName(pmt);
    filePath += GetTestFolderName(test);
    
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

void SetStyle(){
  
  TStyle     *binToRootStyle  = new TStyle("binToRootStyle",
					   "My Root Styles");
  
  const Int_t NCont = 255;
  const Int_t NRGBs = 5;
  
  // Color scheme for 2D plotting with a better defined scale 
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };          
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  
  binToRootStyle->SetNumberContours(NCont);
  
  // General
  // OPTIONS - FILL LINE TEXT MARKER
  
  binToRootStyle->SetFillColor(0);
  binToRootStyle->SetTextSize(0.05);
  
  //-----------  Canvas
  
  binToRootStyle->SetCanvasBorderMode(0);
  binToRootStyle->SetCanvasColor(kWhite);
  
  //------------- Pad
  
  binToRootStyle->SetPadBorderMode(0); 
  binToRootStyle->SetPadColor(kWhite);
  
  // Make more room for X and Y titles
  // one pad
  // binToRootStyle->SetPadRightMargin(0.05);  //percentage
  // binToRootStyle->SetPadLeftMargin(0.1);    //percentage
  // binToRootStyle->SetPadBottomMargin(0.12); //percentage
  
  // six sub-pads
  binToRootStyle->SetPadRightMargin(0.16);  //percentage
  binToRootStyle->SetPadLeftMargin(0.2);    //percentage
  binToRootStyle->SetPadBottomMargin(0.14); //percentage
  
  //----------- Histogram
  
  //Histos
  binToRootStyle->SetHistLineWidth(1);
  binToRootStyle->SetMarkerStyle(20);
  
  //  FILL CONTOURS LINE BAR 
  //  Frames
  binToRootStyle->SetFrameBorderMode(0);
  
  //  FILL BORDER LINE
  //  Graphs
  //  LINE ERRORS
  
  //---------  Axis 
  
  binToRootStyle->SetLabelFont(132,"XYZ"); 
  binToRootStyle->SetLabelSize(0.04,"XYZ");
  binToRootStyle->SetLabelOffset(0.01 ,"Y");
  
  //---------  Title
  binToRootStyle->SetOptTitle(1);
  binToRootStyle->SetTitleStyle(0);
  binToRootStyle->SetTitleBorderSize(0);


  binToRootStyle->SetTitleSize(0.03,"t");
  binToRootStyle->SetTitleFont(132,"t"); 

  binToRootStyle->SetTitleFont(132,"XYZ"); 

  binToRootStyle->SetTitleSize(0.05,"XYZ");
  
  binToRootStyle->SetTitleOffset(1.0,"XYZ");
  
  // 6 sub-pads
  binToRootStyle->SetTitleOffset(1.6,"Y");
  
  //----------  Stats
  binToRootStyle->SetOptStat(0);
  binToRootStyle->SetStatStyle(0);

  binToRootStyle->SetOptFit(1);
  
  //----------  Legend
  binToRootStyle->SetLegendBorderSize(0);
  //binToRootStyle->SetLegendFont(132);
  
  gROOT->SetStyle("binToRootStyle");
  gROOT->ForceStyle();
}

Bool_t IsCleanFFTWaveform(TH1D * hWave,
			  Char_t digitiser,
			  Char_t test){
  
  Bool_t waveformIsClean = kFALSE;
  
  TH1F* hWaveFFT  = new TH1F("hWaveFFT","hWaveFFT",
			     GetNSamples(digitiser,
					 test),
			     0,
			     GetNSamples(digitiser,
					 test));
  
  hWave->FFT(hWaveFFT ,"MAG");
  
  // delete zero frequency data
  hWaveFFT->SetBinContent(1,0.) ;
  
  if(hWaveFFT->GetMaximumBin() == 2)
    waveformIsClean = kTRUE;
  
  hWaveFFT->Delete();
  
  return waveformIsClean;
}

// Read in binary file and 
// write to root file
int ProcessBinaryFile(TString inFilePath,
		      TString outFilePath,
		      int run, int loc,
		      int pmt, int hvStep,
		      char test,
		      int  verbosity,
		      char digitiser,
		      bool negPulsePol = true,
		      char samplingSetting = 'S'
		      ){

  inFilePath = GetBinaryFilePath(inFilePath,
				 run, pmt, loc, test, hvStep);
  
  cout << endl;
  cout << " Processing  " << inFilePath << endl;

  //----------------------
  // Variables for testing
  bool  testMode  = false;
  bool  keepGoing = true;
  int   maxEvents = 100000;
  
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
  
  printDAQConfig(digitiser,
		 samplingSetting,
		 negPulsePol);
  
  outFilePath = GetRawRootFilePath(outFilePath,
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
  TString hQ_FilterNameTemp = "", hQ_FilterName = "";
  
  TString hQ_PeakNameTemp  = "", hQ_PeakName  = "";

  TString hPeakTimeNameTemp  = "", hPeakTimeName  = "";
  TString hPeakVoltNameTemp  = "", hPeakVoltName  = "";

  if(hvStep == 0){
    eventTreeNameTemp = "Events_Run_%d_PMT_%d_Loc_%d_Test_%c";
    eventTreeName.Form(eventTreeNameTemp,run,pmt,loc,test);

    canvasNameTemp = "Canvas_Run_%d_PMT_%d_Loc_%d_Test_%c";
    canvasName.Form(canvasNameTemp,run,pmt,loc,test);
    
    hQ_FixedNameTemp = "hQ_Fixed_Run_%d_PMT_%d_Loc_%d_Test_%c";
    hQ_FixedName.Form(hQ_FixedNameTemp,run,pmt,loc,test);

    hQ_FilterNameTemp = "hQ_Filter_Run_%d_PMT_%d_Loc_%d_Test_%c";
    hQ_FilterName.Form(hQ_FilterNameTemp,run,pmt,loc,test);
    
    hQ_PeakNameTemp = "hQ_Peak_Run_%d_PMT_%d_Loc_%d_Test_%c";
    hQ_PeakName.Form(hQ_PeakNameTemp,run,pmt,loc,test);

    hPeakTimeNameTemp = "hPeakTime_Run_%d_PMT_%d_Loc_%d_Test_%c";
    hPeakTimeName.Form(hPeakTimeNameTemp,run,pmt,loc,test);
    
    hPeakVoltNameTemp = "hPeakVolt_Run_%d_PMT_%d_Loc_%d_Test_%c";
    hPeakVoltName.Form(hPeakVoltNameTemp,run,pmt,loc,test);
    
  }
  else{
    eventTreeNameTemp = "Events_Run_%d_PMT_%d_Loc_%d_HV_%d";
    eventTreeName.Form(eventTreeNameTemp,run,pmt,loc,hvStep);

    canvasNameTemp = "Canvas_Run_%d_PMT_%d_Loc_%d_HV_%d";
    canvasName.Form(canvasNameTemp,run,pmt,loc,hvStep);

    hQ_FixedNameTemp = "hQ_Fixed_Run_%d_PMT_%d_Loc_%d_HV_%d";
    hQ_FixedName.Form(hQ_FixedNameTemp,run,pmt,loc,hvStep);

    hQ_FilterNameTemp = "hQ_Filter_Run_%d_PMT_%d_Loc_%d_HV_%d";
    hQ_FilterName.Form(hQ_FilterNameTemp,run,pmt,loc,hvStep);
    
    hQ_PeakNameTemp = "hQ_Peak_Run_%d_PMT_%d_Loc_%d_HV_%d";
    hQ_PeakName.Form(hQ_PeakNameTemp,run,pmt,loc,hvStep);

    hPeakTimeNameTemp = "hPeakTime_Run_%d_PMT_%d_Loc_%d_HV_%D";
    hPeakTimeName.Form(hPeakTimeNameTemp,run,pmt,loc,hvStep);

    hPeakVoltNameTemp = "hPeakVolt_Run_%d_PMT_%d_Loc_%d_HV_%D";
    hPeakVoltName.Form(hPeakVoltNameTemp,run,pmt,loc,hvStep);

  }
  //---------------------
  // Event Level Data
  TTree * eventTree  = new TTree(eventTreeName,
				 eventTreeName);
  
  UInt_t width = 1280, height = 800;
  
  SetStyle();
  
  TCanvas * canvas = new TCanvas(canvasName,
				 canvasName,
				 width,
				 height);
  
  Int_t nBinsX  = 512;
  Float_t rangeQ[2] = {-500.,1500.};
  
  TH1F * hQ_Fixed = new TH1F(hQ_FixedName,
			     "Gate around delay;Charge (mV ns);Counts",
			     nBinsX,rangeQ[0],rangeQ[1]);

  TString label;
  label = "Gate around delay filtered;Charge (mV ns);Counts";
  
  TH1F * hQ_Filter = new TH1F(hQ_FilterName,
			      label,
			      nBinsX,rangeQ[0],rangeQ[1]);
  
  TH1F * hQ_Peak = new TH1F(hQ_PeakName,
			    "Gate around peak;Charge (mV ns);Counts",
			    nBinsX,rangeQ[0],rangeQ[1]);
  
  Float_t rangeV[2] = {-5.,155.};
  
  TH1F * hPeakV_mV = new TH1F(hPeakVoltName,
			      "Waveform peak voltage; Voltage (mV);Counts",
			      nBinsX,
			      rangeV[0],rangeV[1]);
  
  TH1F * hPeakT_ns = new TH1F(hPeakTimeName,
			      "Time of waveform peak; Time (ns);Counts",
			      GetNSamples(digitiser,test),
			      0.,
			      GetWaveformLength(digitiser,
						test,
						samplingSetting));
  
  label =  "Signal charge vs peak voltage;";
  label += "Charge (mV ns);Peak voltage with baseline subtracted (mV)";
  
  TH2F * hQV = new TH2F("hQFixed_PeakV",
			label,
			nBinsX/4,rangeQ[0],rangeQ[1],
			nBinsX/4,rangeV[0],rangeV[1]);
  
				   
  TH2F * hWaveforms = new TH2F("hWaveforms",
			       "Subset of raw waveforms;Sample;VDC",
			       GetNSamples(digitiser,test),0,
			       (GetNSamples(digitiser,test)-1),
			       GetNVDCBins(digitiser),0.,
			       (GetNVDCBins(digitiser)-1));

  label = "Subset of calibrated waveforms;Time (ns);Voltage (mV)";
  
  rangeQ[0] = 0;
  rangeQ[1] = GetNSamples(digitiser,test)-1;
  rangeQ[1] = rangeQ[1] * GetnsPerSample(digitiser,
					 samplingSetting);
  
  TH2F * hTV = new TH2F("hTV",label,
			GetNSamples(digitiser,test),
			0.,GetWaveformLength(digitiser,
					     test,
					     samplingSetting),
			GetNVDCBins(digitiser),
			0.,GetVoltageRange(digitiser)*1.0e3);

  TH1D * hWave = new TH1D("hWave","hWaveform;Time (ns);ADC counts",
			  GetNSamples(digitiser,
				      test), 0., 
			  GetWaveformLength(digitiser,
					    test,
					    samplingSetting));
  
  int   event  = -1;
  
  // Note that the sign is preserved for VDC, therefore 
  // minima (maxima) will be where the signal peaked 
  // for negative (positive) signal pulses
  short minVDC = 32767, maxVDC = -32768;  
  // sample numbers (indices) corresponding to minVDC and maxVDC
  short minT   = -1, maxT   = -1;  
  
  float peakT_ns = -1., peakV_mV = 0., baselineV_mV =0.;
  
  // accumulators for integrating the sample values
  int   intVDCfixed = 0, intVDCpeak = 0, intVDCbaseline = 0;

    // read in samples per event
  // vector may be better
  short waveform[GetNSamples(digitiser,test)];
  
  eventTree->Branch("event",&event,"event/I");
  eventTree->Branch("minVDC",&minVDC,"minVDC/S");
  eventTree->Branch("maxVDC",&maxVDC,"maxVDC/S");
  eventTree->Branch("minT",&minT,"minT/S");
  eventTree->Branch("maxT",&maxT,"maxT/S");
  
  eventTree->Branch("peakT_ns",&peakT_ns,"peakT_ns/F");
  eventTree->Branch("peakV_mV",&peakV_mV,"peakV_mV/F");
  
  TString arrayString = "";
  arrayString.Form("waveform[%d]/S",GetNSamples(digitiser,test));
  
  eventTree->Branch("waveform",waveform,arrayString);
  
  //---------------------
  // Sample Level Data

  short VDC = 0, sample = 0;
  int   fileHeader = 0;
  float floatVDC   = 0.;
  
  // waveform time in ns
  float waveTime = 0.;
  
  // waveform time with delay subtracted
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
    intVDCfixed = 0, intVDCpeak = 0, intVDCbaseline = 0;
    
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
    
    // read in waveform which comes 
    // in 2 (VME) or 4 (Desktop) bit chunks
    for (short iSample = 0; iSample < GetNSamples(digitiser,test); iSample++){
  
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

      hWave->SetBinContent(iSample+1,
			   (float)(8700 - waveform[iSample]));
      
      
      // for writing
      sample = iSample;
      waveform[sample] = VDC;
      
      if     ( VDC < minVDC ){
	minVDC = VDC;
	minT   = sample;
      }
      else if( VDC > maxVDC ) {
	maxVDC = VDC;
	maxT   = sample;
      }

      waveTime = sample * GetnsPerSample(digitiser,
					 samplingSetting);

      
      // time with delay subtracted 
      // to allow common fixed time windows per run
      time = waveTime - GetDelay(run);
      
      // fixed window accumulations
      // add / subtract / skip sample VDC value
      intVDCfixed += Accumulate_Fixed(VDC,time);
      intVDCbaseline += Accumulate_Baseline(VDC,time);
      
      if(verbosity > 1){
	cout << " VDC(" << iSample << ") = " << VDC << endl;
	cout << " intVDCbaseline = " << intVDCbaseline << endl;
	cout << " time           = " << time << endl;
      }
    } // end: for (short iSample = 0; iSa
    
    



    // time of the signal peak
    
    float nSamplesPerGate = gateWidth() / 1000.;
    
    nSamplesPerGate = nSamplesPerGate * GetSampleRateInMHz(digitiser,
							   samplingSetting);
    
    baselineV_mV = (float)intVDCbaseline / nSamplesPerGate;
    
    if(negPulsePol){
      peakT_ns = minT * GetnsPerSample(digitiser,
				       samplingSetting);
      
      peakV_mV = baselineV_mV - waveform[minT];
      
      
    }
    else{
      peakT_ns = maxT * GetnsPerSample(digitiser,
				       samplingSetting);
      
      peakV_mV = waveform[maxT] - baselineV_mV;
    }

    peakV_mV = peakV_mV * GetmVPerBin(digitiser);

    float timeRelPeak = 0;
    
    // Loop over waveform again
    // necessary for when event-level variables are used
    // such as waveform minimum
    for (short iSample = 0; iSample < GetNSamples(digitiser,test); iSample++){
      
      waveTime = iSample*GetnsPerSample(digitiser,
					samplingSetting);
      
      time = waveTime - GetDelay(run);

      timeRelPeak = waveTime - peakT_ns;
	
      // accumulations wrt waveform peak (minimum)      
      if( isPeakInRange(peakT_ns,digitiser,test,samplingSetting) )
	intVDCpeak += Accumulate_Peak(waveform[iSample],
				      timeRelPeak);
      
      // sample vs VDC and time vs voltage plots
      // for checking signals (delay etc) 
      // plot pulses for 1000 events
      // from first 100,000
      if( ( event < 100000    )  &&
	  ( event % 100  == 0 )
	  ){
	
	hWaveforms->Fill(iSample,waveform[iSample]);
	
	hTV->Fill(time + GetDelay(run),
		  waveform[iSample]*GetmVPerBin(digitiser));
	
      }
      
    } // end: for (short iSample = 0; iS...
    
    hQ_Fixed->Fill(GetCharge(intVDCfixed,digitiser,
			     samplingSetting,negPulsePol));
    
    hPeakT_ns->Fill(peakT_ns);
    hPeakV_mV->Fill(peakV_mV);
    hQV->Fill(GetCharge(intVDCfixed,digitiser,
			samplingSetting,negPulsePol),
	      peakV_mV);
    
    if( isPeakInRange(peakT_ns,digitiser,test,samplingSetting) )
      hQ_Peak->Fill(GetCharge(intVDCpeak,digitiser,
			      samplingSetting,negPulsePol));
    
    if( IsCleanFFTWaveform(hWave,digitiser,test) )
      hQ_Filter->Fill(GetCharge(intVDCfixed,digitiser,
				samplingSetting,negPulsePol));

    //--------------------------------
    // Write event by event data here
    
    eventTree->Fill();
    
    if(verbosity > 0){
      cout << endl;
      cout << " Event " << event << endl;
      cout << " minVDC(" << minT << ") = " << minVDC << endl;
      cout << " maxVDC(" << maxT << ") = " << maxVDC << endl;
      cout << " peakV_mV               = " << peakV_mV << endl;
      cout << " peakT_ns               = " << peakT_ns << endl;
      cout << endl;
    }
    
    if( testMode &&
	(event+1) == maxEvents )
      keepGoing = false;
    
  } // end: while loop

  // close wavedump file
  fileStream.close();	

  canvas->Divide(3,2);

  canvas->cd(1);
  //=================================
  //   Calibrated Pulses
  float minY = GetVoltageRange(digitiser)*(16 - 2)/32*1.0e3;
  float maxY = GetVoltageRange(digitiser)*(16 + 1)/32*1.0e3 ;

  if(!negPulsePol){
    minY = GetVoltageRange(digitiser)*(16 - 1)/32*1.0e3;
    maxY = GetVoltageRange(digitiser)*(16 + 2)/32*1.0e3;
  }
  
  hTV->SetAxisRange(minY,maxY,"Y");

  float minX = 0.;
  float maxX = 500.; // 220.
  
  if(test=='A'){
    hTV->SetAxisRange(minX,maxX,"X");
  }
  
  hTV->Draw("colz");
  
  float lineYMin = minY * (16 - 1)/(16 - 2); 
  float lineYMax = maxY * (16 + 0.25)/(16 + 1); 
  
  float lineXMin = GetDelay(run) - gateWidth();
  float lineXMax = GetDelay(run) + gateWidth();
  
  TLine *lPedMin = new TLine(lineXMin,lineYMin,
			     lineXMin,lineYMax); 
  lPedMin->SetLineColor(kRed);
  
  TLine *lSigMin = new TLine(GetDelay(run),lineYMin,
			     GetDelay(run),lineYMax); 
  lSigMin->SetLineColor(kBlue);
  
  TLine *lSigMax = new TLine(lineXMax,lineYMin,
			     lineXMax,lineYMax); 
  lSigMax->SetLineColor(kBlue);
  
  if(negPulsePol){
    lPedMin->Draw("same");
    lSigMin->Draw("same");
    lSigMax->Draw("same");
  }
  
  //=================================
  canvas->cd(2);
  gPad->SetTicks();
  //=================================
  //  Gate around Delay
  hQ_Fixed->SetAxisRange(-500., 2500.,"X");
  hQ_Filter->SetLineColor(kRed);

  gPad->SetLogy(1);
  hQ_Fixed->Draw();
  hQ_Filter->Draw("same");
  
  //=================================
  canvas->cd(3);
  gPad->SetTicks();
  //=================================
  //  Gate around Peak
  gPad->SetLogy(1);
  
  hQ_Peak->SetAxisRange(-500., 2500.,"X");
  hQ_Peak->Draw();
  
  //=================================
  canvas->cd(4);
  gPad->SetTicks();
  //=================================
  //  Time of Peak
  
  hPeakT_ns->Draw();
  
  //=================================
  canvas->cd(5);
  gPad->SetTicks();
  //=================================
  // Charge vs Peak
  gPad->SetLogz(1);

  hQV->Draw("colz");

  //=================================
  canvas->cd(6);
  //=================================
  //  Raw Waveforms
  
  minY = GetNVDCBins(digitiser)*10/16 ;
  maxY = GetNVDCBins(digitiser)*6/16 ;
  
  hWaveforms->SetAxisRange(minY,maxY,"Y");
  
  hWaveforms->Draw("colz");
  //=================================
  
  canvasName = "./Plots/";
  canvasName += canvas->GetName();
  canvasName += ".png";
  canvas->SaveAs(canvasName);
  
  //--------------------------------
  // Write file info here

  hWaveforms->Delete();
  hTV->Delete();
  hQV->Delete();
  hWave->Delete();
  
  eventTree->Write();
  eventTree->Delete();

  outFile->Write();
  outFile->Close();
  
  return (event+1);
}


string GetDefaultFilePath(){
  return  "../../Data/wave_0.dat";
}

bool GetNegPulsePolUser(){
  
  bool negPulsePol = true;
  
  cout << " Enter pulse polarity: ";
  cout << endl;
  cout << " 1  - negative  "; 
  cout << endl;
  cout << " 0  - positive  ";
  cin >> negPulsePol;
  cout << endl;
  
  return negPulsePol;
}

bool GetNegPulsePol(char digitiser, int run){
  
  if(digitiser == 'V')
    return true;
  else if( run > 9999)
    return false;
  else 
    return GetNegPulsePolUser();
      
 
}

char GetSamplingSettingUser(){
  
  char samplingSetting = 'S';

  cout << " Enter sampling frequency: " << endl;
  cout << " 'S' (Standard)  5 GHZ "    << endl;
  cout << " 'L' (Low)       1 GHZ "    << endl;
  cin  >> samplingSetting;
  cout << endl;
  
  return samplingSetting;
}  

char GetSamplingSetting(char digitiser,
			int run){

  if  (digitiser == 'V')
    return 'S'; 
  if  (run > 9999)
    return 'L';
  else 
    return GetSamplingSettingUser();
}

int GetRunUser(){
  int run;
  cout << " Enter run number: ";
  cin >> run;
  cout << endl;
  return run;
}

int GetPMTUser(){
  int pmt;
  cout << " Enter pmt number: ";
  cin >> pmt;
  cout << endl;
  return pmt;
}

int GetLocUser(){
  int loc;
  cout << " Enter rig location: ";
  cin >> loc;
  cout << endl;
  return loc;
}

char GetTestUser(){
  char test;
  cout << " Enter test type : ";
  cin >> test;
  cout << endl;
  return test;
}

TString GetInDirUser(){
  char * inDir = new char[128];
  cout << " Enter input path: ";
  cin >> inDir;
  cout << endl;
  return inDir;
}

TString GetInDir(int run){

  if( run < 30 ){
    return "/Disk/ds-sopa-group/PPE/Watchman/BinaryData/";
  }
  else if( run == 1000)
    return "/Volumes/GS_External/BinaryData/";
  else if( run >= 1001)
    return "/Users/gsmith23/Desktop/Watchman/Testing/Wavedump_Wrapper/BinaryData/";
  else 
    return GetInDirUser();
}


TString GetOutDirUser(){
  char * outDir = new char[128];
  cout << " Enter output path: ";
  cin >> outDir;
  cout << endl;
  return outDir;
}


TString GetOutDir(int run){

  if( run < 30 ){
    return "/Disk/ds-sopa-group/PPE/Watchman/RawRootData/";
  }
  else if( run == 1000)
    return "/Volumes/GS_External/RawRootData/";
  if( run == 1001)
    return "/Users/gsmith23/Desktop/Watchman/Testing/Wavedump_Wrapper/RawRootData/";
  else 
    return GetOutDirUser();
}

char GetDigitiserUser(){
  char digitiser;
  cout << " Enter digitiser type ('V', or 'D'): ";
  cin >> digitiser;
  cout << endl;
  return digitiser;
}

char GetDigitiser(int run){

  if   (run >= 1000)
    return 'D';
  else 
    return 'V';
}


int GetHVStep(char test ){
  if( test == 'G')
    return 1;
  else
    return 0;
}

int GetNSteps(char test ){
  if( test == 'G')
    return 5;
  else
    return 0;
}

void ExecuteProcessing(int run = 0, int pmt = 0,
		       int loc = 0 ,char test = 'A',
		       TString inDir = "./", 
		       TString outDir = "./",
		       char digitiser = 'V',
		       char verbosity = 0
		       ){

  bool negPulsePol = GetNegPulsePol(digitiser,run);
  char samplingSetting = GetSamplingSetting(digitiser,run);


  int  nEvents = -1;

  for (int hvStep = GetHVStep(test) ; 
       hvStep <= GetNSteps(test)    ; 
       hvStep++ ){
    
    nEvents = ProcessBinaryFile(inDir,
				outDir, 
				run, loc, 
				pmt, hvStep,
				test,
				verbosity, 
				digitiser,
				negPulsePol,
				samplingSetting);
    
    cout << endl;
    cout << " Output file contains " << nEvents << " events " << endl; 
  }  
  
  cout << endl;
  cout << " Completed Processing    " << endl;
  cout << " ----------------------- " << endl;
}

int GetPMT(int run, int iPMT){
  
  if(run == 1){
    int  pmtList[100] = {83 , 88,108,107, // Tent A
			 73 , 76, 84, 87,
			 66 , 78, 82,103,
			 104,106,112,141,
			 61 , 65, 75,105,
			 74 ,111,140,142,
			 143,145,146,147,
			 63 , 67,158,160,
			 139,161,164,165,
			 90 ,159,166,171,
			 81 ,167,169,170,
			 50 , 53,162,163,
			 55 , 56, 92, 94,
			 57 , 51, 54, 59,
			 96 , 97, 98, 99,
			 153,148,154,157,
			 1  ,  3,  6,  7,
			 34 , 37, 39, 42,
			 26 , 27, 28, 29,
			 130,131,132,133, // Tent A
			 102,149,150,152, // Tent B
			 9  , 10, 12, 14,
			 43 , 47, 48, 49,
			 30 , 31, 32, 33,
			 134,135,136,138};// Tent B	   
    
    return pmtList[iPMT];
  }
  else if(run == 1001){
    return 16;
  }
  else
    return 1;
}

int GetLoc(int run, int iPMT){
  
  if(run == 1)
    if( iPMT < 80 )
      return (iPMT%4);
    else
      return (iPMT%4 + 4);
  else
    return 0;

}

int GetNPMTs(int run){

  if(run == 1)
    return 100;
  else
    return 1;
}

char GetTest(char choice, int iTest){

  char testList[5] = {'S','N','G','D','A'};
  
  if( choice != 'E' )
    return choice;
  else
    return testList[iTest]; 

}

int GetNTests(char choice){
  if( choice != 'E' )
    return 1;
  else
    return 5;
}


char GetChoiceUser(){

  char choice;
  cout << " Which test/s to process? " << endl;
  cout << " 'S' - SPEtest "           << endl;
  cout << " 'N' - Nominal "           << endl;
  cout << " 'A' - APTest "            << endl;
  cout << " 'D' - DarkRateTest "      << endl;
  cout << " 'G' - GainTest     "      << endl;
  cout << " 'E' - Every Test     "    << endl;
  cin  >> choice;

  return choice;
}

char GetChoice(int run){

  if( run > 100 ){
    return 'S';
  }
  else
    return GetChoiceUser();
}



void ExecuteProcessingRun(int run){
    
  char choice = GetChoice(run);

  for ( int iPMT = 0 ; iPMT < GetNPMTs(run) ; iPMT++)
    for( int iTest = 0 ; iTest < GetNTests(choice) ; iTest++)
      ExecuteProcessing(run,
			GetPMT(run,iPMT),
			GetLoc(run,iPMT),
			GetTest(choice,iTest),
			GetInDir(run),
			GetOutDir(run),
			GetDigitiser(run)
			);
  
}

int main(int argc, char **argv)
{
  
  cout << endl;
  cout << " BinToRoot " << endl;
  cout << endl;

  int run;
  int pmt;
  int loc;
  char test;
  TString inDir;
  TString outDir;
  char digitiser;

  switch ( argc ){
  case ( 1 ):

    run  = GetRunUser();
    pmt  = GetPMTUser();
    loc  = GetLocUser();
    test = GetTestUser();
    inDir  = GetInDirUser();
    outDir = GetOutDirUser();
    digitiser =  GetDigitiserUser();
    
    ExecuteProcessing(run,
    		      pmt,
    		      loc,
    		      test,
    		      inDir,
    		      outDir,
    		      digitiser);
    break;
  case ( 2 ):
    char * endPtr;
    
    ExecuteProcessingRun(strtol(argv[1], &endPtr, 10));
    
    break;
  case ( 8 ):

    ExecuteProcessing(strtol(argv[1], &endPtr, 10),
		      strtol(argv[2], &endPtr, 10),
		      strtol(argv[3], &endPtr, 10),
		      *argv[4],
		      argv[5],
		      argv[6],
		      *argv[7]);
    
    break;
  default:
    cerr << " Invalid option " << endl;
    return 0;
  }
    
  //
  //---------------
    
  return 1;  
}
