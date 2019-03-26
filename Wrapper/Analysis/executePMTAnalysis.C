/*****************************************************
 * A program to process raw root files from wavedump
 *
 * Author 
 *  gary.smith@ed.ac.uk
 *  01 01 2019
 *
 * Purpose
 *  This program reads a TTree
 *
 * How to build
 *  $ make 
 *
 * How to run
 *
 *  Example One - analyse one raw root file
 *  $ ./executePMTAnalysis ../../RawRootData/Run_11_PMT_90_Loc_0_Test_D.root  
 *
 *  Example Two - analyse two raw root files (two HV steps)
 *  $ ./executePMTAnalysis ../../RawRootData/Run_11_PMT_90_Loc_0_HV_1.root ../../RawRootData/Run_11_PMT_90_Loc_0_HV_2.root 
 *
 *  Example Three - analyse multiple raw root files ( all run numbers )
 *  $ ./executePMTAnalysis ../../RawRootData/Run_*_PMT_90_Loc_0_Test_D.root  
 *
 * Dependencies
 *  root.cern
 *  PMTAnalyser
 *  DataInfo
 *  ShippingData
 *
 * Known issues
 *   Under linux you will likely be required to add 
 *   the shared (.so) file to your library path 
 *
 */ 

#include <iostream>

#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "PMTAnalyser.h"
#include "ShippingData.h"
#include "fileNameParser.h"
    
#include "TCanvas.h"

using namespace std;

int main(Int_t argc, Char_t *argv[]){
  
  TFile * file = nullptr;
  TTree * tree = nullptr;
  
  PMTAnalyser * PMT = nullptr;

  ShippingData * shipData = nullptr;

  Char_t  digitiser = 'V';

  // Old style BinToRoot output
  // or new BinToRoot output?
  // (pulse[] -> waveform[] e.g.)
  Bool_t oldRootFileVersion = kFALSE;

  // Dark Rate
  Float_t thresh_mV  = 10.0;
  Int_t   darkRate   = 8000.;
  
  // argv should be a path to a file
  // or list of files ( wildcards work )
  for( int iFile = 1 ; iFile < argc ; iFile++){
     
     file = new TFile(argv[iFile],"READ");
     
     if (!file || !file->IsOpen()) {
       file = new TFile();
       cerr << " Error, Check File: " << argv[iFile] << endl; 
       return -1;
     }

     //-------
     // use parser to extract test ID
     cout << endl;
     cout << " Run      " << run(argv[iFile])      << endl;
     cout << " PMT      " << pmtID(argv[iFile])    << endl;
     cout << " Location " << location(argv[iFile]) << endl;
     cout << " Test     " << test(argv[iFile]) << endl;
     
     if( test(argv[iFile])=='G')
       cout << " HV step  " << HVStep(argv[iFile]) << endl; 
     
     cout << endl;
     shipData = new ShippingData(pmtID(argv[iFile]));

     // connect to tree in input file
     file->GetObject((TString)GetTreeName(argv[iFile]),tree); 
     
     // initalise analysis object using tree 
     PMT = new PMTAnalyser(tree,
			   digitiser,
			   test(argv[iFile]),
			   oldRootFileVersion);
     
     // Set plot attributes to bespoke TStyle 
     PMT->SetStyle();
     
     // Limit to subset of entries for quicker testing
     //PMT->SetTestMode(kTRUE);
     
     cout << " Hamamatsu Dark Rate = " << shipData->GetDR() << endl;

     Bool_t investigateDarkRate = kFALSE;

     if(investigateDarkRate){
       if(oldRootFileVersion){
	 cout << " Dark Rate method only applicable to new BinToRoot files " << endl;
       }
       else{
	 darkRate   = PMT->DarkRate(thresh_mV);
	 cout << " PMT Test  Dark Rate = " << darkRate          << endl;
       }
     }

     Bool_t investigateFFT = kTRUE;
     // Make Filtered Histograms
     if(investigateFFT){ 
       TCanvas * canvas = PMT->Make_FFT_Canvas();
       
       TString canvasName;
       canvasName.Form("./Plots/FFT_Run_%d_PMT_%d_Test_%c.pdf",
		       run(argv[iFile]),
		       pmtID(argv[iFile]),
		       test(argv[iFile]));
       
       if(canvas){
	 canvas->SaveAs(canvasName);
	 
	 canvasName.Form("./Plots/FFT_Run_%d_PMT_%d_Test_%c.root",
			 run(argv[iFile]),
			 pmtID(argv[iFile]),
			 test(argv[iFile]));
	 
	 canvas->SaveAs(canvasName);
       }
       else{
	 cout << endl;
	 cout << " No canvas produced " << endl;
       }
     }
       
  }
  
  // To Do:
  // TFile * resultsFile = new TFile("Results","RECREATE");
  
  return 1;
}
