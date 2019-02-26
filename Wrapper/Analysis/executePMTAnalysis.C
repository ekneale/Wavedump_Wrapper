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

using namespace std;

int main(Int_t argc, Char_t *argv[]){
  
  TFile * file = nullptr;
  TTree * tree = nullptr;
  
  PMTAnalyser * PMT = nullptr;

  ShippingData * shipData = nullptr;

  Char_t  digitiser = 'V';
  
  //Float_t thresh_mV  = 10.0;
  //Int_t   darkRate   = 8000.;
  
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
			   test(argv[iFile]));
     
     // Double Chooz dark rate threshold (x10 to include amplification)
     // darkRate   = PMT->DarkRate(thresh_mV);
     // cout << " Hamamatsu Dark Rate = " << shipData->GetDR() << endl;
     // cout << " PMT Test  Dark Rate = " << darkRate          << endl;
     
     // Make Filtered Histograms
     PMT->Make_FFT_Histos();
     
     //PMT->Make_Fixed_Gate_Filtered();

  }
  
  // To Do:
  // TFile * resultsFile = new TFile("Results","RECREATE");
  
  return 1;
}
