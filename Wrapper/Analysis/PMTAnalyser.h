#ifndef PMTAnalyser_h
#define PMTAnalyser_h

#include <iostream>

#include "TTree.h"
#include "TFile.h"

#include "TH2.h"

#include "DataInfo.h"

using namespace std;

class PMTAnalyser {
 public :
  TTree          *rawRootTree;  
  Int_t           fCurrent;

  char testType;
  char digitiser;
  
  DataInfo      *dataInfo;
  Int_t          NSamples;
  Int_t          NVDCBins;
  Float_t        VoltageRange;
  Float_t        mVPerBin;
  Float_t        nsPerSample;

  Int_t          event;
  
  Short_t        minVDC;
  Short_t        maxVDC;
  Short_t        minT;
  Short_t        maxT;
  
  Float_t        peakT_ns;
  Float_t        peakV_mV;
  
  Short_t        waveform[5100];
  
  TBranch        *b_event;
  TBranch        *b_minVDC;
  TBranch        *b_maxVDC;
  TBranch        *b_minT;  
  TBranch        *b_maxT;  
  TBranch        *b_peakT_ns;  
  TBranch        *b_peakV_mV;  
  TBranch        *b_waveform; 
  
  PMTAnalyser(TTree *tree=0, Char_t digitiser='V', Char_t userTest='S');
  virtual ~PMTAnalyser();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree,Char_t digitiser,Char_t userTest);
  virtual Int_t    DarkRate(Float_t);
  virtual Int_t    Make_FFT_Histos();
  virtual Int_t    Make_hFixed_Filtered();
  virtual Bool_t   IsCleanFFTWaveform(TH1D *);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PMTAnalyser_cxx

PMTAnalyser::PMTAnalyser(TTree *tree,
			 Char_t digitiser,
			 Char_t userTest) : rawRootTree(0) 
{
  Init(tree,digitiser,userTest);
}

PMTAnalyser::~PMTAnalyser()
{
  if (!rawRootTree) return;
  delete rawRootTree->GetCurrentFile();
}

Int_t PMTAnalyser::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!rawRootTree) return 0;
   return rawRootTree->GetEntry(entry);
}
Long64_t PMTAnalyser::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!rawRootTree) return -5;
   Long64_t centry = rawRootTree->LoadTree(entry);
   if (centry < 0) return centry;
   if (rawRootTree->GetTreeNumber() != fCurrent) {
      fCurrent = rawRootTree->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PMTAnalyser::Init(TTree *tree,
		       Char_t digitiser,
		       Char_t userTest)
{
  
  testType = userTest;

  dataInfo = new DataInfo();

  NSamples     = dataInfo->GetNSamples(testType,
				       digitiser);
  VoltageRange = dataInfo->GetVoltageRange(digitiser);
  NVDCBins     = dataInfo->GetNVDCBins(digitiser);
  mVPerBin     = dataInfo->GetmVPerBin(digitiser);  
  nsPerSample  = dataInfo->GetnsPerSample(digitiser);  
  
  if (!tree) return;
  rawRootTree = tree;
  fCurrent = -1;
  rawRootTree->SetMakeClass(1);
  
  rawRootTree->SetBranchAddress("event", &event, &b_event);
  rawRootTree->SetBranchAddress("minVDC", &minVDC, &b_minVDC);
  rawRootTree->SetBranchAddress("maxVDC", &maxVDC, &b_maxVDC);
  rawRootTree->SetBranchAddress("minT", &minT, &b_minT);
  rawRootTree->SetBranchAddress("maxT", &maxT, &b_maxT);

  //!!!!!!!!
  // temporary
  Bool_t newFile = kFALSE;
  
  if(newFile){
    rawRootTree->SetBranchAddress("waveform", waveform, &b_waveform);
    rawRootTree->SetBranchAddress("peakT_ns", &peakT_ns, &b_peakT_ns);
    rawRootTree->SetBranchAddress("peakV_mV", &peakV_mV, &b_peakV_mV);
  }
  else
    rawRootTree->SetBranchAddress("pulse",waveform, &b_waveform);

  Notify();
}

Bool_t PMTAnalyser::Notify()
{
   return kTRUE;
}

void PMTAnalyser::Show(Long64_t entry)
{
   if (!rawRootTree) return;
   rawRootTree->Show(entry);
}

#endif // #ifdef PMTAnalyser_cxx
