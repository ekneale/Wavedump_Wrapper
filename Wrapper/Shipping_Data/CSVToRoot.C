#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"

void CSVToRoot(){
  
  TFile *outFile = new TFile("ShippingData_AllPMTs.root","RECREATE");
  
  TH1F *h1 = new TH1F("h1","x distribution",100,-4,4);
  
  TTree *T = new TTree("ShippingData","Arrg, this be Shipping Data from CSV file");
  Long64_t nlines;
  TString lineFormat = "PMT/I:Sk/F:Skb/F:Sp/F:Idb/F:EBB/I:DR/I:TTS/F:PTV/F";
  nlines = T->ReadFile("ShippingData_AllPMTs.csv",lineFormat,',');
  
  printf(" Shipping Data for %lld PMTs \n",nlines);

  Int_t cSize = 500;
  
  TCanvas * c1 =  new TCanvas("c1","c1",0,0,cSize,cSize);
  c1->Divide(2,2);

  c1->cd(1);
  T->Draw("EBB");

  c1->cd(2);
  T->Draw("DR");

  c1->cd(3);
  T->Draw("PTV");

  c1->cd(4);
  T->Draw("TTS");

  TCanvas * c2 =  new TCanvas("c2","c2",cSize,0,cSize,cSize);
  c2->Divide(2,2);

  c2->cd(1);
  T->Draw("Sk");

  c2->cd(2);
  T->Draw("Skb");

  c2->cd(3);
  T->Draw("Sp");

  c2->cd(4);
  T->Draw("Idb");
  
  T->Write();
}
