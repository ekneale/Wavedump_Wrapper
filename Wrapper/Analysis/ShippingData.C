#define ShippingData_cxx
#include "ShippingData.h"

#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void ShippingData::SetAll(int userPMT)
{
   if (fTree == 0) return;

   Long64_t nentries = fTree->GetEntriesFast();

   for (Long64_t jentry = 0; jentry < nentries ; jentry++ ) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     
     fTree->GetEntry(jentry);

     if(userPMT == PMT){
       //Show();
       userSk  = Sk;
       userSkb = Skb;
       userSp  = Sp;
       userIdb = Idb;
       userEBB = EBB;
       userDR  = DR;
       userTTS = TTS;
       userPTV = PTV;
       return;
     }
     
     if( PMT > userPMT || 
	 ientry == (nentries-1)){
       std::cerr << " No Shipping Data for PMT " << userPMT << std::endl;
       break;
     }
   
   }
}

void ShippingData::SetNewPMT(int userPMT){
  SetAll(userPMT);
}

Float_t ShippingData::GetSk(){
  return userSk;
}
Float_t ShippingData::GetSkb(){
  return userSkb;
}
Float_t ShippingData::GetSp(){
  return userSp;
}
Float_t ShippingData::GetIdb(){
  return userIdb;
}
Int_t   ShippingData::GetEBB(){
  return userEBB;
}
Int_t   ShippingData::GetDR(){
  return userDR;
}
Float_t ShippingData::GetTTS(){
  return userTTS;
}
Float_t ShippingData::GetPTV(){
  return userPTV;
}
