#ifndef ShippingData_h
#define ShippingData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ShippingData {
public :
   TTree          *fTree;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           PMT;
   Float_t         Sk;
   Float_t         Skb;
   Float_t         Sp;
   Float_t         Idb;
   Int_t           EBB;
   Int_t           DR;
   Float_t         TTS;
   Float_t         PTV;

   // List of branches
   TBranch        *b_PMT;   //!
   TBranch        *b_Sk;   //!
   TBranch        *b_Skb;   //!
   TBranch        *b_Sp;   //!
   TBranch        *b_Idb;   //!
   TBranch        *b_EBB;   //!
   TBranch        *b_DR;   //!
   TBranch        *b_TTS;   //!
   TBranch        *b_PTV;   //!

   ShippingData(int PMT = 0);
   virtual ~ShippingData();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Show(Long64_t entry = -1);
   
   void    SetNewPMT(int userPMT);
   
   Float_t GetSk();
   Float_t GetSkb();
   Float_t GetSp();
   Float_t GetIdb();
   Int_t   GetEBB();
   Int_t   GetDR();
   Float_t GetTTS();
   Float_t GetPTV();
   
 private: 

   Float_t userSk;
   Float_t userSkb;
   Float_t userSp;
   Float_t userIdb;
   Int_t   userEBB;
   Int_t   userDR;
   Float_t userTTS;
   Float_t userPTV;

   void  SetAll(int userPMT);

};

#endif

#ifdef ShippingData_cxx
ShippingData::ShippingData(int PMT)
{
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ShippingData.root");
  if (!f || !f->IsOpen()) {
    f = new TFile("ShippingData.root");
  }

  TTree * tree = nullptr;
  f->GetObject("ShippingData",tree);
  Init(tree);
  SetAll(PMT);
  
}

ShippingData::~ShippingData()
{
   if (!fTree) return;
   delete fTree->GetCurrentFile();
}

Int_t ShippingData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fTree) return 0;
   return fTree->GetEntry(entry);
}
Long64_t ShippingData::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fTree) return -5;
   Long64_t centry = fTree->LoadTree(entry);
   if (centry < 0) return centry;
   if (fTree->GetTreeNumber() != fCurrent) {
      fCurrent = fTree->GetTreeNumber();
   }
   return centry;
}

void ShippingData::Init(TTree *tree)
{
   if (!tree) return;
   fTree = tree;
   fCurrent = -1;
   fTree->SetMakeClass(1);

   fTree->SetBranchAddress("PMT", &PMT, &b_PMT);
   fTree->SetBranchAddress("Sk", &Sk, &b_Sk);
   fTree->SetBranchAddress("Skb", &Skb, &b_Skb);
   fTree->SetBranchAddress("Sp", &Sp, &b_Sp);
   fTree->SetBranchAddress("Idb", &Idb, &b_Idb);
   fTree->SetBranchAddress("EBB", &EBB, &b_EBB);
   fTree->SetBranchAddress("DR", &DR, &b_DR);
   fTree->SetBranchAddress("TTS", &TTS, &b_TTS);
   fTree->SetBranchAddress("PTV", &PTV, &b_PTV);
}

void ShippingData::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fTree) return;
   fTree->Show(entry);
}

#endif // #ifdef ShippingData_cxx
