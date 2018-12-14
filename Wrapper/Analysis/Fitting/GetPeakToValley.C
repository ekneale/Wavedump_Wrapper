#include<iostream>
#include<fstream>
#include<vector>
#include <string>

#include"TMath.h"

#include"TFile.h"
#include"TTree.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "TString.h"

#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDerivative.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "RooPolynomial.h"
#include "RooBifurGauss.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooExponential.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

using namespace std;

// Result

typedef struct {

double value;
double error;


} ValueWithError;

typedef struct {

  ValueWithError peak;
  ValueWithError peakToValley;
  int id;
  
} Result;

void fillValueWithError(ValueWithError* val,TH1F* h){
  val->value = h->GetMean();
  val->error = h->GetRMS();
}

void fillValueWithError(ValueWithError* val,RooRealVar* var){
  val->value = var->getVal();
  val->error = var->getError();
}


#include <tuple>
typedef std::tuple<double,double,double,double,double,double> InitParams;

using namespace RooFit;

InitParams initializeFit(TH1F* h){

  TSpectrum *s = new TSpectrum(2,3);
  Int_t nfound = s->Search(h,2,"goff",0.0002);
  std::cout << "found peaks " << nfound << std::endl;

 
  
  /*** returns positions of Pedestal and Signal approx.****/
  Double_t *peaks;
  peaks = s->GetPositionX();
   std::cout << peaks[0] << " " << peaks[1]  << std::endl;

   if (peaks[1] < 100) peaks[1] = 400; 
   
  /*** Find the valley ***/
  Int_t SigSig = h->FindBin(peaks[1]); //Signal bin approx **used to find valley position
  Int_t PedPed  = h->GetMaximumBin(); //pedestal bin
   int PedMax = h->GetMaximum(); // pedestal events


   
  int Compare[2] = {PedMax,PedMax}; // array to store 2 numbers for comparison to determine which is larger
  Int_t ValleyBin = 0; //to hold valley bin number


  /****Finding the Valley*****/
  for(int i = PedPed; i< SigSig ;i++){
    Compare[0] = h->GetBinContent(i);

    if(Compare[0] < Compare[1]){
      Compare[1] = Compare[0]; //finds minimum number of events
      ValleyBin = i; //sets bin number for minimum
    }//end if
  }//end for

  //  if (peaks[1] < 500) peaks[1] = 400;
  
  double ssignal = (peaks[1] -  h->GetBinCenter(ValleyBin))/3.;
  double sped = (h->GetBinCenter(ValleyBin) - peaks[0])/10.;
  
  // hack
  //  if (sped > 30) sped = 30;
  
   /*****Finding Pedestal events****/
  int Noise = 0; //to hold pedestal events
  for(int i = 0;i< ValleyBin;i++){
    Noise += h->GetBinContent(i); // total noise events
  }//end for

  /** Approximate fraction of photons per event**/
  int Events = h->GetEntries();  // # of events
  double Ratio = (double) Noise / (double) Events; // Ratio of Random events to photon events

  std::cout << "Peaks " << peaks[0] << " " << peaks[1] << " " << sped <<   " " << ssignal << " "   <<  Ratio << std::endl;

  double valleyPos = h->GetBinCenter(ValleyBin) ;

 
  //return InitParams(5,45, 2,5,0.4);
  return InitParams(peaks[0],peaks[1]- peaks[0],sped,ssignal,Ratio,valleyPos);
}

//  std::get<0>(t) 



TH1F* h2h(TH1D* hold ){

  TH1F* h = new TH1F("histo", "histo", hold->GetNbinsX(),hold->GetXaxis()->GetXmin(),  hold->GetXaxis()->GetXmax());
  h->Sumw2();
  double sum =0;
  for (auto i = 1 ; i <= h->GetNbinsX(); ++i ){
    sum += hold->GetBinContent(i);
  }
  int entries = hold->GetEntries();
  for (auto i = 1 ; i <= h->GetNbinsX(); ++i ){
    // sum += hold->GetBinContent(i);
    int bincontent = TMath::Nint(entries*hold->GetBinContent(i)/sum);
    double valx = hold->GetBinCenter(i);
    for (int i = 0; i < bincontent; ++i){
      h->Fill(valx);
    }
  
  }
  std::cout << "histo " << h->GetEntries()<< " " <<  hold->GetEntries() <<  " " << sum <<  std::endl;
  return h;
}



TH1D* findHisto(TFile *input, int histoVersion = 0){

  if (histoVersion == 0 ){
    TH1D* histo = (TH1D*)input->Get("SPE0");
    if (histo) return histo;
    histo = (TH1D*)input->Get("SPE1");
    if (histo) return histo;
    histo = (TH1D*)input->Get("SPE2");
    if (histo) return histo;
    histo = (TH1D*)input->Get("SPE3");
    if (histo) return histo;
  }
  else{
   
    TString hName = "";
    hName.Form("hCharge%d",histoVersion);
    cout << " hName " << hName << endl;
    
    TH1D* histo = (TH1D*)input->Get(hName);

    return histo;
  }
  std::cout << "no histo " << std::endl;

  return 0;
}

Result* fitModel(string file = "PMT_NB0066_HV1820_Analysis.root", 
		 int histoVersion = 0,
		 double minval = -500,
		 double maxval = 2000 ){
  
  
  
  Result* res = new Result();
  
  TCanvas * canvas = new TCanvas("Canvas", "Canvas");
  
  //gROOT -> ProcessLine( ".x ./mattStyle.C" );
  
  TFile *input=new TFile(file.c_str());
  TH1D* histo = findHisto(input,histoVersion);
  if (histo ==0) return 0;
  
  TH1F* fhisto = h2h(histo);
  fhisto->GetXaxis()->SetTitle("charge [mV ns]");
  fhisto->GetYaxis()->SetTitle("Counts");
  fhisto->GetYaxis()->SetTitleOffset(1.3);
  fhisto->GetYaxis()->SetTitleFont(132);
  
  // intial guesses with TSpectrum
  InitParams params = initializeFit(fhisto);
  double valleyPos = std::get<5>(params);
  double peakPos = std::get<0>(params)+ std::get<1>(params);

  // gaussian fit to max
  fhisto->Fit("gaus", "","", peakPos - 100, peakPos +120 );
  TF1* gf = fhisto->GetFunction("gaus");
  TF1* f = (TF1*)gf->Clone();
  double sval = f->Eval(f->GetParameter(1));
  double esval = TMath::Abs(f->Eval(gf->GetParError(1) + f->GetParameter(1) ) - sval) ;
  f->SetLineColor(2);
  f->SetLineWidth(3);
  f->Draw("SAME");


  TFitResultPtr pres= fhisto->Fit("pol2", "S", "",valleyPos - 50 , valleyPos +50);
  TFitResult* polres = pres.Get();
  double ea = polres->Error(2); double eb = polres->Error(1);

  // get correlation
  TMatrixDSym cm = polres->GetCorrelationMatrix();
  TMatrixDSym cov = polres->GetCovarianceMatrix();
  double rho = cov[1][2];

 
  // find min of parabola
  TF1* f2 = fhisto->GetFunction("pol2");
  double a = f2->GetParameter(2); double b = f2->GetParameter(1);
  double vval = f2->GetMinimum( valleyPos -50, valleyPos +50);
  double sum = std::pow(ea/a,2) +  std::pow(eb/b,2.) - (2*rho/(a*b));
  double xmin = -b/(2*a);
  double x2 = xmin *(1 + sum);
  double evval = f->Eval(x2) - f->Eval(xmin) ;

  double ef = (sval/vval)*sqrt( std::pow(evval/vval,2)+ std::pow(esval/sval,2)); 

  // drawing stuff
  f2->SetLineColor(2);
  f2->SetLineWidth(3);
  f2->Draw("SAME");
  fhisto->Draw("HISTO");
  f2->Draw("SAME");
  f->Draw("SAME");

  canvas->SaveAs("./Plots/PeakToValley.pdf");

  // fill result
  res->peak.value = f->GetParameter(1); 
  res->peak.error = f->GetParError(1); 
  res->peakToValley.value = sval/vval;
  res->peakToValley.error = ef;
  //res->peak = sval;

  
 return res;
}

void GetPeakToValley(string rootFileName = "./outputFile.root",
		     int histoVersion = 1){
  
  
  // 
  // histoVersion = 0;
  // rootFileName = "PMT_NB0066_HV1820_Analysis.root";
  
  Result * res = fitModel(rootFileName,
			  histoVersion);
  cout << endl;
  cout << "peak           = " << res->peak.value << " (" << res->peak.error << ") " << endl;
  cout << "peak to valley = " << res->peakToValley.value << " (" << res->peakToValley.error << ") " << endl;
  
  
}

