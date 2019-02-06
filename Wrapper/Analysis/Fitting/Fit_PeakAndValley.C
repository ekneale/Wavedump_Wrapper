#include<iostream>
#include<fstream>
#include<vector>
#include"TMath.h"

#include"TFile.h"
#include"TTree.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

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

// Result

typedef struct {

double value;
double error;

  

} ValueWithError;

typedef struct {

  ValueWithError peak;
  ValueWithError peakToValley;
  ValueWithError noise;
  ValueWithError ped;
  ValueWithError mu;
  ValueWithError valley;
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



TH1D* findHisto(TFile *input){
  TH1D* histo = (TH1D*)input->Get("SPE0");
  if (histo) return histo;
  histo = (TH1D*)input->Get("SPE1");
  if (histo) return histo;
  histo = (TH1D*)input->Get("SPE2");
  if (histo) return histo;
  histo = (TH1D*)input->Get("SPE3");
  if (histo) return histo;
  std::cout << "no histo " << std::endl;
  return 0;
}

//Result* fitModel(TH1F* fhisto, double minval = -500 , double maxval = 2000 ){
Result* Fit_PeakAndValley(TH1F*  fhisto, 
			  double minval = -500 , 
			  double maxval = 2000 ){


  TCanvas* can = new TCanvas("can","can",800,600);
  Result* res = new Result();
  //gROOT -> ProcessLine( ".x ./mattStyle.C" );
  
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //
  //TFile *input=new TFile(file.c_str());
  //TH1D* histo = findHisto(input);
  //if (histo ==0) return 0;
  
  //TH1F* fhisto = h2h(histo);
  fhisto->GetXaxis()->SetTitle("charge [mV ns]");
  fhisto->GetYaxis()->SetTitle("Counts");
  fhisto->GetYaxis()->SetTitleOffset(1.3);
  fhisto->GetYaxis()->SetTitleFont(132);
  
  // intial guesses with TSpectrum
  InitParams params = initializeFit(fhisto);
  double valleyPos = std::get<5>(params);
  double peakPos = std::get<0>(params)+ std::get<1>(params);

  // gaussian fit to max
  TH1F* thisto = (TH1F*)fhisto->Clone("thisto");
  thisto->Fit("gaus", "","", peakPos - 100, peakPos +120 );
  TF1* gf1 = thisto->GetFunction("gaus");
  TF1* f1 = (TF1*)gf1->Clone();
  double sval1 = f1->GetParameter(1);
  std::cout << " s val " << peakPos <<  " " << sval1 << std::endl; 
  
  fhisto->Fit("gaus", "","", sval1 - 100, sval1 +100 );

  TF1* gf = fhisto->GetFunction("gaus");
  TF1* f = (TF1*)gf->Clone("gaus");
  double sval = f->Eval(f->GetParameter(1));
  double esval = TMath::Abs(f->Eval(gf->GetParError(1) + f->GetParameter(1) ) - sval) ;
  f->SetLineColor(2);
  f->SetLineWidth(3);
  f->Draw("SAME");

  
  fhisto->Fit("pol2", "S", "",valleyPos - 40 , valleyPos +50);
  TF1* f2 = fhisto->GetFunction("pol2");
  double a = f2->GetParameter(2); double b = f2->GetParameter(1);
  double xmin = -b/(2*a);
  
  TFitResultPtr pres= fhisto->Fit("pol2", "S", "",xmin - 30 , xmin +50);
  TFitResult* polres = pres.Get();
  double ea = polres->Error(2); double eb = polres->Error(1);

  // get correlation
  TMatrixDSym cm = polres->GetCorrelationMatrix();
  TMatrixDSym cov = polres->GetCovarianceMatrix();
  double rho = cov[1][2];

 
  // find min of parabola
  f2 = fhisto->GetFunction("pol2");
  a = f2->GetParameter(2); b = f2->GetParameter(1);
  double vval = f2->GetMinimum( valleyPos -50, valleyPos +50);
  double sum = sqrt(std::pow(ea/a,2) +  std::pow(eb/b,2.) - (2*rho/(a*b)));
  xmin = -b/(2*a);
  double x2 = xmin *(1 + sum);
  double evval = f->Eval(x2) - f->Eval(xmin) ;

  double ef = (sval/vval)*sqrt( std::pow(evval/vval,2)+ std::pow(esval/sval,2)); 

  double pzero = fhisto->Integral(1, fhisto->FindBin(f2->GetMinimumX( valleyPos -50, valleyPos +50)))/fhisto->GetEntries();
  
  
  // drawing stuff
  f2->SetLineColor(2);
  f2->SetLineWidth(3);
  // f2->Draw("SAME");
  fhisto->SetMaximum(100e3);
  fhisto->Draw("HISTO");
  f2->Draw("SAME");
  f->Draw("SAME");
  

  // fill result
  res->peak.value = f->GetParameter(1); 
  res->peak.error = f->GetParError(1); 
  res->peakToValley.value = sval/vval;
  res->peakToValley.error = ef;
  res->mu.value = -log(pzero);
  res->valley.value = f2->GetMinimumX( valleyPos -50, valleyPos +50);
  
  // noise gauss
  std::cout << "fit the noise " << std::endl;
  TF1* noise = new TF1("ngaus","gaus", -40, 40);
  noise->SetParameter(1,0);
  noise->SetParameter(2,10);
  double maxN = TMath::Max( std::get<0>(params) + 5, f2->GetMinimumX( valleyPos -50, valleyPos +50) - 20 );
  fhisto->Fit(noise, "","", -40,  f2->GetMinimumX( valleyPos -50, valleyPos +50) - 20 );
  res->noise.value = noise->GetParameter(2);
  res->noise.error = noise->GetParError(2);

  res->ped.value = noise->GetParameter(1);
  res->ped.error = noise->GetParError(1);
  
  std::string name = std::string(fhisto->GetName()) + ".png" ;
  
  name = "./Plots/" + name;
  
  can->Print(name.c_str());

  std::string epsname = std::string(fhisto->GetName()) + ".eps" ;
  
  epsname = "./Plots/" + epsname;
  
  can->Print(epsname.c_str());
    
 return res;
}


