/*
 * A program to calculate the gain and peak-to-valley ratio
 * at applied voltage
 *
 * Author
 * Liz Kneale
 * eskneale1@sheffield.ac.uk
 * 
 * Adapted from SPE_Fit.cpp by Tomi Akindele
 * and fit code by Matt Needham
 *
 * This program reads in the charge spectra output
 * from BintoRoot.cpp in Wavedump_Wraper/RawRootData, 
 * and outputs SPE charge, gain and peak-to-valley.
 *
 * Outputs data to screen and plots of SPE spectrum 
 * fit. Requires ./Plots/ in SPE_Test directory.
 *
 *
 * To build:
 *
 * First build the RooExpWindow roofit class:
 * make -f Makefile-1
 *
 * then:
 * make -f Makefile SPE_Fit_singlePMT
 *
 * To run:
 *
 * ./SPE_Fit_singlePMT
 *
 * then input the data as prompted.
 *
 * Data required: pmt, rig location, run, applied voltage
 *
 * Dependencies:
  * ROOT
 */





#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <tuple>
#include <vector>
#include <string>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TNtuple.h"
#include "TH1D.h"

#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDerivative.h"
#include "RooPolynomial.h"
#include "RooBifurGauss.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooExponential.h"
#include "RooExpWindow.h"
#include "RooMsgService.h"

#include "TObject.h"
#include "TLegend.h"


using namespace std;
using namespace TMath;
using namespace RooFit;

/****************** Set up for charge spectrum fit ********************/
typedef struct {

double value;
double error;

} ValueWithError;

typedef struct {

  ValueWithError ped;
  ValueWithError fped;
  ValueWithError pemean;
  ValueWithError pewidth;
  ValueWithError npe;
  ValueWithError peak;
  ValueWithError valley;
  ValueWithError peakToValley;
  ValueWithError fvalley;
  double rate;
  int id;

} Result;

void fillValueWithError(ValueWithError* val,TH1D* h){
  val->value = h->GetMean();
  val->error = h->GetRMS();
}

void fillValueWithError(ValueWithError* val,RooRealVar* var){
  val->value = var->getVal();
  val->error = var->getError();
}


typedef std::tuple<double,double,double,double,double,double> InitParams;

/******************* Find approximate peak, pedestal and valley values for fit **********************/

InitParams initializeFit(TH1F* h){

  TSpectrum *spec = new TSpectrum(3,3);
  Int_t nfound = spec->Search(h,3,"goff",0.0002);
  std::cout << "found peaks " << nfound << std::endl;

  /*** Find the approximate charge at the maximum of the SPE peak ***/
  Float_t *peaks;
  peaks = spec->GetPositionX();
  float sigPeak = 0; float pedPeak = 0;
  std::cout << peaks[0] << " " << peaks[1]  << " " << peaks[2] << std::endl;
  
  /*** Apply workarounds if TSpectrum fails for the varying HV steps ***/

  if (nfound == 1){
      double XMin = 28.;
      double XMax = 1200.;
      h->GetXaxis()->SetRange(XMin,XMax);
      int binMax = h->GetMaximumBin();
      sigPeak = h->GetBinCenter(binMax);
      pedPeak = peaks[0];
      std::cout << "charge is " << sigPeak << " for 1-peak spectrum \n" << std::endl;
  }

  if (nfound ==2) {
      sigPeak = peaks[1];
      pedPeak = peaks[0];
      std::cout << "charge is " << sigPeak << " for 2-peak spectrum \n" << std::endl;

  }

  if ( (nfound == 3) && (peaks[0] < -5)) {
      sigPeak = peaks[2];
      pedPeak = peaks[1];
      std::cout << "charge is " << sigPeak << " for 3-peak spectrum \n" << std::endl;

  }

  if (nfound == 3 && peaks[0] > -5 && peaks[1] > 100){
        sigPeak = peaks[1];
        pedPeak = peaks[0];
        std::cout << "charge is " << sigPeak << " for 3-peak spectrum \n" << std::endl;
  }

  if (nfound == 3 && peaks[0] > -5 && peaks[1] < 100){
        sigPeak = peaks[2];
        pedPeak = peaks[0];
        std::cout << "charge is " << sigPeak << " for 3-peak spectrum \n" << std::endl;
  }
  
  if (sigPeak < 20){
        sigPeak = 35;
        pedPeak = peaks[0];
        std::cout << "charge is " << sigPeak << " for very low-gain peak \n" << std::endl;
  }



  /*** Find the valley ***/
  Int_t SigSig = h->FindBin(sigPeak); //signal bin approx
  Int_t PedPed  = h->FindBin(pedPeak); //pedestal bin approx
  int PedMax = h->GetMaximum(); // pedestal events

  int Compare[2] = {PedMax,PedMax}; // array to store 2 numbers for comparison to determine which is larger
  Int_t ValleyBin = 0; //to hold valley bin number

  for(int i = PedPed; i< SigSig ;i++){
    Compare[0] = h->GetBinContent(i);

    if(Compare[0] < Compare[1]){
      Compare[1] = Compare[0]; //finds minimum number of events
      ValleyBin = i; //sets bin number for minimum
    }//end if
  }//end for

  std::cout << "Ped pos " << PedPed << "Sig pos " << SigSig << "Valley bin " << ValleyBin << std::endl;

  double ssignal = (sigPeak -  h->GetBinCenter(ValleyBin))/2.; // TODO 2 or 3?
  double sped = (h->GetBinCenter(ValleyBin) - pedPeak)/10.;

   /*****Finding Pedestal events****/
  int Noise = 0;
  for(int i = 0;i<ValleyBin;i++){
    Noise += h->GetBinContent(i);
  }//end for 
 
  /** Approximate fraction of photons per event**/
  int Events = h->GetEntries();  // # of events
  double Ratio = (double) Noise / (double) Events; // Ratio of Random events to photon events


  double valleyPos = h->GetBinCenter(ValleyBin) ;
  std::cout << "Peaks " << pedPeak << " " << sigPeak-pedPeak << " " << sped <<   " " << ssignal << " "   <<  Ratio << " " << valleyPos << std::endl;
  
  return InitParams(pedPeak,sigPeak-pedPeak,sped,ssignal,Ratio,valleyPos);
} // initializeFit


const char* varname(std::string header, std::string var){
  return (header + var).c_str();
}

/************************ Build the model ****************************/

// Fuction to model the PE peaks as sums of Gaussians - not currently used
RooAddPdf* makePGaussPDF(int n, RooRealVar* counts, RooFormulaVar* k,
                        RooRealVar* m2, RooRealVar* s2,  RooRealVar* pedm,
                        RooRealVar* peds,RooRealVar* m0, 
                        double mval, double sval){


  RooArgList pdfs; RooArgList coeffs ;
  int nlast = 10*n + 20;

  for (int m = 1; m < nlast; ++m ){


    std::string mform =  "pedmean + " + std::to_string(m) + "*kmean2";
    std::string mname = "kmean2_" + std::to_string(100*n)+ "_"+std::to_string(m);
    RooFormulaVar* mvar =  new RooFormulaVar(mname.c_str(),mname.c_str(),mform.c_str(),RooArgSet(*m2,*pedm));

    std::string sform =  "sqrt(pedsigma*pedsigma + " + std::to_string(m)+ "*ksigma2*ksigma2)";
    std::string sname = "ksigma2_" + std::to_string(100*n) +"_"+std::to_string(m);
    RooFormulaVar* svar =  new RooFormulaVar(sname.c_str(),sname.c_str(),sform.c_str(),RooArgSet(*s2,*peds));

    std::string gname = "gauss_k2_" +  std::to_string(100*n) + "_" +std::to_string(m);
    RooGaussian* gauss = new RooGaussian(gname.c_str(),gname.c_str(),*counts, *mvar, *svar);
    pdfs.add(*gauss);

    std::string fform = "TMath::Poisson(" + std::to_string(m) + "," +std::to_string(n) + "*k)"  ;
    std::string cname = "f_" + std::to_string(100*n)  + "_"+std::to_string(m);
    RooFormulaVar* cvar = new RooFormulaVar(cname.c_str(),cname.c_str(),fform.c_str(), *k);
    if (m != nlast -1) coeffs.add(*cvar);

  }
  std::string astring = "comp_" + std::to_string(100*n);
  RooAddPdf* apdf = new RooAddPdf(astring.c_str(),astring.c_str(),pdfs,coeffs);
  // apdf->Print();
  return apdf;
}
 
/*** construct the pdf ***/
RooAddPdf* makePMTPDF(RooRealVar* counts,double pmval, double psval, double psval2,  double mval, double sval, double f1peval,double vmval, int expvar, double vaval = 0.03, double fvval = 0.3){// counts,pedPeak,sigPeak-pedPeak,sped,ssignal,Ratio,valleyPos
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  std::cout << "pedestal mean " << pmval << "pedestal sigma " << psval << "SPE mean "  << mval << "SPE sigma "  << sval  << std::endl;
  std::cout << "SPE fraction " << f1peval << "valley min " << vmval << "valley alpha " << vaval << std::endl;
  
  /*** construct the pedestal pdf ***/ 
  RooRealVar* pedm = new RooRealVar("pedmean","pedmean",pmval, pmval - 5*psval, pmval + 5*psval );   // pedestal position
  RooRealVar* peds = new RooRealVar("pedsigma","pedsigma",psval,0,5*psval ); // pedestal sigma
  RooRealVar* peds2 = new RooRealVar("pedsigma2","pedsigma2",psval2,psval,10*psval ); // pedestal sigma
  
  RooRealVar* fped = new RooRealVar("fped","fped",0.9,0,1 );
  RooGaussian* pedgauss = new  RooGaussian("pedgauss","pedgauss", *counts, *pedm, *peds);
  RooGaussian* pedgauss2 = new  RooGaussian("pedgauss2","pedgauss2", *counts, *pedm, *peds2);

  RooAddPdf* pedpdf = new RooAddPdf("pedpdf", "pedpdf",RooArgList(*pedgauss,*pedgauss2), *fped);

  /*** construct the exponential fit to the thermionic emission ***/
  RooRealVar* vm = new RooRealVar("vmin","vmin",pmval+psval*expvar,0,14); 
  RooRealVar* va = new RooRealVar("valpha","valpha",vaval,0,0.06);
  RooExpWindow* vexp = new RooExpWindow("vexp","vexp",*counts,*va,*vm);

  /***   construct the single- and multiple-pe pdf ***/
  // first construct the Gaussians
  RooRealVar* m0 = new RooRealVar("mean0","mean0",mval,  mval - 2*sval,  mval  + 10*sval);   // peak position
  RooRealVar* s0 = new RooRealVar("sigma0","sigma0", sval,0, 10*sval);  // peak sigma 

  RooRealVar* npe = new RooRealVar ("npe","npe", f1peval,0.01,1); // number of photo electrons

  /*** create the sum-of-Gaussians fit - not currently used ***/
//  RooRealVar* km2 = new RooRealVar("kmean2","kmean2",mval/5,  0, mval  );
//  RooRealVar* ks2 = new RooRealVar("ksigma2","ksigma2",sval/sqrt(10.),1, sval);
//  RooFormulaVar* k = new RooFormulaVar("k","k", "mean0/kmean2",RooArgSet(*km2,*m0));

//  RooAddPdf* pG1 = makePGaussPDF(1,counts,k,km2,ks2,pedm,peds, m0 ,mval,sval);
//  RooAddPdf* pG2 = makePGaussPDF(2,counts,k,km2,ks2,pedm,peds,m0,mval,sval);

  std::string form1 = "mean0 + pedmean";
  RooFormulaVar* m =  new RooFormulaVar("mean","mean",form1.c_str(),RooArgSet(*m0,*pedm));
  std::string form2 =  "sqrt(sigma0*sigma0 + pedsigma*pedsigma)" ;
  RooFormulaVar* s =  new RooFormulaVar("sigma","sigma",form2.c_str(),RooArgSet(*s0,*peds));
  RooGaussian* gauss = new  RooGaussian("gauss1","gauss1", *counts, *m, *s);

  /*** option to fit 2pe and 3pe peaks ***/
  
  std::string form3 = "2*mean0 + pedmean";
  RooFormulaVar* m2 =  new RooFormulaVar("mean2","mean2",form3.c_str(),RooArgSet(*m0,*pedm));
  std::string form4 =  "sqrt(2*sigma0*sigma0 + pedsigma*pedsigma)" ;
  RooFormulaVar* s2 =  new RooFormulaVar("sigma2","sigma2",form4.c_str(),RooArgSet(*s0,*peds));
  RooGaussian* gauss2 = new  RooGaussian("gauss2","gauss2", *counts, *m2, *s2);
  /*
  std::string form5 = "3*mean0 +pedmean";
  RooFormulaVar* m3 =  new RooFormulaVar("mean3","mean3",form5.c_str(),RooArgSet(*m0,*pedm));
  std::string form6 =  "sqrt(3*sigma0*sigma0 + pedsigma*pedsigma)" ;
  RooFormulaVar* s3 =  new RooFormulaVar("sigma3","sigma3",form6.c_str(),RooArgSet(*s0,*peds));
  RooGaussian* gauss3 = new  RooGaussian("gauss3","gauss3", *counts, *m3, *s3);
  */

  /*** then get the coefficients ***/
  RooFormulaVar*  frac_ped1 = new RooFormulaVar("frac_ped1", "frac_ped1", "TMath::Exp(-npe)", *npe); //Poisson for 0 Photo electrons Poisson(0,npe)
  RooFormulaVar* frac_pe1 = new RooFormulaVar("frac_pe1", "frac_pe1",  "TMath::Poisson(1,npe) ", *npe); // ""   for 1   ""
  /*** option to add coefficients for the 2pe and 3pe peaks ***/
  
  RooFormulaVar* frac_pe2 = new RooFormulaVar("frac_pe2", "frac_pe2",  "TMath::Poisson(2,npe) ", *npe);
  /*
  RooFormulaVar* frac_pe3 = new RooFormulaVar("frac_pe3", "frac_pe3",  "TMath::Poisson(3,npe) ", *npe);
  */
  RooRealVar* frac_v = new RooRealVar("frac_v","frac_v",fvval, 0,1); //fraction of events in the valley region

  RooAddPdf* smodel = new RooAddPdf("smodel","smodel",RooArgList(*pedpdf,*gauss,*gauss2), RooArgList(*frac_ped1,*frac_pe1,*frac_pe2));
  RooAddPdf* model = new RooAddPdf("model","model",RooArgList(*smodel,*vexp),*frac_v);
  
  return model;
}

/*** construct the pdf with parameters extracted from the fit ***/
RooAddPdf* makePMTPDF(RooRealVar* counts,const RooArgList& fitpars){
 
 gErrorIgnoreLevel = kWarning;
 RooRealVar* fpedmean = (RooRealVar*)fitpars.find("pedmean");
 RooRealVar* pedsigma = (RooRealVar*)fitpars.find("pedsigma");
 RooRealVar* pedsigma2 = (RooRealVar*)fitpars.find("pedsigma2");
 RooRealVar* fmean = (RooRealVar*)fitpars.find("mean0");
 RooRealVar* fsigma = (RooRealVar*)fitpars.find("sigma0");
 RooRealVar* f1pe = (RooRealVar*)fitpars.find("npe");
 RooRealVar* vm = (RooRealVar*)fitpars.find("vmin");
 RooRealVar* va = (RooRealVar*)fitpars.find("valpha");
 RooRealVar* fv = (RooRealVar*)fitpars.find("frac_v");

 return makePMTPDF( counts,fpedmean->getVal(),pedsigma->getVal(), pedsigma2->getVal(),fmean->getVal(), fsigma->getVal(), f1pe->getVal(), vm->getVal(), va->getVal(),fv->getVal()); 
}

/*** construct the pdf with initparams ***/
RooAddPdf* makePMTPDF(RooRealVar* counts,InitParams& params,int expvar){
  return makePMTPDF(counts,std::get<0>(params),std::get<2>(params),2*std::get<2>(params), std::get<1>(params),std::get<3>(params),-log(std::get<4>(params)),std::get<5>(params),expvar);
}

/*** Get the results of the fit ***/
Result* propagateAndFill(RooRealVar* counts,RooAddPdf* model ,RooFitResult* fres){

 Result* res = new Result();
 
 /*** get the final fit parameters ***/
 
 const RooArgList& fitpars = fres->floatParsFinal();
 RooRealVar* fpedmean = (RooRealVar*)fitpars.find("pedmean");
 RooRealVar* fmean = (RooRealVar*)fitpars.find("mean0");
 RooRealVar* fsigma = (RooRealVar*)fitpars.find("sigma0");
 RooRealVar* fnpe = (RooRealVar*)fitpars.find("npe");

 if (fpedmean) fillValueWithError(&res->ped,fpedmean);
 if (fmean) fillValueWithError(&res->pemean,fmean);
 if (fsigma) fillValueWithError(&res->pewidth,fsigma);
 if (fnpe) fillValueWithError(&res->npe,fnpe);

 /*** now get the complicated ones that require sampling the fitted pdf/covariance ***/
 TH1D* histo = new TH1D("valley", "valley", 200, res->ped.value,  res->pemean.value ); histo->Sumw2();
 TH1D* histo2 = new TH1D("peak", "peak", 200,  res->pemean.value - res->pewidth.value ,  res->pemean.value +5* res->pewidth.value ); histo2->Sumw2();

 RooArgSet nset(*counts) ;


 for (int i = 0; i < 100; ++i){
   
   const RooArgList& sample = fres->randomizePars();
   RooAddPdf* pdf = makePMTPDF(counts,sample);
   RooMsgService::instance().setSilentMode(true);
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

   TF1* fmodel = pdf->asTF( *counts,fitpars,*counts);
   double vpos = fmodel->GetMinimumX(res->ped.value,res->pemean.value);
   double ppos = fmodel->GetMaximumX(res->pemean.value - res->pewidth.value, res->pemean.value + 5*res->pewidth.value);
   histo->Fill(vpos);
   histo2->Fill(ppos);

   counts->setRange("signal",vpos, 1000) ;

 }

 fillValueWithError(&res->valley,histo);
 fillValueWithError(&res->peak,histo2);

 return res;
}

/******************* Create the charge spectra **************************/

TH1F* h2h(TH1D* hold ){

  TH1F* h = new TH1F("histo", "histo", hold->GetNbinsX(),hold->GetXaxis()->GetXmin(),  hold->GetXaxis()->GetXmax());
  h->Sumw2();
  double sum =0;
  for (auto i = 1 ; i <= h->GetNbinsX(); ++i ){
    sum += hold->GetBinContent(i);
  }
  int entries = hold->GetEntries();
  for (auto i = 1 ; i <= h->GetNbinsX(); ++i ){
    int bincontent = TMath::Nint(entries*hold->GetBinContent(i)/sum);
    double valx = hold->GetBinCenter(i);
    for (int i = 0; i < bincontent; ++i){
      h->Fill(valx);
    }

  }
  std::cout << "histo " << h->GetEntries()<< " " <<  hold->GetEntries() <<  " " << sum <<  std::endl;
  return h;
} // h2h

/******************** Carry out the fit ******************************/

Result* fitModel(TH1F* fhisto, int pmt, int hv,
            double minval = -100,
            double maxval = 1800,
            double max = 10000){
//Result* res = new Result();
  TCanvas * canvas = new TCanvas("Canvas","Canvas");
  fhisto->GetXaxis()->SetTitle("charge [mV ns]");
  fhisto->GetYaxis()->SetTitle("Counts");
  fhisto->GetYaxis()->SetTitleOffset(1.5);
  fhisto->GetYaxis()->SetTitleFont(132);

  /*** make intial guess of parameters with TSpectrum ***/
  InitParams params = initializeFit(fhisto);

  RooRealVar* counts = new RooRealVar("charge", "charge", 0,minval, maxval);
  RooDataHist data("data", "dataset", *counts , fhisto);

  /*** fix the start of the exponential fit depending on HV step ***/
  int expvar;
  if (hv < 4){
    expvar = 4;
  }
  else {
    expvar = 2;
  }
  
  /*** build and fit the model starting from the initial parameters ***/
  RooAddPdf* model = makePMTPDF(counts,params,expvar);

  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  RooFitResult* fres = model->fitTo(data,Save()); //,RooFit::PrintLevel(-1),Warnings(kFALSE),PrintEvalErrors(-1));


  RooPlot* frame = counts->frame();
  data.plotOn(frame);

  model->plotOn(frame,Components("pedpdf"),LineColor(kOrange),LineStyle(2));
  model->plotOn(frame, Components("gauss1"),LineColor(kBlue), LineStyle(2));
  model->plotOn(frame, Components("gauss2"),LineColor(kGreen), LineStyle(2));
  model->plotOn(frame, Components("gauss3"),LineColor(kGreen), LineStyle(2));
  model->plotOn(frame, Components("vexp"),LineColor(kBlack),LineStyle(2));

  model->plotOn(frame, LineColor(kRed));

  std::cout << " plotting " << std::endl;
  TAxis* xachse = frame->GetXaxis();  TAxis* yachse = frame->GetYaxis();
  xachse->SetTitleFont (132);
  yachse->SetTitleFont (132);
  xachse->SetLabelFont (132);
  yachse->SetLabelFont (132);
  yachse->SetTitleOffset(1.1);
  xachse->SetTitleSize(0.06); 
  yachse->SetTitleSize(0.06); 
  xachse->SetLabelSize(0.06); 
  yachse->SetLabelSize(0.06); 
  yachse->SetTitle("Entries");
  xachse->SetTitle("Relative Charge");
  //frame->SetMaximum(max);
  frame->Draw();
//  canvas->SaveAs(Form("./Plots/FullFit/Fit_Run_1_PMT_%d_HV_%d.C",pmt,hv));
  gPad->SetLogy();  
  canvas->SaveAs(Form("./Plots/Fit_Run_1_PMT_%d_HV_%d_Test_S.png",pmt,hv));

  Result* res = propagateAndFill(counts,model,fres);
 
  return res;
} //fitModel


/**************** Gain calculation **********************************/

float GainCalc(double mVnsval){ 

  float amplification = 10.0; //amplification via amplifier
  float splitter = 2.0; //correction for use of splitter box
  float impedence = 50.0;
  return mVnsval*10e-12*splitter/amplification/impedence/(1.602*10e-19)/1e7; //conversion from charge in mVns to gain
}


/*********************************************************************/

int main(int argc,char **argv){	
 
	
  /*** Read in the HV data ***/

  int pmt;
  int loc;
  int run;
  int hv;
  float gain;
  float gainError;
  float peak2valley;
  float peak2valleyError;
  char histname[300]= "";



  /*** Determine the PMT number and applied voltage for each step ***/

  cout << "Input the PMT number \n" ;
  cout << "Note: please neglect the NB and the zeros before the number \n" <<endl;
  cin  >> pmt;
  cout <<endl;

  cout << "Input the location number \n" ;
  cin  >> loc;
  cout <<endl;

  cout << "Input the run number \n";
  cin  >> run;
  cout << endl;

  cout << "Input the applied voltage \n";
  cin  >> hv;
  cout << endl;


  /*** Read in and fit the charge Spectrum ***/
    
  
  sprintf(histname, "../../RawRootData/Run_%d_PMT_%d_Loc_%d_Test_S.root",run,pmt,loc); 
  TFile s(histname);
  s.ls();

  char root_name[50];
  sprintf(root_name, "hQ_Filter_Run_%d_PMT_%d_Loc_%d_Test_S",run,pmt,loc);
  TH1D *speData = (TH1D*)s.Get(root_name);

  TH1F* fhisto = h2h(speData);
  printf("Getting data from SPE spectrum...\n");

  /*** Find the SPE charge output ***/
  Result * res = fitModel(fhisto,pmt,hv);

  float signal = res->pemean.value - res->ped.value;
  float signalError = res->pemean.error;
        

  /*** Calculate the gain ***/
  gain = GainCalc(signal);
  gainError = GainCalc(signalError); 
  peak2valley = res->peak.value/res->valley.value;
  peak2valleyError = sqrt(pow(res->peak.error/res->peak.value,2) + pow(res->valley.error/res->valley.value,2));

  cout << endl;
  cout << "peak           = " << signal << " (" << signalError << ") " << endl;
        
  printf(" charge is %f, gain is %f x 10^7 +/- %f, peak to valley is  %f +/- %f for pmt %d at %dV \n\n\n\n",signal,gain,gainError,peak2valley,peak2valleyError,pmt,hv); 


  return 0;
}






