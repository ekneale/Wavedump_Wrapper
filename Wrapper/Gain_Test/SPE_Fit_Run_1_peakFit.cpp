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


#include "TObject.h"
#include "TLegend.h"


using namespace std;
using namespace TMath;
using namespace RooFit;

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


typedef std::tuple<double,double,double,double,double,double> InitParams;

//========== Find approximate peak, pedestal and valley values for fit paramaters ====================

InitParams initializeFit(TH1F* h){

  TSpectrum *spec = new TSpectrum(3,3);
  Int_t nfound = spec->Search(h,3,"goff",0.0002);
  std::cout << "found peaks " << nfound << std::endl;

  /*** Find the approximate charge at the maximum of the SPE peak ***/
  Float_t *peaks;
  peaks = spec->GetPositionX();
  float sigPeak; float pedPeak;
  std::cout << peaks[0] << " " << peaks[1]  << " " << peaks[2] << std::endl;
  
  /*** Apply workarounds if TSpectrum fails for the varying HV steps ***/

  if (nfound == 1){
      double XMin = 28.;
      double XMax = 1500.;
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

  /****Finding the Valley*****/
  for(int i = PedPed; i< SigSig ;i++){
    Compare[0] = h->GetBinContent(i);

    if(Compare[0] < Compare[1]){
      Compare[1] = Compare[0]; //finds minimum number of events
      ValleyBin = i; //sets bin number for minimum
    }//end if
  }//end for
  std::cout << "Ped pos " << PedPed << "Sig pos " << SigSig << "Valley bin " << ValleyBin << std::endl;

  double ssignal = (sigPeak -  h->GetBinCenter(ValleyBin))/3.;
  double sped = (h->GetBinCenter(ValleyBin) - peaks[0])/10.;

   /*****Finding Pedestal events****/
  int Noise = 0; //to hold pedestal events
  for(int i = 0;i< ValleyBin;i++){
    Noise += h->GetBinContent(i); // total noise events
  }//end for

  /** Approximate fraction of photons per event**/
  int Events = h->GetEntries();  // # of events
  double Ratio = (double) Noise / (double) Events; // Ratio of Random events to photon events

  std::cout << "Peaks " << peaks[0] << " " << sigPeak << " " << sped <<   " " << ssignal << " "   <<  Ratio << std::endl;

  double valleyPos = h->GetBinCenter(ValleyBin) ;
  std::cout << "Valley = " << valleyPos << std::endl;
  
  return InitParams(pedPeak,sigPeak-pedPeak,sped,ssignal,Ratio,valleyPos);
} // initializeFit


//========== Create output spectra ====================

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

//========== Carry out fit ====================

Result* fitModel(TH1F* fhisto,
            double minval = -500,
            double maxval = 2000){

  Result* res = new Result();

  TCanvas * canvas = new TCanvas("Canvas","Canvas");
  fhisto->GetXaxis()->SetTitle("charge [mV ns]");
  fhisto->GetYaxis()->SetTitle("Counts");
  fhisto->GetYaxis()->SetTitleOffset(1.5);
  fhisto->GetYaxis()->SetTitleFont(132);

  // intial guesses with TSpectrum
  InitParams params = initializeFit(fhisto);
  double valleyPos = std::get<5>(params);
  double peakPos = std::get<0>(params)+ std::get<1>(params);
  // gaussian fit to max
  std::cout << "Fitting gaussian to max..." << std::endl;
  double gaussMin = valleyPos ;
  double gaussMax = peakPos + 120;
  fhisto->Fit("gaus", "","", gaussMin, gaussMax );
  TF1* gf = fhisto->GetFunction("gaus");
  TF1* f = (TF1*)gf->Clone();
  double sval = f->Eval(f->GetParameter(1));
  double esval = TMath::Abs(f->Eval(gf->GetParError(1) + f->GetParameter(1) ) - sval) ;
  f->SetLineColor(2);
  f->SetLineWidth(3);
  f->Draw("SAME");

  //polynomial fit to valley
  std::cout << "Fitting polynomial to valley..." << std::endl;
  TFitResultPtr pres= fhisto->Fit("pol2", "S", "",valleyPos - 50 , valleyPos +50);
  TFitResult* polres = pres.Get();
  double ea = polres->Error(2); double eb = polres->Error(1);

  // get correlation
  TMatrixDSym cm = polres->GetCorrelationMatrix();
  TMatrixDSym cov = polres->GetCovarianceMatrix();
  double rho = cov[1][2];


  // find min of parabola
  std::cout << "Getting value for valley" << std::endl;
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

  /*** Apply workaround if fit fails ***/
  if (res->peak.value < 0){
    double XMin = 28.;
    double XMax = 1500.;
    fhisto->GetXaxis()->SetRange(XMin,XMax);
    int binMax = fhisto->GetMaximumBin();
    res->peak.value = fhisto->GetBinCenter(binMax);
    res->peak.error = 0;
  }

 return res;

} // fitModel


//========== Gain curve fit ====================

inline double PowerFunc(double x, double k, double n)
{
  return pow((k*x),n)*10.00;
} // PowerFunc

double fitPow(double *x, double *k)
{
  return PowerFunc(x[0], k[0], k[1]);
} // fitPow



//============================================================

int main(int argc,char **argv){	
  
	
  /*** Read in the HV data ***/

  int nominalHV;
  int pmt;
  int loc;
  int run =1;
  int pmtHV[5];
  
  string hvfile = "../HVScan.txt";
  ifstream file(hvfile.c_str());
  string hvdat;
 
  vector<int> PMT_number(125,0), HV(125,0);
  vector< vector <int> > HVstep;
  vector<int> step(5,0);
  for (int i=0; i<125; i++)
    HVstep.push_back(step);
  
  for (int i=0; i<125; i++){
    for (int j=0; j<7; j++){
      file >> hvdat;
      int pmt_info =atof(hvdat.c_str());
      if (j==0){
	PMT_number[i]=pmt_info;
      }
      if (j!=0 && j!=6){
	HVstep[i][j-1]=pmt_info;
      }
      if (j==6){
	HV[i]=pmt_info;
      }
    }

  }


  /*** Assign PMT numbers and locations of PMTs in the rig ***/
  static const int nPMTsA = 80;
  static const int nPMTsB = 20;
  static const int nPMTs  = nPMTsA + nPMTsB;

  int  pmtAList[nPMTsA] = {83 , 88,108,107,
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
               130,131,132,133};


  int  pmtBList[nPMTsB] = {102,149,150,152,
               9  , 10, 12, 14,
               43 , 47, 48, 49,
               30 , 31, 32, 33,
               134,135,136,138};

  int  locAList[4] = {0,1,2,3};
  int  locBList[4] = {4,5,6,7};

  
  char histname[200]= "";
 

  for (int iPMT = 0 ; iPMT < nPMTs ; iPMT++){

    if( iPMT < nPMTsA ){
        pmt = pmtAList[iPMT];
        loc = locAList[iPMT%4];
    }
    else{
        pmt = pmtBList[iPMT-nPMTsA];
        loc = locBList[(iPMT-nPMTsA)%4];
    }


    /*** Assign HV values for each PMT ***/ 
    for (int i=0;i<125; i++){
     
      if (pmt == PMT_number[i]){
        nominalHV = HV[i];
        printf("nominal HV is %d \n",nominalHV);
      }
    
      for(int h=0; h<5; h++){
        if (pmt==PMT_number[i]){
          pmtHV[h] =HVstep[i][h];
  	      printf("HV %d location  %d Test %d \n",pmtHV[h],loc,h);
	    }
	  }
    }



    //========== Fit the SPE Spectrum ====================
    
    double hvVals[5]; double gainVals[5]; double gainValsError[5];
  
    for (int r=0;r<5;r++){ 
    
      int hv = r+1; // gain test number
    
      sprintf(histname, "/data/kneale/Wavedump_Wrapper/RawRootData/Run_%d_PMT_%d_Loc_%d_HV_%d.root",run,pmt,loc,hv); 
      TFile s(histname);
      s.ls();

      char root_name[30];
      sprintf(root_name, "hQ_Peak_Run_%d_PMT_%d_Loc_%d_HV_%d",run,pmt,loc,hv);
      TH1D *speData = (TH1D*)s.Get(root_name);

      TH1F* fhisto = h2h(speData);

      printf("Getting data from SPE spectrum...\n");

	  
  	  /***Find the value at the maximum***/
      Result * res = fitModel(fhisto);
      float signal = res->peak.value;
      float signalError = res->peak.error;
        
      cout << endl;
      cout << "peak           = " << signal << " (" << signalError << ") " << endl;
        
      printf(" voltage is  %d , charge is %f \n\n\n\n",pmtHV[r],signal); 


      /*** Calculate the gain ***/
  	  float amplification = 10.0; //amplification via amplifier
      float splitter = 2.0; //correction for use of splitter box
	  float impedence = 50.0;
	  gainVals[r] = signal*10e-12*splitter/amplification/impedence/(1.602*10e-19)/1e7; //conversion from charge in mVns to gain
	  gainValsError[r] = signalError*10e-12*splitter/amplification/impedence/(1.602*10e-19)/1e7; //conversion from charge in mVns to gain
      hvVals[r] = pmtHV[r];

    } 
  
    cout << "out" << endl;
    
    // ========== Make the HV fit ====================

    TString canvasNameTemp    = "", canvasName    = "";

    canvasNameTemp = "Canvas_Run_%d_PMT_%d_Loc_%d_HV_G";
    canvasName.Form(canvasNameTemp,run,pmt,loc);

    TCanvas *canvas=new TCanvas(canvasName,canvasName,1);

    TGraph *Gain;

    /*** Plot gain vs voltage for each PMT ***/
    Gain = new TGraph(5,hvVals,gainVals);

    Gain->SetMarkerStyle(2);
    Gain->SetMarkerSize(1);
		
    /*** Fit the gain curve ***/
    double fitMin = hvVals[0];
    double fitMax = hvVals[4];
    TF1 *f14 = new TF1("f14",fitPow,fitMin,fitMax,2);
    f14->SetParameter(0,10);
    f14->SetParameter(1,10);
    TFitResultPtr tfrp14=Gain->Fit("f14","RSE");

    Gain->Draw("AP");
    Gain->GetYaxis()->SetTitle("Gain (10^7)");
    Gain->GetXaxis()->SetTitle("Applied Voltage (V)");
	
	
    /*** Calculate the operating voltage for 10^7 gain ***/
    double operatingHV = pow(1.0/10.,(1./f14->GetParameter(1)))/(f14->GetParameter(0)); //inverse of fit function with y=1 (1e7 gain)
    double operatingHVError = abs(pow(1.0/10.,(1./f14->GetParameter(1)))/(f14->GetParameter(0)) - pow(1.0/10.,(1./(f14->GetParameter(1)+f14->GetParError(1))))/(f14->GetParameter(0)+f14->GetParError(0)));
    printf("\n\n\n\n\n Operating voltage for 10^7 Gain for PMT%d: %f  +/- %f \n\n\n\n\n",pmt, operatingHV, operatingHVError );

    //check if voltages.root exists - if not, create it and set up tree and branches
    //if it exists, add to the branches

    Double_t chi2 = f14->GetChisquare();
    Double_t NDf = f14->GetNDF();
    Double_t prob = f14->GetProb();


    /*** Write ntuples to file ***/
    ifstream fileStream("voltages_peakFit.root");
    if(!fileStream.good()){
      TFile *outfile = new TFile("voltages_peakFit.root","RECREATE");

      TNtuple *voltages = new TNtuple("voltages","voltages","pmt:operatingHV:operatingHVError:nominalHV:chi2:NDf:prob");
      voltages->Fill(pmt,operatingHV,operatingHVError,nominalHV,chi2,NDf,prob);
      voltages->Write();
      outfile->Close();
      delete outfile;
    }  

    else{
      TFile *outfile = new TFile("voltages_peakFit.root","UPDATE");
      TNtuple *voltages = (TNtuple*)outfile->Get("voltages");
      voltages->Fill(pmt,operatingHV,operatingHVError,nominalHV,chi2,NDf,prob);
      voltages->Write("",TObject::kOverwrite);
      outfile->Close();
      delete outfile;
    }
  

    canvasName = "./Plots/";
    canvasName += canvas->GetName();
    canvasName += ".png";
    canvas->SaveAs(canvasName);
  }
  return 0;
}





