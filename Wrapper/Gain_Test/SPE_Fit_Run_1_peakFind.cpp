#include <cmath>
#include <cstdlib>
#include <fstream>

#include "TApplication.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TObject.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TString.h"
#include "TSpectrum.h"
#include "TNtuple.h"

#include "TH1D.h"
#include "TFile.h"

using namespace std;
using namespace TMath;






//======= Gain Curve Fitters =============
inline double PowerFunc(double x, double k, double n)
{
  return pow((k*x),n)*10.00;
}

double fitPow(double *x, double *k)
{
  return PowerFunc(x[0], k[0], k[1]);
}

//=========================================================================================================================================

int main(int argc,char **argv){	
  
  // ******************
  // * Initialization *
  // ******************
	
  //Read in the HV data ====================================================================================
  string hvfile = "../HVScan.txt";
  ifstream file(hvfile.c_str());
  string hvdat;
 
  vector<int> PMT_number(125,0), HV(125,0);
  vector< vector <int> > HVstep;
  vector<int> step(5,0);
  for (int i=0; i<125; i++)
    HVstep.push_back(step);
  

  int pmtHV[5];
  
 
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
  
  //========================================================================================================
  int nominalHV;
  int pmt;
  int loc;
  int run = 1;

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
    //Determing the PMT number, location and the applied Voltage=====================================================

  for (int iPMT = 0 ; iPMT < nPMTs ; iPMT++){

    if( iPMT < nPMTsA ){
        pmt = pmtAList[iPMT];
        loc = locAList[iPMT%4];
    }
    else{
        pmt = pmtBList[iPMT-nPMTsA];
        loc = locBList[(iPMT-nPMTsA)%4];
    }
  
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


    
  //======================================================================================================
  
  
  // ***************
  // * Set up ROOT *
  // ***************

    double hvVals[5]; double gainVals[5];
    double gainValues[5]; //creates an empty array for the gain value for each PMT and each step
    // Fitting the SPE Spectrum =======================================================================
  
    for (int r=0;r<5;r++){ 
    
      TCanvas *tc=new TCanvas("SPE","SPE",1);
      //tc->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
      tc->SetGrid();
    
      // Create canvas, allowing for window close.
    
    
      // *************************
      // * Create output spectra *
      // *************************
    
      int hv = r+1; // gain test number
    
      sprintf(histname, "/data/kneale/Wavedump_Wrapper/RawRootData/Run_%d_PMT_%d_Loc_%d_HV_%d.root",run,pmt,loc,hv); 
    
      //printf("Voltage: 
      TFile s(histname);
      s.ls();
    
      char root_name[30];
      sprintf(root_name, "hQ_Fixed_Run_%d_PMT_%d_Loc_%d_HV_%d",run,pmt,loc,hv);
      TH1D *speData = (TH1D*)s.Get(root_name);

    printf("Getting data from SPE spectrum...\n");
    speData->GetYaxis()->SetTitle("Counts ");
    speData->GetYaxis()->SetTitleOffset(1.5);
    speData->GetXaxis()->SetTitle("Charge (mv*ns)");	
//	speData->GetXaxis()->SetRange(-150,1500);

	// Draw spectrum.
	speData->Draw("same");

	// Update canvas.
	tc->Update();
	tc->Paint();
	tc->Draw();
	tc->Modified();
		
		
		
	/***Find the value at the maximum***/
    TSpectrum *spec = new TSpectrum(3,3);
    Int_t nfound = spec->Search(speData,2,"goff",0.0002);
    std::cout << "found peaks " << nfound << std::endl;

    // returns positions of Pedestal and Signal approx
    Float_t *peaks;
    peaks = spec->GetPositionX();
    float signalMax;


    std::cout << peaks[0] << " " << peaks[1]  << " " << peaks[2] << std::endl;

    // returns the charge at the maximum of the SPE peak
    if (nfound == 1){
      double XMin = 28.;
      double XMax = 1500.;
      speData->GetXaxis()->SetRange(XMin,XMax);
      int binMax = speData->GetMaximumBin();
      signalMax = speData->GetBinCenter(binMax);
      std::cout << "charge is " << signalMax << " for 1-peak spectrum \n" << std::endl;
    }
    

    if (nfound == 3 && peaks[0] < -5){
      signalMax = peaks[2];
      std::cout << "charge is " << signalMax << " for 3-peak spectrum \n" << std::endl;
      
    }

    if (nfound == 3 && peaks[0] > -5){
        signalMax = peaks[1];
        std::cout << "charge is " << signalMax << " for 3-peak spectrum \n" << std::endl;
    }

    if (nfound ==2) {
      signalMax = peaks[1];
      std::cout << "charge is " << signalMax << " for 2-peak spectrum \n" << std::endl;
      
    } 

    printf(" voltage is  %d , charge is %f \n\n\n\n",pmtHV[r],signalMax); 
	float amplification = 10.0; //amplification via amplifier
	float splitter = 2.0; //correction for use of splitter box
	float impedence = 50.0;
	double pmtGain = signalMax*10e-12*splitter/amplification/impedence/(1.602*10e-19)/1e7; //conversion from charge in mVns to gain
	gainValues[r] = pmtGain; //create a list of the gain values

    hvVals[r] = pmtHV[r];
    gainVals[r] = gainValues[r];


    /*
    // save the SPE plots - removed this step as the histograms are generated by BinToRoot	
	TString plotName  = "./Plots/PMT_%d_%d.png";
	plotName.Form("./Plots/PMT_%d_%d.png",pmt,pmtHV[r]);
		
	cout << " pmt     = " << pmt << endl;
	cout << " pmtHV[" << r << "] = " << pmtHV[r] << endl;
	cout << " plotName = \n \n \n \n " << plotName << endl;

	tc->SaveAs(plotName);
    */
  }
  cout << "out" << endl;
    
  TString canvasNameTemp    = "", canvasName    = "";

  canvasNameTemp = "Canvas_Run_%d_PMT_%d_Loc_%d_HV_G";
  canvasName.Form(canvasNameTemp,run,pmt,loc);

 

  TCanvas *canvas=new TCanvas(canvasName,canvasName,1);
  
  // Making the HV fit ========================================================================

  TGraph *Gain;

  //Plot gain vs voltage for each PMT
  Gain = new TGraph(5,hvVals,gainVals);

  Gain->SetMarkerStyle(2);
  Gain->SetMarkerSize(1);
		
  //fit a curve of form y = 10*(kx)^n to the data
  double fitMin = hvVals[0];
  double fitMax = hvVals[4];
  TF1 *f14 = new TF1("f14",fitPow,fitMin,fitMax,2);
  f14->SetParameter(0,10);
  f14->SetParameter(1,10);
  TFitResultPtr tfrp14=Gain->Fit("f14","RSE");

  //Draw the plot plus fit
  Gain->Draw("AP");
  Gain->GetYaxis()->SetTitle("Gain (10^7)");
  Gain->GetXaxis()->SetTitle("Applied Voltage (V)");
  //  Gain->SetTitle(Gain);
		
  //Bias for 10^7 GAIN ==========================================================================================================
  //not based on on a fit to the peak

  //calculate the operating voltage for 10^7 gain
  double operatingHV = pow(1.0/10.,(1./f14->GetParameter(1)))/(f14->GetParameter(0)); //inverse of fit function with y=1 (ie y= 1e7 gain)
  double operatingHVError = abs(pow(1.0/10.,(1./f14->GetParameter(1)))/(f14->GetParameter(0)) - pow(1.0/10.,(1./(f14->GetParameter(1)+f14->GetParError(1))))/(f14->GetParameter(0)+f14->GetParError(0)));
  printf("\n\n\n\n\n Operating voltage for 10^7 Gain for PMT%d: %f  +/- %f \n\n\n\n\n",pmt, operatingHV, operatingHVError );

  //check if voltages.root exists - if not, create it 
  //and set up tree and branches
  //if it exists, add to the branches
  Double_t chi2 = f14->GetChisquare();

//======= Write ntuples to file ===========

  ifstream fileStream("voltages.root");
  if(!fileStream.good()){
    TFile *outfile = new TFile("voltages.root","RECREATE");

    TNtuple *voltages = new TNtuple("voltages","voltages","pmt:operatingHV:operatingHVError:nominalHV:chi2");
    voltages->Fill(pmt,operatingHV,operatingHVError,nominalHV,chi2);
    voltages->Write();
    outfile->Close();
    delete outfile;
  }

  else{
    TFile *outfile = new TFile("voltages.root","UPDATE");
    TNtuple *voltages = (TNtuple*)outfile->Get("voltages");
    voltages->Fill(pmt,operatingHV,operatingHVError,nominalHV,chi2);
    voltages->Write("",TObject::kOverwrite);
    outfile->Close();
    delete outfile;
  }
	//=============================================================================================================================
  

  canvasName = "./Plots/";
  canvasName += canvas->GetName();
  canvasName += ".png";
  canvas->SaveAs(canvasName);
 }
  return 0;
}





