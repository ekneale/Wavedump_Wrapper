#define PMTAnalyser_cxx
#include "PMTAnalyser.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TBranch.h"

// !!! IN PROGRESS
Int_t PMTAnalyser::DarkRate(Float_t threshold = 10)
{
  if (rawRootTree == 0) return -1;
  
  int verbosity = 0;
  
  Long64_t ientry;
  
  Long64_t nentries = rawRootTree->GetEntriesFast();
  Long64_t nDark    = 0 ;
  
  cout << endl;
  cout << " Calculating Dark Rate " << endl;
  cout << " " << nentries << " entries " << endl;

  for (Long64_t jentry=0; jentry < nentries; jentry++) {
    
    ientry = LoadTree(jentry);
    
    if (ientry < 0) break;
    rawRootTree->GetEntry(jentry);   
    
//     for( int iSample = 0 ; iSample < NSamples; iSample++){
      
//     }

    if( peakV_mV > threshold )
      nDark++;

    if(verbosity > 0){
      cout << endl;
      cout << " jentry        = " << jentry        << endl;
    }
    
  }  

  //!!!! Temporary
  Float_t darkRate = (Float_t)nDark/nentries/220.*1.0e9;
  
  if(verbosity > 0){  
    cout << endl;
    cout << " nentries = " << nentries << endl;
    cout << " nDark    = " << nDark    << endl;
    cout << " rate     = " << darkRate << endl;
  }
  
  return darkRate;
}


// !!! IN PROGRESS
Int_t PMTAnalyser::Make_hFixed_Filtered(){

  if (rawRootTree == 0) return -1;

  TCanvas * canvas = new TCanvas();
  
  int verbosity = 1;

  double waveformDuration = NSamples * nsPerSample;

  Int_t nBinsX  = 512;
  Float_t rangeQ[2] = {-500.,1500.};
  
  TH1D * hWave = new TH1D("hWave","hWaveform;Time /ns;ADC counts",
			NSamples, 0., waveformDuration);
  
  Long64_t ientry;
  Long64_t nentries = rawRootTree->GetEntriesFast();

  //!!! Test
  nentries = 300000;
  
  double aoff = 8700;
  
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    
    ientry = LoadTree(jentry);
    if (ientry < 0) break;
    rawRootTree->GetEntry(jentry);   
    
    if( jentry%100000==0 )
      cout << " entry = " << jentry << endl; 
    
    for( int iSample = 0 ; iSample < NSamples; iSample++){
      
      hWave->SetBinContent(iSample+1,
			   (double)(aoff - waveform[iSample]));
      
    } // end: for( int iSample = ...


    // Fill unfiltered histogram

    // Check this entry's waveform
    if( IsCleanFFTWaveform(hWave) ){
      
      // Fill filtered histograms
      
    }
  } //end of: for (Long64_t jentry = 0; jentry

  // Draw plots and append to root file
  
  return 1;

}


Bool_t PMTAnalyser::IsCleanFFTWaveform(TH1D * hWave){
  
  Bool_t waveformIsClean = kFALSE;

  TH1F* hWaveFFT  = new TH1F("hWaveFFT","hWaveFFT",
			    NSamples, 0, NSamples);
  
  hWave->FFT(hWaveFFT ,"MAG");
  
  // delete zero frequency data
  hWaveFFT->SetBinContent(1,0.) ;
  
  if(hWaveFFT->GetMaximumBin() == 2)
    waveformIsClean = kTRUE;
  
  hWaveFFT->Delete();
    
  return waveformIsClean;
}


 // 
 // Member Function to test filtering noisy waveforms 
 //
 // Author 
 //  Steve Quillin
 //
 // Translated  22 02 19 
 //  gary.smith@ed.ac.uk
//
//
// Description 
//
//  Smallest non-zero frequency component in bin #2 
//  at a frequency of  1/110 samples * 2ns/sample ~ 4.55 MHz
//  bins 2 - 54 hold magnitudes at multiples of this frequency
//  bin 55  Nyquist frequency component at 250 MHz = 1/2*2ns  
//
//  Note that the magnitude of the zero frequency FFT bin is 
//  a measure of the trace baseline value (multiplied by the 
//  number of samples)
//  Propose using only data traces where the peak FFT component 
//  falls in the lowest non-zero bin (bin #2)

Int_t PMTAnalyser::Make_FFT_Histos()
{
  if (rawRootTree == 0) return -1;

  TCanvas * canvas = new TCanvas();
  
  int verbosity = 1;

  double waveformDuration = NSamples * nsPerSample;
  
  TH1D* hWave = new TH1D("hWave","hWaveform;Time /ns;ADC counts",
			NSamples, 0., waveformDuration);
  
  TH1F* hWaveFFT_Accumulate  = new TH1F("hWaveFFT_Accumulate","hWaveFFT_Accumulate",
					NSamples, 0, NSamples);
  
  TH1D* hMaxADC = new TH1D("hMaxADC","Waveform ADC Maximum;maxADC;counts",
			   1000, 200.,1200.); 
  
  TH1D* hMaxADC_Filtered = new TH1D("hMaxADC_Filtered",
				    "hMaxADC_Filtered;maxADC;counts",
				    1000, 200.,1200.);
  
  Long64_t ientry;
  Long64_t nentries = rawRootTree->GetEntriesFast();
  
  double aoff   = 8700.; // arbitrary offset / baseline
  double fftoff = 0;     // fft offset / baseline

  cout << endl;
  cout << " Performing FFT " << endl;
  
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    
    ientry = LoadTree(jentry);
    if (ientry < 0) break;
    rawRootTree->GetEntry(jentry);   
    
    if( jentry%100000==0 )
      cout << " entry = " << jentry << endl; 
    
    for( int iSample = 0 ; iSample < NSamples; iSample++){
      
      hWave->SetBinContent(iSample+1,
			  (double)(aoff - waveform[iSample]));
      
      if( verbosity > 1){
	cout << endl;
	cout << " waveform[" << iSample << "] = " 
	     <<   waveform[iSample] << endl;
      }
	
    } // end: for( int iSample = ...
    
    hMaxADC->Fill( hWave->GetMaximum() );

    // Accumulate FFT of waveform in hWaveFFT_Accumulate
    hWave->FFT(hWaveFFT_Accumulate ,"MAG");
    
    // Check this entry's waveform
    if( IsCleanFFTWaveform(hWave) ){
      hMaxADC_Filtered->Fill( hWave->GetMaximum());
    }
    
    fftoff = hWaveFFT_Accumulate->GetBinContent(1)/(double)NSamples; 

    // delete zero frequency data
    hWaveFFT_Accumulate->SetBinContent(1,0.) ;
    
    if(verbosity > 1){
      cout << endl;
      cout << " jentry = " << jentry << endl;
      cout << " fftoff = " << fftoff << endl;
    }
    
  } //end of: for (Long64_t jentry = 0; jentry
  
  if(verbosity > 0){  
    cout << endl;
    cout << " nentries = " << nentries << endl;
  }

  canvas->Divide(2,1);

  canvas->cd(1);
  hMaxADC->Draw();
  
  hMaxADC_Filtered->SetLineColor(kRed);
  hMaxADC_Filtered->Draw("same");

  canvas->cd(2);  
  hWaveFFT_Accumulate->Draw();
    
  canvas->SaveAs("canvas.pdf");
  
  return 1;
}
