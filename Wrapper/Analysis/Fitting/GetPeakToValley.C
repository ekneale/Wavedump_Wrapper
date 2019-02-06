#include "Fit_PeakAndValley.C"
#include "TFile.h"
#include "TString.h"

using namespace std;

Result* Fit_PeakAndValley(TH1F*  fhisto, 
			  double minval, 
			  double maxval);

void GetPeakToValley(TString rootFileName = "Run_200_PMT_15_Loc_0_Test_S"){

  TString hName = "hQ_Fixed_" + rootFileName;
      
  rootFileName = rootFileName + ".root";

  //!!!!!
  // set to relevant local path
  TString RawRootDataDIR = "/Users/gsmith23/Desktop/Watchman/Testing/Wavedump_Wrapper/RawRootData/";
  
  rootFileName = RawRootDataDIR + rootFileName;
  
  TFile * rootFile = new TFile(rootFileName);
  
  Result * results = Fit_PeakAndValley((TH1F*)rootFile->Get(hName),
				       -500,
				       2000);
  
  cout << endl;
  cout << "peak           = " << results->peak.value         << " (" << results->peak.error         << ") " << endl;
  cout << "peak to valley = " << results->peakToValley.value << " (" << results->peakToValley.error << ") " << endl;
  
  
}

