#include <string>
#include <vector>
#include<iostream>
#include<fstream>
#include<sstream>
#include<ostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"


bool dataLine(std::string& line){

  bool ok = true;
  std::string val = std::string(1,line.at(0));
  if (val == "(" || val == "," || val == "S" ||  val ==  "\""){
    ok = false;
  }
  return ok;
}

void getFiles(std::string input, std::vector<std::string>& vec){

  std::ifstream in;
  in.open(input.c_str());
  
  if(in.fail()){
    std::cout << "Failed to find file" <<std::endl;
    return;
  }
  
  // read it;
  std::string line;
  while(getline(in,line)){
    if (dataLine(line)){
      vec.push_back(line);
    }  //if
  } // while
}



void extractHistos(std::string filename, std::string dir , TFile* ofile){

  std::string fullpath =  dir + filename;
  std::cout << "Processing: " << fullpath << std::endl;
  TFile* file  = new TFile(fullpath.c_str());
  
  std::string histoname1 = "hQ_Fixed_" + filename.substr(0,filename.size() - 5); 
  std::string histoname2 = "hQ_Peak_" + filename.substr(0,filename.size() - 5);

  TH1F* h1 = (TH1F*) file->Get(histoname1.c_str()); 
  TH1F* h2 = (TH1F*) file->Get(histoname2.c_str());

  ofile->cd();
  h1->Write();
  h2->Write();
  file->Close();
}

void extract(std::string input = "list.txt",std::string dir = "/Disk/ds-sopa-group/PPE/Watchman/RawRootData/" , std::string oname = "SPEHistos.root"){

 TFile* ofile = new TFile(oname.c_str(), "RECREATE");
 std::vector<std::string> files; getFiles(input,files);
 for (auto ifile : files ){
   extractHistos(ifile.c_str(),dir , ofile);
 }

 ofile->Close();  
}


