#ifndef _fileNameParse_h
#define _fileNameParse_h 1
char parse(std::string f, std::string s1, std::string s2){
 
  int pos1 = f.find(s1,0);
  int pos2 = f.find(s2,pos1);

  // "Test" not found for case "HV"
  if     (pos1 == -1){
    s1 = "HV";
    pos1 = f.find(s1,0);
  }
  else if(pos2 == -1){
    s2 = "HV";
    pos2 = f.find(s2,pos1);
  }

  std::string str = f.substr(pos1,pos2 - pos1);
  
  int first = str.find("_");
  str = str.substr(first+1);
  str = str.substr(0, str.size()-1);
  
  return str[0];
}

string GetFileName(string filePath){

  int pos1 = filePath.find("Run_",0);
  string fileName = filePath.substr(pos1,filePath.size() - pos1);  
  return fileName;
}

string GetFileID(string fileName){
  return fileName.substr(0,fileName.size() - 5);
}

string GetTreeName(string filePath){
  
  string treeName = "Events_" + GetFileID(GetFileName(filePath));
  
  return treeName; 
}

char pmtID(string filename){
  return parse(filename,"PMT_","Loc");
}

char run(string filename){
  return parse(filename,"Run_","PMT");
}

char location(string filename){
  return parse(filename,"Loc_","Test");
}

char test(string filename){ 
  
  char ctr = parse(filename,"Test_","root");

  if( ctr > '0' && ctr < '6' ) 
    return 'G';
  else
    return ctr;
} 

char HVStep(string filename){ 
  
  char step = parse(filename,"Test_","root");

  if( step > '0' && step < '6')
    return step;
  else 
    return '0';
} 

#endif

/* void stringmess(){

  string runname = "Run_4_PMT_90_Loc_0_Test_N.root";
  cout <<pmtID(runname)  <<  " " << run(runname) <<  " " << location(runname)<< endl; 
}
*/
