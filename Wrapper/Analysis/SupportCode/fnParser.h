#ifndef _fnparse_h
#define _fnparse_h 1
int parse(std::string f, std::string s1, std::string s2){
 
  int pos1 = f.find(s1,0);
  int pos2 = f.find(s2,pos1);
  std::string str = f.substr(pos1,pos2 - pos1);
  int first = str.find("_");
  str = str.substr(first+1);
  str = str.substr(0, str.size()-1);
  return std::stoi(str);
  
}

int pmtID(std::string filename){
  return parse(filename,"PMT_","Loc");
}

int run(std::string filename){
   return parse(filename,"Run_","PMT");
}

int location(std::string filename){
  return parse(filename,"Loc_","Test");
}
#endif

/* void stringmess(){

  std::string runname = "Run_4_PMT_90_Loc_0_Test_N.root";
 
  
  std::cout <<pmtID(runname)  <<  " " << run(runname) <<  " " << location(runname)<< std::endl; 
}
*/
