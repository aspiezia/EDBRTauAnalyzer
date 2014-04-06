#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "TString.h"
typedef struct{
  TString name;
  int idx;
  int subidx;
  long long size;
} filedata;


bool decision(filedata f1, filedata f2){
  return (f1.size<f2.size); 
}

void Order(TString filelist, TString filename, bool print_uncut, TString FileName){

  ifstream input(filelist.Data());
  std::vector<filedata> file_vec; 
  filedata aux_data; 

  int name_sz = (int)filename.Sizeof();
  int n=0;
  while(!input.eof()){
  
    TString line;
    input >> line;

    if(line.IsAlnum()){
      aux_data.size = line.Atoll(); 
      file_vec.push_back(aux_data);  
      //if(aux_data.size!=0) std::cout << aux_data.size << " ";
      if(aux_data.size!=0) std::cout << aux_data.name << " ";
    }  
    if(!line.IsAlnum()){
      //aux_data.name = line.Data(); 
      //std::cout << aux_data.name << " ";
    } 
    if(n%2!=0){
      //cout<<n<<" "<<aux_data.size<<" "<<aux_data.name<<endl;
      //file_vec.push_back(aux_data);  
    }
    n=n+1;
    //cout<<endl; 

    if(line.Contains(filename.Data())){
      //file_vec.push_back(aux_data);  
      ////std::cout << line.Data() << std::endl;
      //TString str_aux = line;
      //int idx_chop = line.Index(" "); 
      //line.Remove(0,idx_chop);
      ////cout<<line<<endl;   
      //aux_data.name = line;
      //int idx_aux = str_aux.Index(filename.Data());
      //idx_aux+=name_sz;
      //str_aux.Remove(0,idx_aux-1);
      //str_aux.ReplaceAll(".root","");
      //TString str_aux2 = str_aux;
      //str_aux2.Remove(str_aux.Index("_"));
      //str_aux.Remove(0,str_aux.Index("_")+1);
      //
      //aux_data.idx = str_aux2.Atoi();
      //aux_data.subidx = str_aux.Atoi();
      
      //file_vec.push_back(aux_data);  
    }
  }

  std::sort(file_vec.begin(), file_vec.end(), decision);
  int vec_sz = (int)file_vec.size();
  
  ofstream myfile;
  myfile.open (FileName); 

  for(int i=0; i<vec_sz;i++){
    //std::cout << file_vec.at(i).size << " " << file_vec.at(i).name.Data() <<std::endl;
    if(file_vec.at(i).size!=0) myfile<<file_vec.at(i).size << /*" " << file_vec.at(i).name.Data()<<*/std::endl;
  }
  cout<<endl;
}
