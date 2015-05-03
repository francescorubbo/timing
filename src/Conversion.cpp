#include "TimingTracker.h"

int main(int argc, char* argv[]){

  float eta,phi;
  if(argc >= 2){
    eta=atof(argv[1]);
    phi=atof[argv[2]);
  }
  else{
    cerr << "Enter Eta/Phi" << endl;
  }

  cout << eta << "\t" << phi << "\t-\t";
  auto converted=EtaPhi_to_xy(eta,phi,1.2);
  cout << converted.first << "\t" << converted.second << "\t-\t";
  auto back=xy_to_EtaPhi(converted.first,converted.second,1.2);
  cout << back.first << "\t" << back.second << endl;

  return 0;
}
