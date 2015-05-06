#include "TimingTracker.h"

int main(int argc, char* argv[]){

  float eta,phi;
  if(argc >= 2){
    eta=atof(argv[1]);
    phi=atof(argv[2]);
  }
  else{
    cerr << "Enter Eta/Phi" << endl;
    return 1;
  }
  double z=7.26;

  cout << eta << "\t" << phi << "\t-\t";
  double etaSign = (eta < 0) ? -1 : 1;
  auto converted=EtaPhi_to_xy(etaSign*eta,phi,z);
  cout << converted.first << "\t" << converted.second << "\t-\t";
  auto back=xy_to_EtaPhi(converted.first,converted.second,z);
  cout << etaSign*back.first << "\t" << back.second << endl;

  return 0;
}
