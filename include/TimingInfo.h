#ifndef TIMINGINFO_H
#define TIMINGINFO_H

#include "fastjet/PseudoJet.hh"

using namespace fastjet;

class TimingInfo : public PseudoJet::UserInfoBase{
 public:
 TimingInfo(const int & pdg_id_in,const int & pythia_id_in,  const double & pv_in, const bool & pileup_in, const double & time_in) :
  _pdg_id(pdg_id_in),_pythia_id(pythia_id_in), _pv(pv_in),_pileup(pileup_in),_time(time_in){}
  int pdg_id() const { return _pdg_id;}
  int pythia_id() const {return _pythia_id;}
  double pv() const { return _pv;}
  bool pileup() const { return _pileup;}
  double time() const { return _time;}
 protected:
  int _pdg_id;         // the associated pdg id
  int _pythia_id;  // index in pythia.event
  double _pv;  // the particle pv
  bool _pileup;
  double _time;  //
};

#endif
