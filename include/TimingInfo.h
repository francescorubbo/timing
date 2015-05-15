#ifndef TIMINGINFO_H
#define TIMINGINFO_H

#include "fastjet/PseudoJet.hh"

using namespace fastjet;

class TimingInfo : public PseudoJet::UserInfoBase{
 public:
 
 TimingInfo(const int & pdg_id_in,
	    const int & charge_in,
	    const int & pythia_id_in,  
	    const double & pv_in, 
	    const bool & pileup_in,
	    const double & pt_in,
	    const double & time_in,
	    const double & abstime_in,
	    const double & pixel_id_in=0,
	    const int & pixel_num_in=-1) 
   :_pdg_id(pdg_id_in),
    _charge(charge_in),
    _pythia_id(pythia_id_in), 
    _pv(pv_in),
    _pileup(pileup_in),
    _pt(pt_in),
    _time(time_in),
    _abstime(abstime_in),
    _pixel_id(pixel_id_in),
    _pixel_num(pixel_num_in){}
  
  int pdg_id() const { return _pdg_id;}
  int charge() const { return _charge;}
  int pythia_id() const {return _pythia_id;}
  double pixel_id() const { return _pixel_id;}
  int pixel_num() const {return _pixel_num;}
  bool pileup() const { return _pileup;}
  double pv() const { return _pv;}  
  double pt() const { return _pt;}
  double time() const { return _time;}
  double abstime() const { return _abstime;}
  bool isGhost() const { return (_pixel_id != 0); }
 protected:
  int _pdg_id;         // the associated pdg id
  int _charge;         // the particle charge
  int _pythia_id;  // index in pythia.event
  double _pv;  // the particle pv
  bool _pileup; //true if pileup, false if truth
  double _pt;  //pt
  double _time;  //corrected time
  double _abstime; //time from vertex
  double _pixel_id; //Pixel particle found in
  int _pixel_num; //Number of particle in pixel
};

#endif
