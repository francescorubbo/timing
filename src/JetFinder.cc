#include "JetFinder.h"

JetFinder::JetFinder(double etaMin, double etaMax, double pixelSize=1e-4, double R=0.4){
  
  _etaMin= (etaMin > 0) ? etaMin : 0;
  if(etaMax <= etaMin){
    cerr << "Invalid Eta Limits " << etaMin << " -> " << etaMax << "Passed to JetFinder" << endl;
    exit(20);
  }
  _etaMax= etaMax;
  PixelSize(pixelSize);
  _R=R;

  const double grid_spacing(0.6);

  //suppress fastjet banner
  fastjet::ClusterSequence::set_fastjet_banner_stream(NULL);

  jetDef.reset(new JetDefinitio(fastjet::antikt_algorithm, _R, fastjet::E_scheme, fastjet::Best));
  active_area.reset(new AreaDefinition(fastjet::active_area));
  bge.reset(new GridMedianBackgroundEstimator(_etaMax, grid_spacing));
  select_fwd.reset(new Selector(SelectorAbsRapRange(_etaMin,_etaMax)));
}

void JetFinder::PixelSize(double pixelSize){
  _pixelSize = (pixelSize > 1e-6) ? pixelSize : 1e-6;
}

void JetFinder::selectJets(JetVector &particlesForJets, JetVector &selectedJets){
  
  fastjet::ClusterSequenceArea clustSeq(particlesForJets, (*jetDef), (*active_area));
  bge->set_particles(particlesForJets);
  fastjet::Subtractor subtractor(bge.get());
  
  JetVector inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(10.));
  JetVector subtractedJets = subtractor(inclusiveJets);

  selectedJets.clear();
  selectedJets = (*select_fwd)(subtractedJets);
}

void JetFinder::selectSegmentedJets(JetVector particlesForJets, JetVector selectedJets){
  selectJets(particlesForJets,selectedJets);
}
