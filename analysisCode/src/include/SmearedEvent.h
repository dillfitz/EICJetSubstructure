#ifndef SMEAREDEVENT_H
#define SMEAREDEVENT_H

#include <utility>
#include <tuple>

#include <eicsmear/smear/EventS.h>
#include <eicsmear/erhic/EventPythia.h>

#include <TLorentzVector.h>

#include "JetDef.h"
#include "SoftDropJetDef.h"
#include "BreitFrame.h"

using PseudoJetVec = std::vector<fastjet::PseudoJet>;
using TLorentzVectorVec = std::vector<TLorentzVector>;
using JetConstPair = std::pair<TLorentzVector, std::vector<TLorentzVector>>;
using JetConstVec = std::vector<JetConstPair>;
using TLorentzPair = std::pair<TLorentzVector, TLorentzVector>;
using TLorentzPairVec = std::vector<TLorentzPair>;

class SmearedEvent {
  
 public:
  SmearedEvent(){}
 SmearedEvent(erhic::EventPythia &truthEvent, Smear::Event &smearEvent, std::vector<int> chadChildIndices)
   : m_truthEvent(&truthEvent)
    , m_smearEvent(&smearEvent)
    , m_chadChildIndices(chadChildIndices)
    {}

  ~SmearedEvent(){}

  /// Main workhorse function, which is called from event loop
  void processEvent();
  void setVerbosity(int verb) { m_verbosity = verb; }
  TLorentzVector getExchangeBoson();

  void setMaxPartEta(double eta){m_maxPartEta = eta;}
  void setMinPartPt(double pt){m_minPartPt = pt;}
  void setScatteredLepton();
  void setSmearedParticles();
  TLorentzPairVec getMatchedParticles();

  PseudoJetVec getRecoJets(fastjet::ClusterSequence *cs,
			   JetDef jetDef);
  PseudoJetVec getRecoSoftDropJets(PseudoJetVec recoJets, 
				   SoftDropJetDef sdJetDef);

  std::vector<PseudoJetVec> matchTruthRecoJets(PseudoJetVec truthjets, 
					       PseudoJetVec recojets);

  void useBreitFrame(bool yesorno) { m_breitFrame = yesorno; }

  // bool D0kpiNoSmearFilter();

 private:
  /// Need truth event for identifying only final state particles
  erhic::EventPythia *m_truthEvent;
  Smear::Event *m_smearEvent;

  double m_maxPartEta;
  double m_minPartPt;

  const Smear::ParticleMCS *m_scatLepton;
  bool m_breitFrame;
  std::vector<PseudoJetVec> m_matchedJets;
  std::vector<int> m_chadChildIndices;
  
  PseudoJetVec m_particles;
  PseudoJetVec m_truthParticles;
  int m_verbosity = 0;

 

};


#endif
