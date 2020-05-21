#include "include/TruthEvent.h"

#include <fastjet/Selector.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include <iostream>
#include <algorithm>
#include <vector>


void TruthEvent::processEvent()
{

  setScatteredLepton();
  setTruthParticles();

}

TLorentzVector TruthEvent::getExchangeBoson()
{
  TLorentzVector *vector = new TLorentzVector( m_truthEvent->ExchangeBoson()->Get4Vector());
  
  if(m_breitFrame){
    
    BreitFrame breit(*m_truthEvent);
    breit.labToBreitTruth( vector );
    
  }
    
  return *vector;

}
bool TruthEvent::passCuts()
{
  double y = m_truthEvent->GetTrueY();
  double x = m_truthEvent->GetTrueX();
  double q2 = m_truthEvent->GetTrueQ2();
  
  if(m_verbosity > 3)
    {
      std::cout << "Truth x, q2, and y: " << x << "," << q2
		<< "," << y << std::endl;
    }
  
  return x > m_minX && y > m_minY && y < m_maxY && q2 > m_minq2;
}
void TruthEvent::setScatteredLepton()
{
  m_scatLepton = m_truthEvent->ScatteredLepton();
  if(m_verbosity > 1)
    {
      std::cout<< "Scattered lepton : "<< m_scatLepton->GetPx() << " "
	       << m_scatLepton->GetPy() << " " << m_scatLepton->GetPz()
	       << " " << m_scatLepton->GetE() << std::endl;
    }
}

void TruthEvent::setTruthParticles()
{

  BreitFrame breit(*m_truthEvent);
  for(int part = 0; part < m_truthEvent->GetNTracks(); ++part)
    {
      /// Skip the beam
      if( part < 3 )
	continue;
      
      const Particle *truthParticle = m_truthEvent->GetTrack(part);
      /// only want final state particles
      if(truthParticle->GetStatus() != 1)
	continue;
      
      /// Skip the scattered electron, since it is special
      if(truthParticle->GetE() == m_scatLepton->GetE())
	continue;

      /// Check that eta is within nominal detector acceptance
      if(fabs(truthParticle->GetEta()) > 3.5)
	continue;
      if(truthParticle->GetPt() < 0.25)
	continue;

      if(m_verbosity > 2)
	{
	  std::cout << "Truth (lab) : " <<truthParticle->Id() 
		    << " " <<truthParticle->GetPx() << " " 
		    << truthParticle->GetPy() << " " << truthParticle->GetPz()
		    << " " << truthParticle->GetE() << std::endl;	  
	}

      // Transform Particle 4 Vectors to the Breit Frame 
      TLorentzVector *partFourVec = new TLorentzVector( truthParticle->PxPyPzE() );
      if(m_breitFrame)
	breit.labToBreitTruth( partFourVec );
      
      if(m_verbosity > 0)
	{
	  std::cout << "Truth  : " <<truthParticle->Id() 
		    << " " <<partFourVec->Px() << " " 
		    << partFourVec->Py() << " " << partFourVec->Pz()
		    << " " << partFourVec->E() << std::endl;	  
	}

      
      m_particles.push_back(fastjet::PseudoJet(partFourVec->Px(),
					       partFourVec->Py(),
					       partFourVec->Pz(),
					       partFourVec->E()));
      
    }

  return;
}

PseudoJetVec TruthEvent::getTruthJets(fastjet::ClusterSequence *cs, 
				       JetDef jetDef)
{
  /// Create the cluster sequence
  cs = new fastjet::ClusterSequence(m_particles, jetDef.getJetDef());
  
  PseudoJetVec allTruthJets = fastjet::sorted_by_pt(cs->inclusive_jets());

  if(m_verbosity > 1)
    {
      std::cout << "Finding jets for jet pT>" << jetDef.getMinJetPt() 
		<< " and |eta|<" << jetDef.getMaxJetRapidity() << std::endl;
    }
  
  /// Make eta/pt selections
  fastjet::Selector selectPt = fastjet::SelectorPtMin(jetDef.getMinJetPt());
  fastjet::Selector selectEta = fastjet::SelectorAbsRapMax(jetDef.getMaxJetRapidity());
  fastjet::Selector select = selectPt and selectEta;
  
  PseudoJetVec selectJets = select(allTruthJets);

  PseudoJetVec charmJets = CharmTagging(selectJets);

  // return selectJets;
  return charmJets;

}

PseudoJetVec TruthEvent::getTruthSoftDropJets(PseudoJetVec truthJets, SoftDropJetDef sdJetDef)
{

  fastjet::contrib::SoftDrop sd(sdJetDef.getSoftDrop());
  
  PseudoJetVec softDropJets;
 
  for(int jet = 0; jet < truthJets.size(); ++jet)
    {
      fastjet::PseudoJet softDropJet = sd(truthJets[jet]);
      softDropJets.push_back(softDropJet);
    }

  return softDropJets;

}



bool TruthEvent::CharmEvent()
{
  std::vector<int> pids; 
  for(int part = 0; part < 11; ++part)
    {    
      const Particle *truthParticle = m_truthEvent->GetTrack(part);
      pids.push_back(abs(truthParticle->GetPdgCode()));
      
    }

  if ( std::find(pids.begin(), pids.end(), 4) != pids.end() )
    return true;
  else 
    return false;
  
}

bool TruthEvent::disCharmEvent()
{
  std::vector<int> pids; 
  if (m_truthEvent->GetProcess() != 99 )
    return false;

  for(int part = 0; part < 11; ++part)
    {    
      const Particle *truthParticle = m_truthEvent->GetTrack(part);
      pids.push_back(abs(truthParticle->GetPdgCode()));
      
    }

  if ( std::find(pids.begin(), pids.end(), 4) != pids.end() )
    return true;
  else 
    return false;
  
}


bool TruthEvent::pgfCharmEvent()
{
  std::vector<int> pids; 
  if ( (m_truthEvent->GetProcess() != 135 ) &&  (m_truthEvent->GetProcess() != 136 ) )
    return false;

  for(int part = 0; part < 11; ++part)
    {    
      const Particle *truthParticle = m_truthEvent->GetTrack(part);
      pids.push_back(abs(truthParticle->GetPdgCode()));
      
    }

  if ( std::find(pids.begin(), pids.end(), 4) != pids.end() )
    return true;
  else 
    return false;
  
}

void TruthEvent::PrintCharmEvent()
{
  cout << "Process ID : " << m_truthEvent->GetProcess() << endl;
  for(int part = 0; part < 11; ++part)
    {
      const Particle *truthParticle = m_truthEvent->GetTrack(part);
      cout << "PID : " << truthParticle->GetPdgCode() << endl;
    }

}

PseudoJetVec TruthEvent::CharmTagging(PseudoJetVec truthJets)
{
  PseudoJetVec charmJets;
  std::vector<PseudoJet> cons;
  //std::vector<std::pair<int, int>> parent_ids;
  std::vector<int> parent_ids;

  for(int jet = 0; jet < truthJets.size(); ++jet)
    {
      cons.clear();
      parent_ids.clear();

      cons = truthJets.at(jet).constituents();
      if ( m_verbosity == -4 ) 
	{
	  cout << "new jet : ";
	  cout << "num constituents : " << cons.size() << endl;
	}
      for (int con = 0; con < cons.size(); ++con)
	{
	  PseudoJet constituent;
	  constituent = cons.at(con);
	  for (int part = 0; part < m_truthEvent->GetNTracks(); ++part)
	    {
	      //cout << cons.size() << endl;

	      /// Skip the beam
	      if( part < 3 )
		continue;
      
	      const Particle *truthParticle = m_truthEvent->GetTrack(part);
	      /// only want final state particles
	      if(truthParticle->GetStatus() != 1)
		continue;

	      if ( truthParticle->GetPx() == constituent.px() &&
		  truthParticle->GetPy() == constituent.py() &&
		  truthParticle->GetPz() == constituent.pz() &&
		  truthParticle->GetE() == constituent.E() )
		{
		  // cout << "We found a match!" << endl;
		  const Particle *parentPart = m_truthEvent->GetTrack(truthParticle->GetParentIndex());
		  /*
		  if ( parentPart->GetNChildren() == 0 )
		    {
		      const Particle *grandparentPart= m_truthEvent->GetTrack(parentPart->GetParentIndex());
		      std::pair<int,int> parent_data(parentPart->GetPdgCode(),
						     parentPart->GetNChildren());

		      
		    }
		  else
		  */

		  //std::pair<int,int> parent_data(parentPart->GetPdgCode(), parentPart->GetNChildren());

	      parent_ids.push_back(abs(parentPart->GetPdgCode()));

		}
	      
	    }
      
	}

      if (m_verbosity == -4 )
	{
	  
	  for (int matchedpart=0; matchedpart < parent_ids.size(); ++matchedpart)
	    cout << parent_ids.at(matchedpart) << endl;
	  bool check = ((std::find(parent_ids.begin(), parent_ids.end(), 411) != parent_ids.end()) ||
			(std::find(parent_ids.begin(), parent_ids.end(), 413) != parent_ids.end()) ||
			(std::find(parent_ids.begin(), parent_ids.end(), 421) != parent_ids.end()) ||
			(std::find(parent_ids.begin(), parent_ids.end(), 423) != parent_ids.end()) ||
			(std::find(parent_ids.begin(), parent_ids.end(), 431) != parent_ids.end()) || 
			(std::find(parent_ids.begin(), parent_ids.end(), 433) != parent_ids.end()) ||
			(std::find(parent_ids.begin(), parent_ids.end(), 4122)!= parent_ids.end()) ||
			(std::find(parent_ids.begin(), parent_ids.end(), 4114)!= parent_ids.end()) ||
			(std::find(parent_ids.begin(), parent_ids.end(), 4214)!= parent_ids.end()) ||
			(std::find(parent_ids.begin(), parent_ids.end(), 4222)!= parent_ids.end()) || 
			(std::find(parent_ids.begin(), parent_ids.end(), 4224)!= parent_ids.end()) || 
			(std::find(parent_ids.begin(), parent_ids.end(), 4132)!= parent_ids.end()) );
	  cout << "Is it in le list?  "   << check << endl;
	  	  
	}

      if  ((std::find(parent_ids.begin(), parent_ids.end(), 411) != parent_ids.end()) ||
	   (std::find(parent_ids.begin(), parent_ids.end(), 413) != parent_ids.end()) ||
	   (std::find(parent_ids.begin(), parent_ids.end(), 421) != parent_ids.end()) ||
	   (std::find(parent_ids.begin(), parent_ids.end(), 423) != parent_ids.end()) ||
	   (std::find(parent_ids.begin(), parent_ids.end(), 431) != parent_ids.end()) || 
	   (std::find(parent_ids.begin(), parent_ids.end(), 433) != parent_ids.end()) ||
	   (std::find(parent_ids.begin(), parent_ids.end(), 4122)!= parent_ids.end()) ||
	   (std::find(parent_ids.begin(), parent_ids.end(), 4114)!= parent_ids.end()) ||
	   (std::find(parent_ids.begin(), parent_ids.end(), 4214)!= parent_ids.end()) ||
	   (std::find(parent_ids.begin(), parent_ids.end(), 4222)!= parent_ids.end()) || 
	   (std::find(parent_ids.begin(), parent_ids.end(), 4224)!= parent_ids.end()) || 
	   (std::find(parent_ids.begin(), parent_ids.end(), 4132)!= parent_ids.end()) ) 
	{
	  charmJets.push_back(truthJets.at(jet));
	}
	
    }

  return charmJets;
}
