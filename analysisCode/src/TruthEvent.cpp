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

  if (m_verbosity == -4)
    cout << "New D0 Event... " << endl;

  BreitFrame breit(*m_truthEvent);
  int chadIndex = -99;
  vector<int> chadChildIndices;
  for(int part = 0; part < m_truthEvent->GetNTracks(); ++part)
    {

      /// Skip the beam
      if( part < 3 )
	continue; 

      const Particle *truthParticle = m_truthEvent->GetTrack(part);

      // Identify a charm hadron from a hard scattered charm quark //
      if (truthParticle->GetParentIndex() == 10 && abs(truthParticle->GetPdgCode()) == 4 && truthParticle->GetChild1Index() != 0 ) 
	{
	  // Note the indexing for GetTrack is different than the Pythia listing (it indexes from 0), so we map accordingly. 
	  for (int child = truthParticle->GetChild1Index() - 1; child < truthParticle->GetChildNIndex(); ++child)
	    {
	 
	      const Particle *chad = m_truthEvent->GetTrack(child);
	      //cout << "PIDs from hard scattered c quark... " << chad->GetPdgCode() << endl;
	      int chadPid = abs(chad->GetPdgCode());
	      if (chadPid == 421)
		{
		  chadIndex = chad->GetIndex();
		  CharmDecayTagger( chad, chadChildIndices );	  
		}
	    }
	}

      // Aside from the charm hadron, we only want final state particles
      if(truthParticle->GetStatus() != 1 && truthParticle->GetIndex() != chadIndex)
	continue;
    
      /// Skip the scattered electron, since it is special
      if(truthParticle->GetE() == m_scatLepton->GetE())
	continue;

      /// Check that eta is within nominal detector acceptance
      if(fabs(truthParticle->GetEta()) > 3.5)
	continue;
      if(truthParticle->GetPt() < 0.25)
      	continue;

      bool chadChildren = false;
      for (int i = 0; i<chadChildIndices.size(); ++i)
      {
        if( truthParticle->GetIndex() == chadChildIndices.at(i) )
	  {
	    chadChildren = true;      
	    m_chadChildIndices.push_back( truthParticle->GetIndex() );
	  }
      }
      if ( chadChildren )
	continue;

      if(m_verbosity == -4)
	{
	  std::cout << "Truth (lab) : " <<truthParticle->Id() 
		    << " " <<truthParticle->GetPx() << " " 
		    << truthParticle->GetPy() << " " << truthParticle->GetPz()
		    << " " << truthParticle->GetE() << " " <<truthParticle->GetIndex() <<std::endl;	  
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
		    << " " << partFourVec->E() << " mass : " << partFourVec->M() << std::endl;	  
	}

      m_partIndices.push_back( truthParticle->GetIndex() );
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

  PseudoJetVec charmJets = CharmJetTagging(selectJets);

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

bool TruthEvent::disD0Event()
{
  if (m_truthEvent->GetProcess() != 99 )
    return false;

  for(int part = 0; part < m_truthEvent->GetNTracks(); ++part)
    {    
      const Particle *truthParticle = m_truthEvent->GetTrack(part);

      // Identify a charm hadron from a hard scattered charm quark //
      if (truthParticle->GetParentIndex() == 10 && abs(truthParticle->GetPdgCode()) == 4 && truthParticle->GetChild1Index() != 0 ) 
	{
	  // Note the indexing for GetTrack is different than the Pythia listing (it indexes from 0), so we map accordingly. 
	  for (int child = truthParticle->GetChild1Index() - 1; child < truthParticle->GetChildNIndex(); ++child)
	    {
	      const Particle *chad = m_truthEvent->GetTrack(child);
 
	      //cout << "PID TEST... " << chad->GetPdgCode() << endl;

	      int chadPid = abs(chad->GetPdgCode());
	      if (chadPid == 421)
		{
		  bool passedCharmChildCuts;
		  passedCharmChildCuts = CharmDecayFilter( chad );

		  if (!passedCharmChildCuts )
		    return false;
		  else
		    return true;
			
		}
	    }
	}
    }

    return false;  
}

bool TruthEvent::disD0kpiEvent()
{
  if (m_truthEvent->GetProcess() != 99 )
    return false;

  for(int part = 0; part < m_truthEvent->GetNTracks(); ++part)
    {    
      const Particle *truthParticle = m_truthEvent->GetTrack(part);

      // Identify a charm hadron from a hard scattered charm quark //
      if (truthParticle->GetParentIndex() == 10 && abs(truthParticle->GetPdgCode()) == 4 && truthParticle->GetChild1Index() != 0 ) 
	{
	  // Note the indexing for GetTrack is different than the Pythia listing (it indexes from 0), so we map accordingly. 
	  for (int child = truthParticle->GetChild1Index() - 1; child < truthParticle->GetChildNIndex(); ++child)
	    {
	      const Particle *chad = m_truthEvent->GetTrack(child);
 
	      //cout << "PID TEST... " << chad->GetPdgCode() << endl;

	      int chadPid = abs(chad->GetPdgCode());
	      if (chadPid == 421)
		{
		  bool passedCharmChildCuts;
		  passedCharmChildCuts = D0kpiDecayFilter( chad );

		  if (!passedCharmChildCuts )
		    return false;
		  else
		    return true;
			
		}
	    }
	}
    }

    return false;  
}

bool TruthEvent::CharmDecayFilter( const Particle *part )
{

  for (int childIndex = part->GetChild1Index() - 1; childIndex < part->GetChildNIndex(); ++childIndex)
    {
      const Particle *child = m_truthEvent->GetTrack(childIndex);		      
      if (child->GetStatus() != 1)
	{
	  
	  CharmDecayFilter( child ); 
	}
      
      else
	{
	  if(fabs(child->GetEta()) > 3.5 && child->GetPt() < 0.25)
	    return false;	
	}
    }

  return true;
} 

bool TruthEvent::D0kpiDecayFilter( const Particle *part )
{

  if (part->GetNChildren() != 2) 
    { return false; }
  else
    for (int childIndex = part->GetChild1Index() - 1; childIndex < part->GetChildNIndex(); ++childIndex)
      {
	const Particle *child = m_truthEvent->GetTrack(childIndex);		      
	if ( abs(child->GetPdgCode())!= 211 &&  abs(child->GetPdgCode()) != 321 ) 
	  {
	    return false;
	  }
	
	else
	  {
	    if(fabs(child->GetEta()) > 3.5 || child->GetPt() < 0.25)
	      return false;	
	  }
      }
  
  return true;
} 

void TruthEvent::CharmDecayTagger( const Particle *part, vector<int> &childIndices )
{
  for (int childIndex = part->GetChild1Index() - 1; childIndex < part->GetChildNIndex(); ++childIndex)
    {
      const Particle *child = m_truthEvent->GetTrack(childIndex);		      
      if (child->GetStatus() != 1)
	{
	  CharmDecayTagger( child, childIndices ); 
	}
      else
	{
	  childIndices.push_back( child->GetIndex() );

	  if (m_verbosity == -4 )
	    std::cout << "Charm Child: " << " PID : " << child->Id()
		      << " " <<child->GetPt() << " " << child->GetEta() 
		      << " " <<" Index : " << child->GetIndex() << std::endl;
	}
    }

  return;
} 

PseudoJetVec TruthEvent::CharmJetTagging(PseudoJetVec truthJets)
{
  PseudoJetVec charmJets;
  std::vector<PseudoJet> cons;
  std::vector<int> consPids;

  for(int jet = 0; jet < truthJets.size(); ++jet)
    {
      cons.clear();
      consPids.clear();

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
	      //if(truthParticle->GetStatus() != 1)
	      //continue;

	      if ( truthParticle->GetPx() == constituent.px() &&
		  truthParticle->GetPy() == constituent.py() &&
		  truthParticle->GetPz() == constituent.pz() &&
		  truthParticle->GetE() == constituent.E() )
		{
		  //cout << "Cons Pid : " << truthParticle->GetPdgCode() << endl;
		   consPids.push_back(abs(truthParticle->GetPdgCode()));		   
	
		}	      
	    }      
	}
            
      // Lets just tag D0s since we can readily reconstruct these
      if  (std::find(consPids.begin(), consPids.end(), 421) != consPids.end())
	{
	  charmJets.push_back(truthJets.at(jet));
	}
      
    }

  return charmJets;
}
