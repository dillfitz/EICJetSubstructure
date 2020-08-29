#include "include/SmearedEvent.h"

#include <fastjet/Selector.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include <iostream>
#include <random>

bool defaultPID = true;
bool requestedPID = false;

void SmearedEvent::processEvent()
{
  setScatteredLepton();
  setSmearedParticles();
  //setSmearedParticlesPid();
}

void SmearedEvent::setScatteredLepton()
{
  m_scatLepton = m_smearEvent->ScatteredLepton();
  if(m_verbosity > 1)
    {
      std::cout<< "Scattered lepton : "<< m_scatLepton->GetPx() << " "
	       << m_scatLepton->GetPy() << " " << m_scatLepton->GetPz()
	       << " " << m_scatLepton->GetE() << std::endl;
    }
}

TLorentzVector SmearedEvent::getExchangeBoson()
{
 
  TLorentzVector init(m_smearEvent->BeamLepton()->Get4Vector());
  TLorentzVector scat(m_smearEvent->ScatteredLepton()->Get4Vector());
  
  TLorentzVector exchangeBoson = init - scat;

  if(m_breitFrame)
    {
      BreitFrame breit(*m_truthEvent, *m_smearEvent);
      breit.labToBreitSmear(&exchangeBoson);
    }
 
  return exchangeBoson;

}
void SmearedEvent::setSmearedParticles()
{
  if (m_verbosity == -4)
    std::cout << "New D0 Smeared Event... " << std::endl;

  //std::cout << "Truth N Tracks " << m_truthEvent->GetNTracks() << " Smeared N Tracks " << m_smearEvent->GetNTracks() << std:: endl;

  int charmFlag = 0;
  double epsilon = 1e-7;
  BreitFrame breit(*m_truthEvent, *m_smearEvent);

  for(int part = 0; part < m_smearEvent->GetNTracks(); ++part)
    {
      
      /// Skip the beam
      if( part < 3 )
	continue;

      
      const Smear::ParticleMCS *particle = m_smearEvent->GetTrack(part);
      const Particle *truthParticle = m_truthEvent->GetTrack(part);
      
      /// only want final state particles -- should be taken care of by truth list, and will still clster D meson //
      // if(truthParticle->GetStatus() != 1)
      //	continue;


      /// only particles that could nominally be in the detector
      if(fabs(truthParticle->GetEta()) > m_maxPartEta)
	continue;
      if(truthParticle->GetPt() < m_minPartPt)
	continue;
      /// If truth particle wasn't smeared (e.g. out of acceptance), skip
      /// This can also manifest itself as E==0 && p==0
      if(particle == NULL || 
	 (fabs(particle->GetE()) <= epsilon && fabs(particle->GetP()) <= epsilon))
	continue;
      
      /// Skip the scattered electron, since it is special 
      /// (and we don't want it in jet finding, anyway)
      if(particle->GetE() == m_scatLepton->GetE())
	continue;
      
      double px = particle->GetPx();
      double py = particle->GetPy();
      double pz = particle->GetPz();
      double e  = particle->GetE();

      double truthPx = truthParticle->GetPx();
      double truthPy = truthParticle->GetPy();
      double truthPz = truthParticle->GetPz();
      double truthE  = truthParticle->GetE();


      if(m_verbosity > 2)
	{
	  std::cout << "Truth (lab) : "<<truthParticle->Id() 
		    << " " <<truthParticle->GetPx() << " " 
		    << truthParticle->GetPy() << " " << truthParticle->GetPz()
		    << " " << truthParticle->GetE() << std::endl;
	  
	  std::cout << "Smeared (lab) : " << px << " " 
		    << py << " " << pz << " " 
		    << e <<  std::endl;
	}

      /// We need to handle edge cases from EICsmear where e.g. a track
      /// is found but not connected to a calorimeter cluster
      double p = sqrt(px * px + py * py + pz * pz);
    
      /// Track found but not connected to calo cluster
      if(fabs(p) > epsilon && fabs(e) <= epsilon)
	{
	  /// The EIC framework ditches hadron/emcal info, so we need 
	  /// to cheat and find out what kind of particle it is
	  auto abspid = std::abs(truthParticle->GetPdgCode());
	  bool EM = abspid == 22 || abspid == 11;
	  if(EM)
	    {
	      /// must be an electron since there is a track
	      e = std::sqrt(p * p + 0.000511 * 0.000511);
	    } 
	  else 
	    {
	      /// Get truth mass
	      double m = truthParticle->GetM();

	      /// Calculate truth energy
	      e = std::sqrt(p * p + m * m);
	      
	      /// Smear energy out based on truth mass
	      if(m_smearHCal)
		{
		  float res = m_HCalRes->Eval(e);
		  /// Smear it by a gaussian with width of this resolution value
		  std::default_random_engine generator;
		  std::normal_distribution<double> distribution(0,res);
		  float smear = distribution(generator);
		  /// add the smearing term
		  e += (smear * e);
		}
	      else
		{
		  /// otherwise just assume it is a pion
		  e = std::sqrt(p * p + 0.139 * 0.139);
		}

	    }
	}
      
      /// Cal cluster only, no momentum info
      if( fabs(p) <= epsilon && fabs(e) > epsilon)
	{
	  /// again, cheat to figure out whether emcal or hcal
	  auto abspid = std::abs(truthParticle->GetPdgCode());
	  bool EM = abspid == 22 || abspid == 11;
	  if(EM){
	    /// assume m = 0 since electron mass is so small
	    p = e;	    
	  } 
	  else 
	    {
	      /// check truth mass
	      double m = truthParticle->GetM();
	      p = std::sqrt(e * e - m * m);
	      /// if particle was smeared such that e*e-m*m is negative, just set to e
	      if(p != p)
		p = e;
	    }

	  auto phi = particle->GetPhi();
	  auto theta = particle->GetTheta();
	  
	  px = p * sin(theta) * cos(phi);
	  py = p * sin(theta) * sin(phi);
	  pz = p * cos(theta);
	  
	}

      if(m_verbosity > 3)
	{
	  std::cout<<"particle to be boosted "<<std::endl;
	  std::cout<<"("<<px<<","<<py<<","<<pz<<","<<e<<")"<<std::endl;
	}


      // Check if smeared particle also passes some nominal pT cut
      if(sqrt(px * px + py * py) < m_minPartPt)
	continue;
	
      
      // This will skip the charm children for clustering and replace them with the charm hadron //
      bool chadChildren = false;
      for (int i = 0; i<m_chadChildIndices.size(); ++i)
      {
        if( truthParticle->GetIndex() == m_chadChildIndices.at(i) )
	  {
	    // if ( m_verbosity == -4 )
	    //  std::cout << "Child index .. " << m_chadChildIndices.at(i) << std::endl;

	    chadChildren = true;      
	  }
     
      }
      
      if ( chadChildren )
	{
	  const Smear::ParticleMCS *child1 = m_smearEvent->GetTrack( m_chadChildIndices.at(0) - 1 );
	  const Particle *truthChild1 = m_truthEvent->GetTrack( m_chadChildIndices.at(0) - 1 );

	  const Smear::ParticleMCS *child2 = m_smearEvent->GetTrack( m_chadChildIndices.at(1) - 1 );
	  const Particle *truthChild2 = m_truthEvent->GetTrack( m_chadChildIndices.at(1) - 1 );

	  charmFlag++;
	  if (charmFlag == 2)
	    {
	      TLorentzVector child1FourVec = child1->Get4Vector();
	      TLorentzVector truthChild1FourVec = truthChild1->PxPyPzE();
	      child1FourVec.SetE(std::sqrt(child1FourVec.P()*child1FourVec.P() + truthChild1FourVec.M2()));

	      TLorentzVector child2FourVec = child2->Get4Vector();
	      TLorentzVector truthChild2FourVec = truthChild2->PxPyPzE();;	      
	      child2FourVec.SetE(std::sqrt(child2FourVec.P()*child2FourVec.P() + truthChild2FourVec.M2()));
	      TLorentzVector *chadFourVec  = new TLorentzVector();
	      TLorentzVector *truthChadFourVec  = new TLorentzVector();

	      TLorentzVector temp = child1FourVec + child2FourVec;
	      TLorentzVector truthTemp = truthChild1FourVec + truthChild2FourVec;

	      chadFourVec->SetPxPyPzE( temp.Px(), temp.Py(), temp.Pz(), temp.E() );
	      truthChadFourVec->SetPxPyPzE( truthTemp.Px(), truthTemp.Py(), truthTemp.Pz(), truthTemp.E() );
	      /*
	      if ( m_verbosity == -4 )
		{
		  std::cout << "truth1 mass : " << truthChild1FourVec.M() << " truth2 mass : " << truthChild2FourVec.M() << " smear1 mass : " << child1FourVec.M() << " smear2 mass : " << child2FourVec.M() << std::endl;

		  std::cout << "truth D0 mass : " << truthChadFourVec->M() << " reco D0 mass : " << chadFourVec->M() << std::endl;
		}
	      */
	      if(m_breitFrame)
		{
		  breit.labToBreitTruth( truthChadFourVec );
		  breit.labToBreitSmear( chadFourVec );
		}
	      
	      m_particles.push_back(fastjet::PseudoJet(chadFourVec->Px(),
						       chadFourVec->Py(),
						       chadFourVec->Pz(),
						       chadFourVec->E()));
	      
	      m_truthParticles.push_back(fastjet::PseudoJet(truthChadFourVec->Px(),
	      						    truthChadFourVec->Py(),
							    truthChadFourVec->Pz(),
							    truthChadFourVec->E()));
	    }
	  
	  continue;
	}
      
      /*
      const Smear::ParticleMCS *testPart = m_smearEvent->GetTrack( part );
      const Particle *truthTestPart = m_truthEvent->GetTrack( part );
      TLorentzVector testPartVec = testPart->Get4Vector();
      TLorentzVector truthTestPartVec = truthTestPart->PxPyPzE();

      std::cout << "truth test part mass : " << truthTestPartVec.M() << " reco test part mass : " << testPartVec.M() << std::endl;
      */



      TLorentzVector *partFourVec = new TLorentzVector();
      TLorentzVector *truthPartFourVec = new TLorentzVector();
      partFourVec->SetPxPyPzE(px,py,pz,e);
      truthPartFourVec->SetPxPyPzE(truthPx,truthPy,truthPz,truthE);
      /*
      std::cout << "Truth (lab) : "<<truthParticle->Id() 
		    << " " <<truthPx << " " 
		    << truthPy << " " << truthPz
		<< " " << truthE << " " << truthPartFourVec->M() << std::endl;
	  
      std::cout << "Smeared (lab) : " << px << " " 
		    << py << " " << pz << " " 
		<< e << " " << partFourVec->M()<< std::endl;
      */
      if(m_breitFrame)
	{
	  breit.labToBreitTruth( truthPartFourVec );
	  breit.labToBreitSmear( partFourVec );
	}


      if(m_verbosity > 0)
	{
	  std::cout << "Smeared : " <<partFourVec->Px() << " " 
		    << partFourVec->Py() << " " << partFourVec->Pz()
		    << " " << partFourVec->E() << " mass : " << partFourVec->M() << std::endl;	      
	}

      m_particles.push_back(fastjet::PseudoJet(partFourVec->Px(),
					       partFourVec->Py(),
					       partFourVec->Pz(),
					       partFourVec->E()));
      /// vectors will have a one-to-one correspondance
      m_truthParticles.push_back(fastjet::PseudoJet(truthPartFourVec->Px(),
						    truthPartFourVec->Py(),
						    truthPartFourVec->Pz(),
						    truthPartFourVec->E()));
    }

  return;
}

void SmearedEvent::setSmearedParticlesPid()
{
  std::vector<int> skipParts;
  if (m_verbosity == -4)
    std::cout << "New D0 Smeared Event... " << std::endl;

  //std::cout << "Truth N Tracks " << m_truthEvent->GetNTracks() << " Smeared N Tracks " << m_smearEvent->GetNTracks() << std:: endl;

  int charmFlag = 0;
  double epsilon = 1e-7;
  BreitFrame breit(*m_truthEvent, *m_smearEvent);

  for(int part = 0; part < m_smearEvent->GetNTracks(); ++part)
    {
      
      /// Skip the beam
      if( part < 3 )
	{
	  skipParts.push_back(part);
	  std:: cout << part << std::endl;
	  continue;
	}

      
      const Smear::ParticleMCS *particle = m_smearEvent->GetTrack(part);
      const Particle *truthParticle = m_truthEvent->GetTrack(part);
      
      /// only want final state particles -- should be taken care of by truth list, and will still clster D meson //
      if(truthParticle->GetStatus() != 1)
	{
	  skipParts.push_back(part);
	  std:: cout << part << std::endl;
	  continue;
	}

      /// only particles that could nominally be in the detector
      if(fabs(truthParticle->GetEta()) > m_maxPartEta)
	{
	  skipParts.push_back(part);
	  std:: cout << part << std::endl;
	  continue;
	}
      if(truthParticle->GetPt() < m_minPartPt)
	{
	  skipParts.push_back(part);
	  std:: cout << part << std::endl;
	  continue;
	}

      /// If truth particle wasn't smeared (e.g. out of acceptance), skip
      /// This can also manifest itself as E==0 && p==0
      if(particle == NULL || 
	 (fabs(particle->GetE()) <= epsilon && fabs(particle->GetP()) <= epsilon))
	{
	  skipParts.push_back(part);
	  std:: cout << part << std::endl;
	  continue;
	}
      /// Skip the scattered electron, since it is special 
      /// (and we don't want it in jet finding, anyway)
      if(particle->GetE() == m_scatLepton->GetE())
	{
	  skipParts.push_back(part);
	  std:: cout << part << std::endl;
	  continue;
	}
      
      double px = particle->GetPx();
      double py = particle->GetPy();
      double pz = particle->GetPz();
      double e  = particle->GetE();

      double truthPx = truthParticle->GetPx();
      double truthPy = truthParticle->GetPy();
      double truthPz = truthParticle->GetPz();
      double truthE  = truthParticle->GetE();

      if(m_verbosity > 2)
	{
	  std::cout << "Truth (lab) : "<<truthParticle->Id() 
		    << " " <<truthParticle->GetPx() << " " 
		    << truthParticle->GetPy() << " " << truthParticle->GetPz()
		    << " " << truthParticle->GetE() << std::endl;
	  
	  std::cout << "Smeared (lab) : " << px << " " 
		    << py << " " << pz << " " 
		    << e <<  std::endl;
	}

      /// We need to handle edge cases from EICsmear where e.g. a track
      /// is found but not connected to a calorimeter cluster
      double p = sqrt(px * px + py * py + pz * pz);
    
      /// Track found but not connected to calo cluster
      if(fabs(p) > epsilon && fabs(e) <= epsilon)
	{
	  /// The EIC framework ditches hadron/emcal info, so we need 
	  /// to cheat and find out what kind of particle it is
	  auto abspid = std::abs(truthParticle->GetPdgCode());
	  bool EM = abspid == 22 || abspid == 11;
	  if(EM)
	    {
	      /// must be an electron since there is a track
	      e = std::sqrt(p * p + 0.000511 * 0.000511);
	    } 
	  else 
	    {
	      /// Get truth mass
	      double m = truthParticle->GetM();

	      /// Calculate truth energy
	      e = std::sqrt(p * p + m * m);
	      
	      /// Smear energy out based on truth mass
	      if(m_smearHCal)
		{
		  float res = m_HCalRes->Eval(e);
		  /// Smear it by a gaussian with width of this resolution value
		  std::default_random_engine generator;
		  std::normal_distribution<double> distribution(0,res);
		  float smear = distribution(generator);
		  /// add the smearing term
		  e += (smear * e);
		}
	      else
		{
		  /// otherwise just assume it is a pion
		  e = std::sqrt(p * p + 0.139 * 0.139);
		}

	    }
	}
      
      /// Cal cluster only, no momentum info
      if( fabs(p) <= epsilon && fabs(e) > epsilon)
	{
	  /// again, cheat to figure out whether emcal or hcal
	  auto abspid = std::abs(truthParticle->GetPdgCode());
	  bool EM = abspid == 22 || abspid == 11;
	  if(EM){
	    /// assume m = 0 since electron mass is so small
	    p = e;	    
	  } 
	  else 
	    {
	      /// check truth mass
	      double m = truthParticle->GetM();
	      p = std::sqrt(e * e - m * m);
	      /// if particle was smeared such that e*e-m*m is negative, just set to e
	      if(p != p)
		p = e;
	    }

	  auto phi = particle->GetPhi();
	  auto theta = particle->GetTheta();
	  
	  px = p * sin(theta) * cos(phi);
	  py = p * sin(theta) * sin(phi);
	  pz = p * cos(theta);
	  
	}

      if(m_verbosity > 3)
	{
	  std::cout<<"particle to be boosted "<<std::endl;
	  std::cout<<"("<<px<<","<<py<<","<<pz<<","<<e<<")"<<std::endl;
	}


      // Check if smeared particle also passes some nominal pT cut
      if(sqrt(px * px + py * py) < m_minPartPt)
	{
	  skipParts.push_back(part);
	  std:: cout << part << std::endl;
	  continue;
	}

      double eta = particle->GetEta();
      
      // Let's implement PID here... //
      int truthPid = truthParticle->GetPdgCode();
      int pid = 0;

      /// Get truth mass
      double mass = truthParticle->GetM();

      pid = pidDetectors(eta, p, mass, e, truthPid);
      /*
      if (defaultPID)
	{
	  if ( abs(truthPid) == 211 || abs(truthPid) == 321 || abs(truthPid) == 2212 )
	    {
	      if (-3.5 < eta && eta < -1.0)
		if (p < 7)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (-1.0 < eta && eta < 1.0)
		if (p < 5)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (1.0 < eta && eta < 2.0)
		if (p < 8)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (2.0 < eta && eta < 3.0)
		if (p < 20)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (3.0 < eta && eta < 3.5)
		if (p < 45)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	    }
	}

      else if (requestedPID)
	{
	  if ( abs(truthPid) == 211 || abs(truthPid) == 321 || abs(truthPid) == 2212 )
	    {
	      if (-3.5 < eta && eta < -1.0)
		if (p < 7.0)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (-1.0 < eta && eta < 0.0)
		if (p < 10.0)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (0.0 < eta && eta < 1.0)
		if (p < 15.0)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (1.0 < eta && eta < 1.5)
		if (p < 30.0)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (1.5 < eta && eta < 2.0)
		if (p < 50.0)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (2.0 < eta && eta < 2.5)
		if (p < 50.0)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (2.5 < eta && eta < 3.0)
		if (p < 30.0)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	      if (3.0 < eta && eta < 3.5)
		if (p < 20.0)
		  {
		    pid = truthPid;
		    e = std::sqrt(p * p + m * m);
		  }
	    }
	}
      */
      std::cout << "Smeared particle PID : " << pid << std::endl;
      std::cout << "Truth particle PID : " << truthPid << std::endl;

      TLorentzVector *partFourVec = new TLorentzVector();
      TLorentzVector *truthPartFourVec = new TLorentzVector();
      partFourVec->SetPxPyPzE(px,py,pz,e);
      truthPartFourVec->SetPxPyPzE(truthPx,truthPy,truthPz,truthE);


      /*
      std::cout << "Truth (lab) : "<<truthParticle->Id() 
		    << " " <<truthPx << " " 
		    << truthPy << " " << truthPz
		<< " " << truthE << " " << truthPartFourVec->M() << std::endl;
	  
      std::cout << "Smeared (lab) : " << px << " " 
		    << py << " " << pz << " " 
		<< e << " " << partFourVec->M()<< std::endl;
      */

      TLorentzVector partVec(px,py,pz,e);
      TLorentzVector pairVec;

      //std::cout<< m_smearEvent->GetTrack(23)->GetPx() << std::endl;

      
      for(int part2 = 0; part2 < m_smearEvent->GetNTracks(); ++part2)
	{
	  // if  (std::find(skipParts.begin(), skipParts.end(), part2) != skipParts.end())
	  //continue;

	  /// Skip the beam
	  if( part2 < 3 )
	    {
	      continue;
	    }
    
	  const Smear::ParticleMCS *particle2 = m_smearEvent->GetTrack(part2);
	  const Particle *truthParticle2 = m_truthEvent->GetTrack(part2);
	  
	  if(truthParticle2->GetStatus() != 1)
	    {
	      continue;
	    }
	  
	  /// only particles that could nominally be in the detector
	  if(fabs(truthParticle2->GetEta()) > m_maxPartEta)
	    {
	      continue;
	    }
	  if(truthParticle2->GetPt() < m_minPartPt)
	    {
	      continue;
	    }
	  
	  /// If truth particle wasn't smeared (e.g. out of acceptance), skip
	  /// This can also manifest itself as E==0 && p==0
	  if(particle2 == NULL || 
	     (fabs(particle2->GetE()) <= epsilon && fabs(particle2->GetP()) <= epsilon))
	    {
	      continue;
	    }
	  /// Skip the scattered electron, since it is special 
	  /// (and we don't want it in jet finding, anyway)
	  if(particle2->GetE() == m_scatLepton->GetE())
	    {
	      continue;
	    }

	  double px2 = particle2->GetPx();
	  double py2 = particle2->GetPy();
	  double pz2 = particle2->GetPz();
	  double e2  = particle2->GetE();

	  /// We need to handle edge cases from EICsmear where e.g. a track
	  /// is found but not connected to a calorimeter cluster
	  double p2 = sqrt(px2 * px2 + py2 * py2 + pz2 * pz2);
	  /// Track found but not connected to calo cluster
	  if(fabs(p2) > epsilon && fabs(e2) <= epsilon)
	    {
	      /// The EIC framework ditches hadron/emcal info, so we need 
	      /// to cheat and find out what kind of particle it is
	      auto abspid2 = std::abs(truthParticle2->GetPdgCode());
	      bool EM2 = abspid2 == 22 || abspid2 == 11;
	      if(EM2)
		{
		  /// must be an electron since there is a track
		  e2 = std::sqrt(p2 * p2 + 0.000511 * 0.000511);
		} 
	      else 
		{
		  /// Get truth mass
		  double m2 = truthParticle2->GetM();
		  
		  /// Calculate truth energy
		  e2 = std::sqrt(p2 * p2 + m2 * m2);
		  
		  /// Smear energy out based on truth mass
		  if(m_smearHCal)
		    {
		      float res = m_HCalRes->Eval(e2);
		      /// Smear it by a gaussian with width of this resolution value
		      std::default_random_engine generator;
		      std::normal_distribution<double> distribution(0,res);
		      float smear = distribution(generator);
		      /// add the smearing term
		      e2 += (smear * e2);
		    }
		  else
		    {
		      /// otherwise just assume it is a pion
		      e2 = std::sqrt(p2 * p2 + 0.139 * 0.139);
		    }
		  
		}
	    }
	  
	  /// Cal cluster only, no momentum info
	  if( fabs(p2) <= epsilon && fabs(e2) > epsilon)
	    {
	      /// again, cheat to figure out whether emcal or hcal
	      auto abspid2 = std::abs(truthParticle2->GetPdgCode());
	      bool EM2 = abspid2 == 22 || abspid2 == 11;
	      if(EM2){
		/// assume m = 0 since electron mass is so small
		p2 = e2;	    
	      } 
	      else 
		{
		  /// check truth mass
		  double m2 = truthParticle2->GetM();
		  p2 = std::sqrt(e2 * e2 - m2 * m2);
		  /// if particle was smeared such that e*e-m*m is negative, just set to e
		  if(p2 != p2)
		    p2 = e2;
		}
	      
	      auto phi2 = particle2->GetPhi();
	      auto theta2 = particle2->GetTheta();
	      
	      px2 = p2 * sin(theta2) * cos(phi2);
	      py2 = p2 * sin(theta2) * sin(phi2);
	      pz2 = p2 * cos(theta2);
	      
	    }

	  // Check if smeared particle also passes some nominal pT cut
	  if(sqrt(px2 * px2 + py2 * py2) < m_minPartPt)
	    {
	      continue;
	    }
	  
	  //const Smear::ParticleMCS *particle2 = m_smearEvent->GetTrack(part2);
	  //std::cout << "second loop counter " << part2 << std::endl;
	  // double 


	  TLorentzVector part2Vec(px2, py2, pz2, e2);
	  pairVec = partVec + part2Vec;
	  std::cout << "Pair Mass : " << pairVec.M() << std::endl;
	  	  
	}

      if(m_breitFrame)
	{
	  breit.labToBreitTruth( truthPartFourVec );
	  breit.labToBreitSmear( partFourVec );
	}


      if(m_verbosity > 0)
	{
	  std::cout << "Smeared : " <<partFourVec->Px() << " " 
		    << partFourVec->Py() << " " << partFourVec->Pz()
		    << " " << partFourVec->E() << " mass : " << partFourVec->M() << std::endl;	      
	}

     m_particles.push_back(fastjet::PseudoJet(partFourVec->Px(),
					       partFourVec->Py(),
					       partFourVec->Pz(),
					       partFourVec->E()));
      /// vectors will have a one-to-one correspondance
      m_truthParticles.push_back(fastjet::PseudoJet(truthPartFourVec->Px(),
						    truthPartFourVec->Py(),
						    truthPartFourVec->Pz(),
						    truthPartFourVec->E()));
    }

  return;
}

PseudoJetVec SmearedEvent::getRecoJets(fastjet::ClusterSequence *cs, 
				       JetDef jetDef)
{
  PseudoJetVec allRecoJets;

  /// Check that there weren't any weird EICSmear edge cases
  for(int i = 0; i< m_particles.size(); i++)
    {
      double px = m_particles.at(i).px();
      double py = m_particles.at(i).py();
      double pz = m_particles.at(i).pz();
      double e = m_particles.at(i).e();
      /// If any of the particles have nan components, the clustering will fail
      /// so we don't want this event anyway because it indicates an edge 
      /// case in EICsmear, so just return an empty vector
      if(px != px || py != py || pz != pz || e!=e)
	return allRecoJets;
    }

  /// Create the cluster sequence
  cs = new fastjet::ClusterSequence(m_particles, jetDef.getJetDef());
  
  allRecoJets = fastjet::sorted_by_pt(cs->inclusive_jets());

  if(m_verbosity > 1)
    {
      std::cout << "Finding jets for jet pT>" << jetDef.getMinJetPt() 
		<< " and |eta|<" << jetDef.getMaxJetRapidity() << std::endl;
    }
  
  /// Make eta/pt selections
  fastjet::Selector selectPt = fastjet::SelectorPtMin(jetDef.getMinJetPt());
  fastjet::Selector selectEta = fastjet::SelectorAbsRapMax(jetDef.getMaxJetRapidity());
  fastjet::Selector select = selectPt and selectEta;
  
  PseudoJetVec selectJets = select(allRecoJets);

  PseudoJetVec charmJets = CharmJetTagging(selectJets);

  // return selectJets;
  return charmJets;

}

TLorentzPairVec SmearedEvent::getMatchedParticles()
{
  TLorentzPairVec pairVec;
  for(int i = 0; i < m_particles.size(); i++)
    {
      fastjet::PseudoJet truth = m_truthParticles.at(i);
      fastjet::PseudoJet reco = m_particles.at(i);
      TLorentzVector tlTruth, tlReco;
      tlTruth.SetPxPyPzE(truth.px(), truth.py(), truth.pz(), truth.e());
      tlReco.SetPxPyPzE(reco.px(), reco.py(), reco.pz(), reco.e());

      std::pair<TLorentzVector, TLorentzVector> pair = 
	std::make_pair(tlTruth, tlReco);

      pairVec.push_back(pair);
    }

  return pairVec;
}

PseudoJetVec SmearedEvent::getRecoSoftDropJets(PseudoJetVec recoJets, SoftDropJetDef sdJetDef)
{

  fastjet::contrib::SoftDrop sd(sdJetDef.getSoftDrop());
  
  PseudoJetVec softDropJets;
 
  for(int jet = 0; jet < recoJets.size(); ++jet)
    {
      fastjet::PseudoJet softDropJet = sd(recoJets[jet]);
      softDropJets.push_back(softDropJet);
    }

  return softDropJets;

}


std::vector<PseudoJetVec> SmearedEvent::matchTruthRecoJets(
	    PseudoJetVec truthjets,
	    PseudoJetVec recojets)
{
  PseudoJetVec matchJetPair;
  std::vector<PseudoJetVec> matchedVector;
  
  for(int reco = 0; reco < recojets.size(); ++reco)
    {
      matchJetPair.clear();
      
      fastjet::PseudoJet recojet = recojets[reco];
      double mindR = 9999;
      fastjet::PseudoJet matchedTruthJet;
      
      /// Find the closest truth jet in DeltaR space
      for(int truth = 0; truth < truthjets.size(); ++truth)
	{
	  fastjet::PseudoJet truthjet = truthjets[truth];
	  double dR = recojet.delta_R(truthjet);
	  if(dR < mindR)
	    {
	      mindR = dR;
	      matchedTruthJet = truthjet;
	      
	      if(m_verbosity > 1)
		std::cout << "Found matched jet with dr " << mindR<<std::endl;
	    }
	}
      
      if(m_verbosity > 1)
	{
	  std::cout << "recojet: " << recojet.px() << "  " << recojet.py() 
		    << "  " << recojet.pz() << "  " << recojet.e() << std::endl;
	  std::cout << "truejet: " << matchedTruthJet.px() << "  " 
		    << matchedTruthJet.py() << "  " << matchedTruthJet.pz() 
		    << "  " << matchedTruthJet.e() << std::endl;
	}
      
      matchJetPair.push_back(matchedTruthJet);
      matchJetPair.push_back(recojet);
      matchedVector.push_back(matchJetPair);
    }

  return matchedVector;

}

/*
bool SmearedEvent::D0kpiNoSmearFilter()
{
  double epsilon = 1e-7;
  const Smear::ParticleMCS *child1 = m_smearEvent->GetTrack(m_chadChildIndices.at(0)-1);
  const Smear::ParticleMCS *child2 = m_smearEvent->GetTrack(m_chadChildIndices.at(1)-1);
  const Particle *truthChild1 = m_truthEvent->GetTrack(m_chadChildIndices.at(0)-1);
  const Particle *truthChild2 = m_truthEvent->GetTrack(m_chadChildIndices.at(1)-1);
  // if (  child1 == NULL ||  fabs(child1->GetE()) <= epsilon || fabs(child1->GetP()) <= epsilon || 
  //	child2 == NULL ||  fabs(child2->GetE()) <= epsilon || fabs(child2->GetP()) <= epsilon  )
  //  return false;
  //else
  //return true;

  if (  child1 == NULL || child2 == NULL  )
    return false;
  else
    return true;
}
*/

PseudoJetVec SmearedEvent::CharmJetTagging(PseudoJetVec recoJets)
{
  PseudoJetVec charmJets;
  std::vector<fastjet::PseudoJet> cons;

  for(int jet = 0; jet < recoJets.size(); ++jet)
    {
      cons.clear();
      cons = recoJets.at(jet).constituents();
      if ( m_verbosity == -4 ) 
	{
	  std::cout << "new jet : ";
	  std::cout << "num constituents : " << cons.size() << std::endl;
	}
      for (int con = 0; con < cons.size(); ++con)
	{
	  fastjet::PseudoJet constituent;
	  constituent = cons.at(con);
	  TLorentzVector const_vect(constituent.px(), constituent.py(), constituent.pz(), constituent.e());

	  if (fabs(const_vect.M() - 1.864) < 0.03 )
	    {
	      
	      if ( m_verbosity == -4 ) 
		{
		  std::cout << " we have a reco level charm tagged jet! " << std::endl;
		}
	      charmJets.push_back(recoJets.at(jet));
	    }      
	}
            
      // Lets just tag D0s since we can readily reconstruct these

    }

  return charmJets;
}


int SmearedEvent::pidDetectors(double eta, double p, double m, double &e, int truthPid )
{
  int pid = 0;
  if (defaultPID)
    {
      if ( abs(truthPid) == 211 || abs(truthPid) == 321 || abs(truthPid) == 2212 )
	{
	  if (-3.5 < eta && eta < -1.0)
	    if (p < 7)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (-1.0 < eta && eta < 1.0)
	    if (p < 5)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (1.0 < eta && eta < 2.0)
	    if (p < 8)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (2.0 < eta && eta < 3.0)
	    if (p < 20)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (3.0 < eta && eta < 3.5)
	    if (p < 45)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	}
    }
  
  else if (requestedPID)
    {
      if ( abs(truthPid) == 211 || abs(truthPid) == 321 || abs(truthPid) == 2212 )
	{
	  if (-3.5 < eta && eta < -1.0)
	    if (p < 7.0)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (-1.0 < eta && eta < 0.0)
	    if (p < 10.0)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (0.0 < eta && eta < 1.0)
	    if (p < 15.0)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (1.0 < eta && eta < 1.5)
	    if (p < 30.0)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (1.5 < eta && eta < 2.0)
	    if (p < 50.0)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (2.0 < eta && eta < 2.5)
	    if (p < 50.0)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (2.5 < eta && eta < 3.0)
	    if (p < 30.0)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	  if (3.0 < eta && eta < 3.5)
	    if (p < 20.0)
	      {
		pid = truthPid;
		e = std::sqrt(p * p + m * m);
	      }
	}
    }
  return pid;
}
