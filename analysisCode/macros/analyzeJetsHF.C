#include "analyzeJetsHF.h"
#include "HistoManagerHF.h"
const float pion_mass = 0.13957;
const float kaon_mass = 0.4936;
int verbosity = 0;


void analyzeJetsHF(std::string file)
{
  std::string filename = file;
  infile = TFile::Open(filename.c_str());

  gROOT->ProcessLine(".L ../src/fastJetLinker.C+");
  setupTree();

  instantiateHistos();

  loop();

  write(file);
}


void recoJetAnalysis(JetConstVec *recojets)
{

  for(int jet = 0; jet < recojets->size(); jet++)
    {
      TLorentzVector jetVec;
      jetVec = recojets->at(jet).first;

      float jetpt = jetVec.Pt();
      //if(jetpt < minjetpt )
      //continue;
      //if(fabs(jetVec.Eta()) > maxjeteta)
      //continue;

      recojetpteta->Fill(jetpt, jetVec.Eta());
      recojetptphi->Fill(jetpt, jetVec.Phi());
      recojetmass->Fill( jetVec.M() );
      recojetnconstpt->Fill( recojets->at(jet).second.size(), jetpt );


      TVector3 jet3;
      jet3.SetXYZ(jetVec.Px(), jetVec.Py(), jetVec.Pz());

      /// Iterate over constituents
      for(int j = 0; j < recojets->at(jet).second.size(); ++j)
	{
	  TLorentzVector con = recojets->at(jet).second.at(j);
	  TVector3 con3;
	  con3.SetXYZ(con.Px(), con.Py(), con.Pz());
	  TVector3 cross = jet3.Cross(con3);
	  
	  float z = jet3.Dot(con3) / (jet3.Mag2());
	  float jt = cross.Mag() / jet3.Mag();
	  float r = sqrt(pow(checkdPhi(jetVec.Phi() - con.Phi()), 2) + pow(jetVec.Eta() - con.Eta(),2));
	  
	  if (con.M() > 1.5 && con.M() < 2.0)
	    {
	      recojetptz->Fill(z, jetpt);
	      recojetptjt->Fill(jt, jetpt);
	      recojetptr->Fill(r, jetpt);
	    }

	  recoconstmass->Fill(con.M());
	  
	}
    }
}

double truthJetAnalysis(JetConstVec *truthjets)
{

  double jetpt = 0;
  for(int jet = 0; jet < truthJets->size(); jet++)
    {
      TLorentzVector jetVec;
      jetVec = truthJets->at(jet).first;
      if(jetVec.Pt() > jetpt)
	jetpt = jetVec.Pt();
       
      // if(jetVec.Pt() < minjetpt)
      //continue;
      //if(fabs(jetVec.Eta()) > maxjeteta)
      //continue;
      
      truejetptphi->Fill(jetVec.Pt(), jetVec.Phi());
      truejetpteta->Fill(jetVec.Pt(), jetVec.Eta());
      truejetmass->Fill( jetVec.M() );
      truejetnconstpt->Fill( truthjets->at(jet).second.size(), jetVec.Pt() );


      TVector3 jet3;
      jet3.SetXYZ(jetVec.Px(), jetVec.Py(), jetVec.Pz());

      TLorentzVector con;
      for(int j = 0; j < truthJets->at(jet).second.size(); ++j)
	{
	  con = truthJets->at(jet).second.at(j);
	  TVector3 con3;
	  con3.SetXYZ(con.Px(), con.Py(), con.Pz());
	  TVector3 cross = jet3.Cross(con3);	  

	  float z = jet3.Dot(con3) / (jet3.Mag2());
	  float jt = cross.Mag() / jet3.Mag();
	  float r = sqrt(pow(checkdPhi(jetVec.Phi() - con.Phi()), 2) + pow(jetVec.Eta() - con.Eta(),2));

	  if (con.M() > 1.5 && con.M() < 2.0)
	    {
	      truejetptz->Fill(z, jetpt);
	      truejetptjt->Fill(jt, jetpt);
	      truejetptr->Fill(r, jetpt);
	    }

	  trueconstmass->Fill(con.M());
	  
	} 
    }  

  return jetpt;
}

void loop()
{

  for(int i=0; i<jettree->GetEntries(); i++)
    {
      if (verbosity == 1)
	cout << "New Event... " << endl;

      if(i%10000 == 0)
	std::cout << "Processed " << i << " events " << std::endl;
      jettree->GetEntry(i);
      
      recoJetAnalysis(recoJets);

      float highestTruthJetPt = truthJetAnalysis(truthJets);

      analyzeMatchedJets(matchedJets, matchedParticles);

      recoSDJetAnalysis(recoSDJets);
      truthSDJetAnalysis(truthSDJets);
      analyzeMatchedSDJets(matchedSDJets);

      compareAKTSDTruthJets(truthJets, truthSDJets);

      /// Event level kinematics
      truerecx->Fill(truex,recx);
      truerecy->Fill(truey,recy);
      truerecq2->Fill(trueq2,recq2);
      trueQ2x->Fill(truex,trueq2);
      trueQ2pT->Fill(trueq2, highestTruthJetPt);
      truenjetevent->Fill(truthJets->size());
      reconjetevent->Fill(recoJets->size());
  

    }

}

void compareAKTSDTruthJets(JetConstVec *truthjets, JetConstVec *truthsdjets)
{
  for(int i=0; i<truthjets->size(); ++i)
    {

      TLorentzVector aktjet = truthjets->at(i).first;
      TLorentzVector sdjet = truthsdjets->at(i).first;
    
      sdenergygroomed->Fill(sdjet.E() / aktjet.E());
    }

}

void analyzeMatchedJets(MatchedJets *matchedjets,
			TLorentzPairVec *matchedparticles)
{
  for(int i = 0; i < matchedjets->size(); i++)
    {
      JetConstPair truthJetConst = matchedjets->at(i).at(0);
      JetConstPair recoJetConst = matchedJets->at(i).at(1);

      TLorentzVector truthJet = truthJetConst.first;
      TLorentzVector recoJet = recoJetConst.first;
      TLorentzVectorVec truthConst = truthJetConst.second;
      TLorentzVectorVec recoConst = recoJetConst.second;

      //  if(truthJet.Pt() < minjetpt || fabs(truthJet.Eta()) > maxjeteta)
      //continue;
      if(recoJet.Pt() >minjetpt && fabs(recoJet.Eta()) < maxjeteta)
	{
	  recojetptetatruejetpt->Fill(truthJet.Pt(), truthJet.Eta());
	}
      matchedJetDr->Fill((float)truthJet.DeltaR(recoJet));
      matchedJetdPhi->Fill((float)truthJet.DeltaPhi(recoJet));
      matchedJetdEta->Fill((float)truthJet.Eta() - recoJet.Eta());

      recotruejetpt->Fill(truthJet.Pt(), recoJet.Pt());
      recotruejeteta->Fill(truthJet.Eta(), recoJet.Eta());
      recotruejetphi->Fill(truthJet.Phi(), recoJet.Phi());
      
      recotruejetp->Fill(truthJet.P(), recoJet.P());
      recotruejete->Fill(truthJet.E(), recoJet.E());

      /// Match constituents up
      for(int j = 0; j< recoConst.size(); j++)
	{
	  TLorentzVector recoCon = recoConst.at(j);
	  TLorentzVector truthMatch;
	  
	  for(int k =0; k< matchedparticles->size(); k++)
	    {
	      TLorentzVector matchreco = matchedparticles->at(k).second;
	      if(matchreco.Px() == recoCon.Px() &&
		 matchreco.Py() == recoCon.Py() &&
		 matchreco.Pz() == recoCon.Pz())
		{
		  truthMatch = matchedparticles->at(k).first;
		  if (verbosity == 1)
		    if (truthMatch.M() > 1.5 && truthMatch.M() < 2.0)
		      cout << " reco D0 mass : " <<  matchreco.M() << " truth D0 mass " << truthMatch.M() << endl;

		}
	    }

	  bool matched = false;
	  /// now check that the truth particle was actually in the jet
	  for(int k = 0; k < truthConst.size(); k++)
	    {
	      TLorentzVector truthCon = truthConst.at(k);
	      // if(truthCon.Px() == truthMatch.Px() &&
	      // truthCon.Py() == truthMatch.Py() &&
	      // truthCon.Pz() == truthMatch.Pz())

	      if(fabs(truthCon.Px() -  truthMatch.Px()) < 0.00001 &&
		 fabs(truthCon.Py() -  truthMatch.Py()) < 0.00001 &&
		 fabs(truthCon.Pz() -  truthMatch.Pz()) < 0.00001)
		{

		if (truthMatch.M() > 1.5 && truthMatch.M() < 2.0)
		  {
		    if (verbosity == 1)
		      cout << " Found a matched D0 " << endl;
		  }
		  matched = true;
		  break;
		}
	      else
		if (verbosity == 1)
		  if (truthMatch.M() > 1.5 && truthMatch.M() < 2.0)
		    {
		      cout << " truthMatchP : (" << truthMatch.Px() <<" ,"<< truthMatch.Py() << " ," << truthMatch.Pz() <<" )" << " truthMatchM : " << truthMatch.M() << endl;
		      cout << " truthConP : (" << truthCon.Px() <<" ,"<< truthCon.Py() << " ," << truthCon.Pz() <<" )" << " truthConM : " << truthCon.M() <<  endl;

		    }
	    }
	  if(matched)
	    {
	      /// Found a matched truth constituent and it was in the truth jet
	      ///recoCon and truthMatch
	      ////truthJet and recoJet
	      TVector3 truthJet3, recoJet3, recoCon3, truthMatch3;
	      truthJet3.SetXYZ(truthJet.Px(), 
			       truthJet.Py(), truthJet.Pz());
	      recoJet3.SetXYZ(recoJet.Px(), 
			      recoJet.Py(), recoJet.Pz());
	      recoCon3.SetXYZ(recoCon.Px(), 
			      recoCon.Py(), recoCon.Pz());
	      truthMatch3.SetXYZ(truthMatch.Px(),
				 truthMatch.Py(), truthMatch.Pz());

	      if (truthMatch.M() > 1.5 && truthMatch.M() < 2.0)
		{
		  float recoz = recoJet3.Dot(recoCon3) / (recoJet3.Mag2());
		  float truthz = truthJet3.Dot(truthMatch3) / (truthJet3.Mag2());
		  TVector3 truecross = truthJet3.Cross(truthMatch3);
		  TVector3 recocross = recoJet3.Cross(recoCon3);
		  float recojt = recocross.Mag() / recoJet3.Mag();
		  float truejt = truecross.Mag() / truthJet3.Mag();
		  float recodphi = checkdPhi(recoJet.Phi() - recoCon.Phi());
		  float truedphi = checkdPhi(truthJet.Phi() - truthMatch.Phi());
		  
		  truthRecoConstdPhi->Fill(checkdPhi(truthMatch.Phi() - recoCon.Phi()));
		  truthRecoConstdEta->Fill(truthMatch.Eta() - recoCon.Eta());
		  truthRecoConstdRap->Fill(truthMatch.Rapidity() - recoCon.Rapidity());
		  
		  float recor = sqrt(pow(recodphi ,2) +
				     pow(recoJet.Rapidity() - recoCon.Rapidity(), 2));
		  float truer = sqrt(pow(truedphi ,2) +
				     pow(truthJet.Rapidity() - truthMatch.Rapidity(),2));
	      
		  truerecoz->Fill(truthz, recoz);
		  truerecojt->Fill(truejt,recojt);
		  truerecor->Fill(truer, recor);

		  if (verbosity == 1)
		    cout << " D0 info written to histo " << endl;

		}
	      else
		{
		}
	    }
	  else
	    {
	      /// If a match couldn't be found, reco jet const was mistakenly
	      /// reconstructed within jet
	    }

	}
    }

  // cout << "number of unmatched jets... " << nomatch << endl
}

void recoSDJetAnalysis(JetConstVec *recojets)
{
  for(int ijet = 0; ijet< recojets->size(); ijet++)
    {
      TLorentzVector jet = recojets->at(ijet).first;
      TLorentzVectorVec constituents = recojets->at(ijet).second;
      ///first two constituents are the subjets
      TLorentzVector subjet1, subjet2;
      subjet1 = constituents.at(0);
      subjet2 = constituents.at(1);
      float zg = std::min(subjet1.Pt(), subjet2.Pt()) / (subjet1.Pt() + subjet2.Pt());
      float Rg = subjet1.DeltaR(subjet2);
      recoSDjetzg->Fill(zg, jet.Pt());
      recoSDjetrg->Fill(Rg, jet.Pt());
    }

}

void truthSDJetAnalysis(JetConstVec *truthjets)
{
  for(int ijet = 0; ijet < truthjets->size(); ijet++)
    {
      TLorentzVector jet;
      jet = truthjets->at(ijet).first;
      
      TLorentzVectorVec constituents;
      constituents = truthjets->at(ijet).second;
      ///first two constituents are the subjets
      TLorentzVector subjet1, subjet2;
      subjet1 = constituents.at(0);
      subjet2 = constituents.at(1);
      float zg = std::min(subjet1.Pt(), subjet2.Pt()) / (subjet1.Pt() + subjet2.Pt());
      float Rg = subjet1.DeltaR(subjet2);
      truthSDjetzg->Fill(zg, jet.Pt());
      truthSDjetrg->Fill(Rg, jet.Pt());
    }

}

void analyzeMatchedSDJets(MatchedJets *matchedjets)
{
  for(int i = 0; i < matchedjets->size(); i++)
    {
      JetConstPair truthJetConst = matchedjets->at(i).at(0);
      JetConstPair recoJetConst = matchedjets->at(i).at(1);
      
      TLorentzVector truthJet = truthJetConst.first;
      TLorentzVector recoJet = recoJetConst.first;
      TLorentzVectorVec truthConst = truthJetConst.second;
      TLorentzVectorVec recoConst = recoJetConst.second;

      TLorentzVector truthSubjet1, truthSubjet2;
      truthSubjet1 = truthConst.at(0);
      truthSubjet2 = truthConst.at(1);

      TLorentzVector recoSubjet1, recoSubjet2;
      recoSubjet1 = recoConst.at(0);
      recoSubjet2 = recoConst.at(1);
      
      float truthzg = std::min(truthSubjet1.Pt(), truthSubjet2.Pt())/(truthSubjet1.Pt() + truthSubjet2.Pt());
      float truthrg = truthSubjet1.DeltaR(truthSubjet2);
      float recozg = std::min(recoSubjet1.Pt(), recoSubjet2.Pt())/(recoSubjet1.Pt() + recoSubjet2.Pt());
      float recorg = recoSubjet1.DeltaR(recoSubjet2);
      truthrecozg->Fill(truthzg,recozg);
      truthrecorg->Fill(truthrg,recorg);

    }


}

void setupTree()
{

  jettree = (TTree*)infile->Get("jettree");
  jettree->SetBranchAddress("recx", &recx);
  jettree->SetBranchAddress("recq2", &recq2);
  jettree->SetBranchAddress("recy", &recy);
  jettree->SetBranchAddress("truex", &truex);
  jettree->SetBranchAddress("truey", &truey);
  jettree->SetBranchAddress("trueq2", &trueq2);
  jettree->SetBranchAddress("truthR1Jets", &truthJets);
  jettree->SetBranchAddress("recoR1Jets", &recoJets);
  jettree->SetBranchAddress("recoR1SDJets", &recoSDJets);
  jettree->SetBranchAddress("matchedR1Jets", &matchedJets);
  jettree->SetBranchAddress("matchedR1SDJets", &matchedSDJets);
  jettree->SetBranchAddress("exchangeBoson", &truthExchangeBoson);
  jettree->SetBranchAddress("matchedParticles", &matchedParticles);
  jettree->SetBranchAddress("smearExchangeBoson", &smearedExchangeBoson);
  jettree->SetBranchAddress("truthR1SDJets", &truthSDJets);
}

float checkdPhi(float dphi)
{
  float newdphi = dphi;
  if(dphi < -1 * PI)
    newdphi += 2. * PI;
  else if(dphi > PI)
    newdphi -= 2. * PI;

  return newdphi;

}
