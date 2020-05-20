#include "analyzeJetsHF.h"
#include "HistoManagerHF.h"

void analyzeJetsHF() 
{
  SetsPhenixStyle(); 

  infileName  = "../dataFiles/labJets_pE275_eE18_oc.root";
  outfileName = "../histos/labJetHistos_pE275_eE18_oc.root";

  gROOT->ProcessLine(".L ../src/fastJetLinker.C+");

  setupTree(infileName);

  instantiateHistos();

  loop();

  write(outfileName);

}

void loop()
{
  for (int nEntry = 0; nEntry < jettree->GetEntries(); ++nEntry)
    {
      jettree->GetEntry(nEntry);

      // Truth Loop //
      for (int nJet = 0; nJet < truthJets->size(); ++nJet)
	{
	  TLorentzVector jetVect(0,0,0,0);
	  jetVect = truthJets->at(nJet).first;
	  truthJetPt->Fill(jetVect.Pt());
	  truthnConstPt->Fill(truthJets->at(nJet).second.size(), jetVect.Pt()); 
	  truthQ2Pt->Fill(trueq2, jetVect.Pt());
	  truthJetPtPhi->Fill(jetVect.Pt(), jetVect.Phi());
	  truthJetPtEta->Fill(jetVect.Pt(), jetVect.Eta());

	  for (int nPart = 0; nPart < truthJets->at(nJet).second.size(); ++nPart)
	    {
	      TLorentzVector partVect(0,0,0,0);
	      partVect = truthJets->at(nJet).second.at(nPart);

	      double r, jt, z;
	      r = 0; z = 0; jt = 0;
	      z = ( jetVect.P() - partVect.P() ) / ( jetVect.P() * jetVect.P() );
	      jt = partVect.Vect().Cross(jetVect.Vect()).Mag()/jetVect.P();
	      r = sqrt ( ( partVect.Phi() - jetVect.Phi() ) *  ( partVect.Phi() - jetVect.Phi() ) +  
                         ( partVect.Rapidity() - jetVect.Rapidity() ) *  ( partVect.Rapidity() - jetVect.Rapidity() ) );


	      truthZ->Fill(z);
	      truthJt->Fill(jt);
	      truthR->Fill(r);
										      
	    }

	}
      
      // Smeared Loop //
      for (int nJet = 0; nJet < recoJets->size(); ++nJet)
	{
	  TLorentzVector jetVect(0,0,0,0);
	  jetVect = recoJets->at(nJet).first;
	  recoJetPt->Fill(jetVect.Pt());
	  reconConstPt->Fill(recoJets->at(nJet).second.size(), jetVect.Pt()); 
	  recoQ2Pt->Fill(recq2, jetVect.Pt());
	  recoJetPtPhi->Fill(jetVect.Pt(), jetVect.Phi());
	  recoJetPtEta->Fill(jetVect.Pt(), jetVect.Eta());
	

	  for (int nPart = 0; nPart < recoJets->at(nJet).second.size(); ++nPart)
	    {
	      TLorentzVector partVect(0,0,0,0);
	      partVect = recoJets->at(nJet).second.at(nPart);

	      double r, jt, z;
	      r = 0; z = 0; jt = 0;
	      z = ( jetVect.P() - partVect.P() ) / ( jetVect.P() * jetVect.P() );
	      jt = partVect.Vect().Cross(jetVect.Vect()).Mag()/jetVect.P();
	      r = sqrt ( ( partVect.Phi() - jetVect.Phi() ) *  ( partVect.Phi() - jetVect.Phi() ) +  
                         ( partVect.Rapidity() - jetVect.Rapidity() ) *  ( partVect.Rapidity() - jetVect.Rapidity() ) );

	      recoZ->Fill(z);
	      recoJt->Fill(jt);
	      recoR->Fill(r);
					
	      
	    }

	}

      // SD Loop //
      for (int nJet = 0; nJet < recoSDJets->size(); ++nJet)
	{
	  TLorentzVector jetVect(0,0,0,0);
	  jetVect = recoSDJets->at(nJet).first;
	  
	  PseudoJet softDropJet(jetVect.Px(), jetVect.Py(), jetVect.Pz(), jetVect.E() );

	  // Unfortunately we do not have access to jet substructure without reclustering.. 
	  // Should we recluser or pass in variables we need (i.e. Rg and zg)?
	  //float  subJetDR = softDropJet.structure_of<contrib::SoftDrop>().delta_R();
	  //float        zg = softDropJet.structure_of<contrib::SoftDrop>().symmetry();
	  //float mass_drop = softDropJet.structure_of<contrib::SoftDrop>().mu();




	  for (int nPart = 0; nPart < recoSDJets->at(nJet).second.size(); ++nPart)
	    {
	      TLorentzVector partVect(0,0,0,0);
	      partVect = recoSDJets->at(nJet).second.at(nPart);					
	      
	    }

	}

      // Matched Loop //
      for (int nJetPair = 0; nJetPair < matchedJets->size(); ++nJetPair)
	{
	  TLorentzVector truthJetVect(0,0,0,0);
	  TLorentzVector recoJetVect(0,0,0,0);
	  truthJetVect = matchedJets->at(nJetPair).at(0).first;
	  recoJetVect = matchedJets->at(nJetPair).at(1).first;
	  matchedJetPt->Fill(truthJetVect.Pt(),recoJetVect.Pt());
	  matchedJetEta->Fill(truthJetVect.Eta(),recoJetVect.Eta());

	  /*
	  for (int nPart = 0; nPart < matchedets->at(nJetPair).second.size(); ++nPart)
	    {
	      TLorentzVector partVect(0,0,0,0);
	      partVect = recoJets->at(nJetPair).second.at(nPart);
	      double r, jt, z;
	      r = 0; z = 0; jt = 0;
	      z = ( jetVect.P() - partVect.P() ) / ( jetVect.P() * jetVect.P() );
	      jt = partVect.Vect().Cross(jetVect.Vect()).Mag()/jetVect.P();
	      r = sqrt ( ( partVect.Phi() - jetVect.Phi() ) *  ( partVect.Phi() - jetVect.Phi() ) +  
                         ( partVect.Eta() - jetVect.Eta() ) *  ( partVect.Eta() - jetVect.Eta() ) );
	      recoZ->Fill(z);
	      recoJt->Fill(jt);
	      recoR->Fill(r);
					
	      
	    }
	  */
	}
      
      truthxQ2->Fill(truex, trueq2);
      recoxQ2->Fill(recx, recq2);
      truthRecX->Fill(truex,recx);
      truthRecY->Fill(truey,recy);
      truthRecQ2->Fill(trueq2,recq2);
      // This should always be -1 in the Breit frame
      truthBosonCosTheta->Fill(truthExchangeBoson->CosTheta());

    }


}


void setupTree(string inf)
{

  infile = TFile::Open(inf.c_str());

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


}
