#include "analyzeJetsHF.h"

void analyzeJetsHF() 
{
  
  infileName  = "../dataFiles/labJets_pE275_eE18_oc.root";
  outfileName = "../histos/labJetHistos_pE275_eE18_oc.root";

  gROOT->ProcessLine(".L ../src/fastJetLinker.C+");

  setupTree(infileName);

  instantiateHistos();

  loop();

  write(outfileName);

}


void write(string outf)
{
  outfile = new TFile(outf.c_str(),"RECREATE");
  
  truthJetPt->Write();
  recoJetPt->Write();
  matchedJetPt->Write();
  matchedJetEta->Write();
  truthJetPtEta->Write();
  truthJetPtPhi->Write();
  recoJetPtEta->Write();
  recoJetPtPhi->Write();
  truthnConstPt->Write();
  reconConstPt->Write();
  truthQ2Pt->Write();
  recoQ2Pt->Write();

  truthZ->Write();
  truthJt->Write();
  truthR->Write();
  recoZ->Write();
  recoJt->Write();
  recoR->Write();

  truthxQ2->Write();
  recoxQ2->Write();
  truthRecX->Write();
  truthRecY->Write();
  truthRecQ2->Write();
  
  truthBosonCosTheta->Write();

  outfile->Write();

  outfile->Close();

}


void instantiateHistos()
{
  truthJetPt = new TH1F("truthJetPt","truth jet p_{T}; p_{T} (GeV); Counts", 40, 0., 20.);
  truthJetPt->Sumw2();

  recoJetPt = new TH1F("recoJetPt","reco jet p_{T}; p_{T} (GeV); Counts", 40, 0., 20.);
  recoJetPt->Sumw2();

  matchedJetPt = new TH2F("matchedJetPt","matched jet p_{T}; true p_{T} (GeV); reco p_{T} (GeV)", 40, 0., 20., 40, 0., 20.);

  matchedJetEta = new TH2F("matchedJetEta","matched jet #eta; true #eta; reco #eta", 40, -4., 4., 40, -4., 4.);

  truthZ = new TH1F("truthZ","truth z; z; Counts", 40, 0., 1.);
  truthZ->Sumw2();

  truthJt = new TH1F("truthJt","truth j_{T}; j_{T}; Counts", 40, 0., 1.);
  truthJt->Sumw2();

  truthR = new TH1F("truthR","truth r; r; Counts", 40, 0., 1.);
  truthR->Sumw2();

  recoZ = new TH1F("recoZ","reco z; z; Counts", 40, 0., 1.);
  recoZ->Sumw2();

  recoJt = new TH1F("recoJt","reco j_{T}; j_{T}; Counts", 40, 0., 1.);
  recoJt->Sumw2();

  recoR = new TH1F("recoR","reco r; r; Counts", 40, 0., 1.);
  recoR->Sumw2();

  truthJetPtEta = new TH2F("truthJetPtEta",";p_{T}^{true} [GeV]; #eta",20,4,24,50,-3,3);

  truthJetPtPhi = new TH2F("truthJetPtPhi",";p_{T}^{true} [GeV]; #phi [rad]",20,4,24,50,-3.14159,3.14159);

  recoJetPtEta = new TH2F("recoJetPtEta",";p_{T}^{rec} [GeV]; #eta",20,4,24,50,-3,3);

  recoJetPtPhi = new TH2F("recoJetPtPhi",";p_{T}^{rec} [GeV]; #phi [rad]",20,4,24,50,-3.14159,3.14159);
  truthnConstPt = new TH2F("truthNConstPt","truth N_{con} p_{T}; N constituents; p_{T} (GeV)", 40, 0., 40., 40, 0., 20.);

  reconConstPt = new TH2F("recoNConstPt","reco N_{con} p_{T}; N constituents; p_{T} (GeV)", 40, 0., 40., 40, 0., 20.);

  truthQ2Pt = new TH2F("truthQ2Pt", "truth Q^{2} p_{T}; Q^{2} (GeV^{2}); p_{T} (GeV)", 40, 0., 1000., 40, 0., 20.);

  recoQ2Pt = new TH2F("recoQ2Pt", "reco Q^{2} p_{T}; Q^{2} (GeV^{2}); p_{T} (GeV)", 40, 0., 1000., 40, 0., 20.);

  truthxQ2 = new TH2F("truthxQ2", "truth x Q^{2}; x; Q^{2} (GeV^{2})", 40, 0., 1., 40, 0., 1000.);

  recoxQ2 = new TH2F("recoxQ2", "reco x Q^{2}; x; Q^{2} (GeV^{2})", 40, 0., 1., 40, 0., 1000.);

  truthBosonCosTheta = new TH1F("bosonCosTheta","exchange boson cos#theta; cos#theta; Counts", 40, -2., 2.);
  truthBosonCosTheta->Sumw2();

  truthRecX = new TH2F("truthRecX",";x_{true}; x_{rec}",nxbins, xbins, nxbins, xbins);

  truthRecY = new TH2F("truthRecY",";y_{true}; y_{rec}",100,0,1,100,0,1);

  truthRecQ2 = new TH2F("truthRecQ2",";Q^{2}_{true} [GeV^{2}]; Q^{2}_{reco} [GeV^{2}]",nq2bins,qbins,nq2bins,qbins);

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
