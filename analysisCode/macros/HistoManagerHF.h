/*
 * A separate header file to keep all of the histogram clutter away from
 * the actual source code
 */

TH1 *truthJetPt, *recoJetPt;
TH2 *matchedJetPt, *matchedJetEta;
TH1 *truthZ, *truthJt, *truthR;
TH1 *recoZ, *recoJt, *recoR;
TH2 *truerecx, *truerecy, *truerecq2;
TH2 *truthxQ2, *recoxQ2, *truthQ2Pt, *recoQ2Pt, *truthnConstPt, *reconConstPt;
TH2 *truthJetPtEta, *truthJetPtPhi;
TH2 *recoJetPtEta, *recoJetPtPhi;
TH2 *truthRecX, *truthRecY, *truthRecQ2;
TH1 *truthBosonCosTheta;

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
  truthJetPt = new TH1F("truthJetPt","; p_{T} [GeV]; Counts", 40, 0., 20.);
  truthJetPt->Sumw2();

  recoJetPt = new TH1F("recoJetPt","; p_{T} [GeV]; Counts", 40, 0., 20.);
  recoJetPt->Sumw2();

  matchedJetPt = new TH2F("matchedJetPt","; true p_{T} [GeV]; reco p_{T} [GeV]", 40, 0., 20., 40, 0., 20.);

  matchedJetEta = new TH2F("matchedJetEta","; true #eta; reco #eta", 40, -4., 4., 40, -4., 4.);

  truthZ = new TH1F("truthZ","; z; Counts", 40, 0., 1.);
  truthZ->Sumw2();

  truthJt = new TH1F("truthJt","; j_{T}; Counts", 40, 0., 1.);
  truthJt->Sumw2();

  truthR = new TH1F("truthR","; r; Counts", 40, 0., 1.);
  truthR->Sumw2();

  recoZ = new TH1F("recoZ","; z; Counts", 40, 0., 1.);
  recoZ->Sumw2();

  recoJt = new TH1F("recoJt","; j_{T}; Counts", 40, 0., 1.);
  recoJt->Sumw2();

  recoR = new TH1F("recoR","; r; Counts", 40, 0., 1.);
  recoR->Sumw2();

  truthJetPtEta = new TH2F("truthJetPtEta",";p_{T}^{true} [GeV]; #eta",20,4,24,50,-3,3);

  truthJetPtPhi = new TH2F("truthJetPtPhi",";p_{T}^{true} [GeV]; #phi [rad]",20,4,24,50,-3.14159,3.14159);

  recoJetPtEta = new TH2F("recoJetPtEta",";p_{T}^{rec} [GeV]; #eta",20,4,24,50,-3,3);

  recoJetPtPhi = new TH2F("recoJetPtPhi",";p_{T}^{rec} [GeV]; #phi [rad]",20,4,24,50,-3.14159,3.14159);
  truthnConstPt = new TH2F("truthNConstPt","; N constituents; p_{T} [GeV]", 40, 0., 40., 40, 0., 20.);

  reconConstPt = new TH2F("recoNConstPt","; N constituents; p_{T} [GeV]", 40, 0., 40., 40, 0., 20.);

  truthQ2Pt = new TH2F("truthQ2Pt", "; Q^[2} [GeV^{2}]; p_{T} [GeV]", 40, 0., 1000., 40, 0., 20.);

  recoQ2Pt = new TH2F("recoQ2Pt", "; Q^{2} [GeV^{2}]; p_{T} [GeV]", 40, 0., 1000., 40, 0., 20.);

  truthxQ2 = new TH2F("truthxQ2", "; x; Q^{2} [GeV^{2}]", 40, 0., 1., 40, 0., 1000.);

  recoxQ2 = new TH2F("recoxQ2", "; x; Q^{2} [GeV^{2}]", 40, 0., 1., 40, 0., 1000.);

  truthBosonCosTheta = new TH1F("bosonCosTheta","; cos#theta; Counts", 40, -2., 2.);
  truthBosonCosTheta->Sumw2();

  truthRecX = new TH2F("truthRecX",";x_{true}; x_{rec}",nxbins, xbins, nxbins, xbins);

  truthRecY = new TH2F("truthRecY",";y_{true}; y_{rec}",100,0,1,100,0,1);

  truthRecQ2 = new TH2F("truthRecQ2",";Q^{2}_{true} [GeV^{2}]; Q^{2}_{reco} [GeV^{2}]",nq2bins,qbins,nq2bins,qbins);

}
