/*
 * A separate header file to keep all of the histogram clutter away from
 * the actual source code
 */

TH2 *truerecoz, *truerecojt, *truerecor;
TH2 *recojetpteta, *recojetptphi;
TH2 *truerecx, *truerecy, *truerecq2;
TH2 *trueQ2x, *trueQ2pT;
TH2 *truejetpteta, *truejetptphi, *trueD0peta;
TH2 *recojetptz, *recojetptjt, *recojetptr;
TH2 *truejetptz, *truejetptjt, *truejetptr;
TH2 *recojetptetatruejetpt;
TH2 *recotruejetpt, *recotruejeteta, *recotruejetphi, 
  *recotruejetp, *recotruejete;
TH1 *matchedJetDr, *matchedJetdEta, *matchedJetdPhi;
TH1 *truthRecoConstdPhi, *truthRecoConstdEta, *truthRecoConstdRap;
TH2 *truthSDjetzg, *truthSDjetrg;
TH2 *recoSDjetzg, *recoSDjetrg;
TH2 *truthrecozg, *truthrecorg;
TH1 *sdenergygroomed;
TH1 *truenjetevent, *reconjetevent;
TH1 *trueconstmass, *truepairmass, *truethreebodymass;
TH1 *recoconstmass;
TH1 *truejetmass, *recojetmass;
TH2 *truejetnconstpt, *recojetnconstpt;
TH2 *trued0decaypartspeta;

void write(std::string infileName)
{
  std::string outfileName = "../histos/";

 if( infileName.find("100Minclusive") != string::npos )
   {
     outfileName += "100Minclusive/";
   }

  if( infileName.find("breit") != string::npos )
    outfileName += "breitJetHistos";
  else 
    outfileName += "labJetHistos";

  if( infileName.find("275") != string::npos )
    outfileName += "_pE275_eE18";
  else if( infileName.find("100") != string::npos )
    outfileName += "_pE100_eE10";

 if( infileName.find("_d0") != string::npos )
    outfileName += "_d0";

 if( infileName.find("kpi") != string::npos )
    outfileName += "kpi";

 if( infileName.find("100Minclusive") != string::npos )
   {
     for (int i=0; i<100; ++i)
       {
	 string index = "";
	 string tag = "";
	 index = std::to_string(i);
	 tag = "_" + index + ".";
	 if( infileName.find( tag )!= string::npos )
	   outfileName += "_" + index;
       }
   }


  outfileName += ".root";
  outfile = new TFile(outfileName.c_str(),"RECREATE");  

  sdenergygroomed->Write();
  truthSDjetzg->Write();
  truthSDjetrg->Write();
  recoSDjetzg->Write();
  recoSDjetrg->Write();
  truthrecozg->Write();
  truthrecorg->Write();

  truerecx->Write();
  truerecy->Write();
  truerecq2->Write();
  trueQ2x->Write();
  
  trueQ2pT->Write();
  truejetpteta->Write();
  truejetptphi->Write();
  recojetptr->Write();
  recojetptjt->Write();
  recojetptz->Write();
  truejetptr->Write();
  truejetptjt->Write();
  truejetptz->Write();
  recojetpteta->Write();
  recojetptphi->Write();
  recotruejetpt->Write();
  recotruejetphi->Write();
  recotruejeteta->Write();
  recotruejetp->Write();
  recotruejete->Write();
  
  truerecoz->Write();
  truerecojt->Write();
  truerecor->Write();
  matchedJetDr->Write();
  matchedJetdPhi->Write();
  matchedJetdEta->Write();
  truthRecoConstdPhi->Write();
  truthRecoConstdEta->Write();
  truthRecoConstdRap->Write();
  recojetptetatruejetpt->Write();
  truenjetevent->Write();
  reconjetevent->Write();
  truejetmass->Write();
  recojetmass->Write();
  trueconstmass->Write();
  recoconstmass->Write();
  //truepairmass->Write();
  //truethreebodymass->Write();
  truejetnconstpt->Write();
  recojetnconstpt->Write();
  trueD0peta->Write();
  trued0decaypartspeta->Write();


  outfile->Write();
  outfile->Close();

}
void instantiateHistos()
{
  sdenergygroomed = new TH1F("sdenergygroomed",";E_{SD}/E_{AKT}",101,0,1.01);
  sdenergygroomed->Sumw2();

  truthrecozg = new TH2F("truthrecozg",";z_{g}^{true}; z_{g}^{reco}",
			 nzgbins, zgbins, nzgbins, zgbins);
  truthrecorg = new TH2F("truthrecorg",";R_{g}^{true}; R_{g}^{reco}",
			 nrbins, rbins, nrbins, rbins);
  truthSDjetzg = new TH2F("truthSDjetzg",
			  ";z_{g}^{true};p_{T}^{Soft Drop, true} [GeV]",
			  nzgbins,zgbins, nptbins, ptbins);
  truthSDjetrg = new TH2F("truthSDjetrg",
			  ";R_{g}^{true};p_{T}^{Soft Drop, true} [GeV]",
			  nrbins,rbins,nptbins,ptbins);
  recoSDjetzg = new TH2F("recoSDjetzg",
			  ";z_{g}^{reco};p_{T}^{Soft Drop, reco} [GeV]",
			  nzgbins,zgbins, nptbins, ptbins);
  recoSDjetrg = new TH2F("recoSDjetrg",
			  ";R_{g}^{reco};p_{T}^{Soft Drop, reco} [GeV]",
			  nrbins,rbins,nptbins,ptbins);



  truerecx = new TH2F("truerecx",";x_{true}; x_{rec}",nxbins, xbins, 
		      nxbins, xbins);
  truerecy = new TH2F("truerecy",";y_{true}; y_{rec}",100,0,1,100,0,1);
  truerecq2 = new TH2F("truerecq2",
		       ";Q^{2}_{true} [GeV^{2}]; Q^{2}_{reco} [GeV^{2}]",
		       nq2bins,qbins,nq2bins,qbins);
  trueQ2x = new TH2F("trueq2x",";x_{true};Q^{2}_{true} [GeV^{2}]",
		     nxbins,xbins,nq2bins,qbins);
  
  trueQ2pT = new TH2F("trueq2pt",";Q_{true}^{2} [GeV^{2}]; p_{T}^{true} [GeV]",
		      nq2bins,qbins,nptbins,ptbins);
  truejetpteta = new TH2F("truejetpteta",";p_{T}^{true} [GeV]; #eta",
			  30,4,34,70,-3.5,3.5);
  trueD0peta = new TH2F("trueD0peta",";p^{true}_{D^{0}} [GeV]; #eta",
			  100,0,100,70,-3.5,3.5);
  truejetptphi = new TH2F("truejetptphi",";p_{T}^{true} [GeV]; #phi [rad]",
			  30,4,34,50,-3.14159,3.14159);
  recojetptz = new TH2F("recojetptz",";z_{reco};p_{T}^{jet,reco} [GeV]",
			nzbins,zbins,nptbins,ptbins);
  recojetptjt = new TH2F("recojetptjt",
			 ";j_{T}^{reco} [GeV]; p_{T}^{jet,reco} [GeV]",
			 njtbins,jtbins,nptbins,ptbins);
  recojetptr = new TH2F("recojetptr",";r_{reco}; p_{T}^{jet,reco} [GeV]",
			nrbins,rbins,nptbins,ptbins);
  truejetptz = new TH2F("truejetptz",";z_{true};p_{T}^{jet,true} [GeV]",
			nzbins,zbins,nptbins,ptbins);
  truejetptjt = new TH2F("truejetptjt",
			 ";j_{T}^{true} [GeV]; p_{T}^{jet,true} [GeV]",
			 njtbins,jtbins,nptbins,ptbins);
  truejetptr = new TH2F("truejetptr",";r_{true}; p_{T}^{jet,true} [GeV]",
			nrbins,rbins,nptbins,ptbins);
  recojetpteta = new TH2F("recojetpteta", ";p_{T}^{jet,reco} [GeV]; #eta",
			  30,4,34,50,-3,3);
  recojetptphi = new TH2F("recojetptphi",";p_{T}^{jet,reco} [GeV]; #phi",
			  30,4,34,50,-3.14159,3.14159);
  recotruejetpt = new TH2F("recotruejetpt",";p_{T}^{jet,true} [GeV]; p_{T}^{jet,reco} [GeV]",
			   nptbins, ptbins, nptbins, ptbins);
  recotruejetphi = new TH2F("recotruejetphi",";#phi_{true} [rad];#phi_{reco} [rad]",
			    200, -3.14159,3.14159, 200,-3.14159, 3.14159);
  recotruejeteta = new TH2F("recotruejeteta",";#eta_{true}; #eta_{reco}",
			    200,-3,3, 200,-3,3);
  recotruejetp = new TH2F("recotruejetp",";p^{jet,true} [GeV];p^{jet,reco} [GeV]", npbins, pbins, npbins, pbins);
  recotruejete = new TH2F("recotruejete",";E^{jet,true} [GeV]; E^{jet,reco} [GeV]", npbins, pbins, npbins, pbins);

  truerecoz = new TH2F("truerecoz",";z_{true};z_{reco}",nzbins,zbins,nzbins,zbins);
  truerecojt = new TH2F("truerecojt",";j_{T}^{true} [GeV];j_{T}^{reco} [GeV]",njtbins,jtbins,njtbins,jtbins);
  truerecor = new TH2F("truerecor",";r_{true}; r_{reco}",nrbins,rbins,nrbins,rbins);

  matchedJetDr = new TH1F("matchedJetDr",";#DeltaR(true,reco)",100,0,1.1);
  matchedJetDr->Sumw2();

  matchedJetdEta = new TH1F("matchedJetdEta",";#Delta#eta(true,reco)",100,-0.5,0.5);
  matchedJetdEta->Sumw2();

  matchedJetdPhi = new TH1F("matchedJetdPhi",";#Delta#phi(true,reco) [rad]",100,-0.5,0.5);
  matchedJetdPhi->Sumw2();

  truthRecoConstdPhi = new TH1F("truthRecoConstdPhi",";#Delta#phi(true_{const},reco_{const}) [rad]", 100,-0.5,0.5);

  truthRecoConstdEta = new TH1F("truthRecoConstdEta",";#Delta#eta(true_{const}, reco_{const})",100,-0.5,0.5);

  truthRecoConstdRap = new TH1F("truthRecoConstdRap",";#Deltay(true_{const}, reco_{const})",100,-0.5,0.5);

  recojetptetatruejetpt = new TH2F("recojetptetatruejetpt",";p_{T}^{jet,true} [GeV]; #eta^{jet,true}",30,4,34,50,-3,3);

  truenjetevent = new TH1I("truenjetevent", "; N^{jet, true}; Counts", 10, 0, 10);

  reconjetevent = new TH1I("reconjetevent", "; N^{jet, true}; Counts", 10, 0, 10);

  trueconstmass = new TH1F("trueconstmass", ";M_{con} [GeV]; ", 120, 1., 3.);

  recoconstmass = new TH1F("recoconstmass", ";M_{con} [GeV]; ", 120, 1., 3.);

  truejetmass = new TH1F("truejetmass", ";M_{jet}^{true} [GeV]; ", 100, 0., 20.);

  recojetmass = new TH1F("recojetmass", ";M_{jet}^{reco} [GeV]; ", 100, 0., 20.);

  truepairmass = new TH1F("truepairmass", ";M_{pair} [GeV]; ", 40, 1, 3.);

  truethreebodymass = new TH1F("truethreebodymass", ";M_{pair} [GeV]; ", 600, 1.8, 2.4);

  recojetnconstpt = new TH2F ("recojetnconstpt", ";N_{const}^{reco}; p_{T}^{reco} [GeV] ", 20, 0, 20, 26, 4., 30.);
  
  truejetnconstpt = new TH2F ("truejetnconstpt", ";N_{const}^{true}; p_{T}^{true} [GeV] ", 20, 0, 20, 26, 4., 30.);

  trued0decaypartspeta = new TH2F("trued0decaypartspeta",";p^{true} [GeV]; #eta",
			  60,0,60,70,-3.5,3.5);
}
