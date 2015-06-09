#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/baselineDef.h"
#include "TStopwatch.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Closure.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"
#include "TauResponse.h"
#include "utils.h"
#include "Efficiency.h"
#include "TRandom3.h"

using namespace std;

static const int nSR = 1;

void passBaselineFunc1(NTupleReader &tr)
{
  bool passBaseline = true;
  //Form TLorentzVector of MET
  TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr.getVar<double>("met"), 0, tr.getVar<double>("metphi"), 0);
  //Calculate number of leptons
  int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsMiniIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsMiniIsoArr);
  int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesMiniIso"), tr.getVec<double>("elesMtw"), tr.getVec<unsigned int>("elesisEB"), AnaConsts::elesArr);
  //  int nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<double>("loose_isoTrks_iso"), tr.getVec<double>("loose_isoTrks_mtw"), AnaConsts::isoTrksArr);

  //Calculate number of jets and b-tagged jets
  int cntCSVS = AnaFunctions::countCSVS(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
  int cntNJetsPt50Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt50Eta24Arr);
  int cntNJetsPt30Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Eta24Arr);
  int cntNJetsPt30 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Arr);

  //Calculate deltaPhi
  std::vector<double> * dPhiVec = new std::vector<double>();
  (*dPhiVec) = AnaFunctions::calcDPhi(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVar<double>("metphi"), 3, AnaConsts::dphiArr);

  std::vector<TLorentzVector> *jetsLVec_forTagger = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger = new std::vector<double>();
  AnaFunctions::prepareJetsForTagger(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), (*jetsLVec_forTagger), (*recoJetsBtag_forTagger));


  //Pass lepton veto?
  bool passLeptVeto = true;
  if( nMuons != AnaConsts::nMuonsSel ){ passBaseline = false; passLeptVeto = false; }
  if( nElectrons != AnaConsts::nElectronsSel ){ passBaseline = false; passLeptVeto = false; }

  //Pass number of jets?
  bool passnJets = true;
  if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){ passBaseline = false; passnJets = false;}
  if( cntNJetsPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 ){ passBaseline = false; passnJets = false;}

  //Pass deltaPhi?
  bool passdPhis = true;
  if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ){ passBaseline = false; passdPhis = false; }

  //Pass number of b-tagged jets?
  bool passBJets = true;
  if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ){ passBaseline = false; passBJets = false; }
  
  //Pass the baseline MET requirement?
  bool passMET = true;
  if( tr.getVar<double>("met") < AnaConsts::defaultMETcut ){ passBaseline = false; passMET = false; }

    //Calculate top tagger related variables.
  //Note that to save speed, only do the calculation after previous base line requirements.
  int nTopCandSortedCnt = -1;
  double MT2 = -1;
  double mTcomb = -1;
  if( passnJets && cntNJetsPt30 >= AnaConsts::nJetsSel ){
      type3Ptr->processEvent((*jetsLVec_forTagger), (*recoJetsBtag_forTagger), metLVec);
      nTopCandSortedCnt = type3Ptr->nTopCandSortedCnt;
      MT2 = type3Ptr->best_had_brJet_MT2;
      //      mTcomb = type3Ptr->mTbJet + 0.5*type3Ptr->mTbestTopJet;
      mTcomb = type3Ptr->best_had_brJet_mTcomb;
  }

  //Pass top tagger requirement?
  bool passTagger = type3Ptr->passNewTaggerReq();
  if( !passTagger ) passBaseline = false;
  bool passNewCuts = type3Ptr->passNewCuts();

  //register new var

  tr.registerDerivedVar("nMuons_CUT2", nMuons);
  tr.registerDerivedVar("nElectrons_CUT2", nElectrons);
  tr.registerDerivedVar("cntNJetsPt30Eta24", cntNJetsPt30Eta24);

  tr.registerDerivedVar("passBaseline", passBaseline);
  tr.registerDerivedVar("passLeptVeto", passLeptVeto);
  tr.registerDerivedVar("passMET", passMET);
  tr.registerDerivedVar("passnJets", passnJets);
  tr.registerDerivedVar("passdPhis", passdPhis);
  tr.registerDerivedVar("passBJets", passBJets);
  tr.registerDerivedVar("passTagger", passTagger);
  tr.registerDerivedVar("cntCSVS", cntCSVS);
  tr.registerDerivedVar("nTopCandSortedCnt", nTopCandSortedCnt);
  tr.registerDerivedVar("MT2_new", MT2);
  tr.registerDerivedVar("mTcomb_new", mTcomb);
  tr.registerDerivedVar("passNewCuts", passNewCuts);
}

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 3)
    {
      std::cerr <<"Please give 3 arguments " << "inputList " << " " <<" "<<"input template"<<" "<< "outputFileName" << std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Closure List1_ttbar.txt HadTau_TauResponseTemplates.root HadTau_Closure.root" << std::endl;
      return -1;
    }
  const char *inputFileList = argv[1];
  const char *respTempl = argv[2];  
  const char *outFileName = argv[3];
  TChain *fChain = new TChain("stopTreeMaker/AUX");
  if(!FillChain(fChain, inputFileList))
    {
      std::cerr << "Cannot get the tree " << std::endl;
    }

  NTupleReader tr(fChain);
  AnaFunctions::prepareTopTagger();
  tr.registerFunction(&passBaselineFunc1);
// Add cleanJets function
  stopFunctions::cjh.setMuonIso("mini");
  stopFunctions::cjh.setElecIso("mini");
  tr.registerFunction(&stopFunctions::cleanJets);

  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(outFileName);

  TauResponse tauResp(respTempl);

  TRandom3 * rndm = new TRandom3(12345);
  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;

  int entries = tr.getNEntries();
  std::cout<<"\nentries : "<<entries<<std::endl;

  std::vector<double> trueVec(nSR), predVec(nSR), pred_from_taumuVec(nSR);
  int cnt_nomtw=0,cnt_mtw=0;

  // Loop over the events (tree entries)
  int k = 0;
  while(tr.getNextEvent()){
    k++;

// Add print out of the progress of looping
    if( tr.getEvtNum()-1 == 0 || tr.getEvtNum() == entries || (tr.getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr.getEvtNum()-1<<"th event ..."<<std::endl;

    const vector<TLorentzVector> &muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
    const vector<double> &muonsRelIso = tr.getVec<double>("muonsRelIso");
    const vector<double> &muonsMiniIso = tr.getVec<double>("muonsMiniIso");
    const vector<double> &muonsMtw = tr.getVec<double>("muonsMtw");
    const vector<TLorentzVector> &genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
    const vector<int> &genDecayIdxVec = tr.getVec<int>("genDecayIdxVec");
    const vector<int> &genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
    const vector<int> &W_emuVec = tr.getVec<int>("W_emuVec");
    const vector<int> &W_tau_emuVec = tr.getVec<int>("W_tau_emuVec");
    const vector<TLorentzVector> &jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
    const vector<double> &recoJetsBtag_0 = tr.getVec<double>("recoJetsBtag_0");
    int nElectrons = tr.getVar<int>("nElectrons_CUT2");
    int nMuons = tr.getVar<int>("nMuons_CUT2");
    double met=tr.getVar<double>("met");
    double metphi=tr.getVar<double>("metphi");
    double ht=tr.getVar<double>("ht");
    vector<int> W_tau_prongsVec = tr.getVec<int>("W_tau_prongsVec");

    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);

    //Expectation part -- do it before prediction & do NOT skipping events
    bool passBaseline_tru = tr.getVar<bool>("passBaseline");

    const int nJets_tru = tr.getVar<int>("cntNJetsPt30Eta24");
    const int nbJets_tru = tr.getVar<int>("cntCSVS");
    const int nTops_tru = tr.getVar<int>("nTopCandSortedCnt");
    const double MT2_tru = tr.getVar<double>("MT2_new");
    const double mTcomb_tru = tr.getVar<double>("mTcomb_new");

    //Select only events where the W decayed into a hadronically decaying tau
    //Note that for ttbar this includes (for two W's)
    // ] W->tau->had, W->tau->had
    // ] W->tau->had, W->qq
    // ] W->tau->had, W->e/mu (e or mu is lost)
    if(W_tau_prongsVec.size() !=0 && passBaseline_tru){
// Currently only one signal region which is the baseline
// In the future, this can extend to more regions

       int iSR = 0;
       trueVec[iSR] ++;

       myBaseHistgram.hTrueHt->Fill(ht);
       myBaseHistgram.hTruemet->Fill(met);
       myBaseHistgram.hTrueNJets->Fill(nJets_tru);
       myBaseHistgram.hTrueNbJets->Fill(nbJets_tru);
       myBaseHistgram.hTrueNTops->Fill(nTops_tru);
       myBaseHistgram.hTrueMT2->Fill(MT2_tru);
       myBaseHistgram.hTruemTcomb->Fill(mTcomb_tru);
    }

    //Prediction part
    //Control sample
    // The kinematic properties of the well-reconstructed, isolated muon                                                                    
    vector<TLorentzVector> isomuonsLVec;
    vector<int> isomuonsIdxVec;
    for(unsigned int im=0; im<muonsLVec.size(); im++){
      if( AnaFunctions::passMuon(muonsLVec.at(im), muonsMiniIso.at(im), muonsMtw.at(im), AnaConsts::muonsMiniIsoArr) ){ isomuonsLVec.push_back(muonsLVec.at(im)); isomuonsIdxVec.push_back(im); }
    }

// Require one and only one muon
// Veto events with additional electrons (same veto criteria as baseline for electrons)
    if( nMuons == 1 && nElectrons == AnaConsts::nElectronsSel ) {
      if( nMuons != isomuonsLVec.size() ){ std::cout<<"ERROR ... mis-matching between veto muon and selected muon! Skipping..."<<std::endl; continue; }
      const TLorentzVector muLVec = isomuonsLVec.at(0);
      // Use only events where the muon is inside acceptance                                                                             
      if( muLVec.Pt() < TauResponse::ptMin() ) continue;
      if( fabs(muLVec.Eta()) > TauResponse::etaMax() ) continue;
      //mtW correction
      const double mtw = calcMT(muLVec, metLVec);
      bool pass_mtw = false;
      if(mtw<100)pass_mtw = true;

      // Find events that contain W->tau->mu
      // Note that any code using gen info should be checked if they work for data or not!
      bool istaumu = false, istaumu_genRecoMatch = false;
      for(unsigned int ig=0; ig<W_tau_emuVec.size(); ig++){
         int genIdx = W_tau_emuVec.at(ig);
         if( std::abs(genDecayPdgIdVec.at(genIdx)) == 13 ){
            istaumu = true;
            const TLorentzVector & genLVec = genDecayLVec.at(genIdx);
            if( muLVec.DeltaR(genLVec) < 0.2 ) istaumu_genRecoMatch = true;
         }
      }
//      if( W_emuVec.empty() && istaumu && !istaumu_genRecoMatch ) std::cout<<"WARNING ... reco muon does NOT match to the tau->mu?!"<<std::endl;
      // "Cross cleaning": find the jet that corresponds to the muon
      const std::vector<TLorentzVector>& cleanJetVec      = tr.getVec<TLorentzVector>("cleanJetVec");
      const std::vector<double>& cleanJetBtag             = tr.getVec<double>("cleanJetBTag");
      // Get the cleaned jet indice (pointing to the jetsLVec) for the corresponding muons
      const std::vector<int>& rejectJetIdx_formuVec = tr.getVec<int>("rejectJetIdx_formuVec");
      const double & cleanHt = tr.getVar<double>("cleanHt");
      const double & cleanMHt = tr.getVar<double>("cleanMHt");
      const double & cleanMHtPhi = tr.getVar<double>("cleanMHtPhi");

      TLorentzVector selMhtLVec; selMhtLVec.SetPtEtaPhiM(cleanMHt, 0, cleanMHtPhi, 0);
// Force the mass to be 0 for met and mht
      TLorentzVector selmetLVec; selmetLVec.SetVectM( (metLVec+ muLVec).Vect(), 0 );
      int selNJetPt30Eta24 = AnaFunctions::countJets(cleanJetVec, AnaConsts::pt30Eta24Arr);
      int selNJetPt50Eta24 = AnaFunctions::countJets(cleanJetVec, AnaConsts::pt50Eta24Arr);

      // rejecting events with nJets less than requirements even adding one more tau jet
      if( selNJetPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 - 1 ) continue;
      if( selNJetPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 - 1 ) continue;

      // Get random number from tau-response template
      // The template is chosen according to the muon pt 
      const double scale = tauResp.getRandom(muLVec.Pt());
      // Scale muon pt and energy with tau response --> simulate tau jet pt and energy
      const double simTauJetPt = scale * muLVec.Pt();
      const double simTauJetE = scale * muLVec.E();
      const double simTauJetEta = muLVec.Eta();
      const double simTauJetPhi = muLVec.Phi();

      TLorentzVector tauJetLVec; tauJetLVec.SetPtEtaPhiE(simTauJetPt, simTauJetEta, simTauJetPhi, simTauJetE);
      // Decide the CSV value for the tau jet -> use the CSV of the associated muon-jet as the tau jet CSV
      // Default set to be 0 (low enough to be NOT a b jet)
      double oriJetCSVS = 0;
      if( rejectJetIdx_formuVec.at(isomuonsIdxVec.at(0)) != -1 ) oriJetCSVS = recoJetsBtag_0[rejectJetIdx_formuVec.at(isomuonsIdxVec.at(0))];
      double mistag = Efficiency::mistag(Efficiency::Ptbin1(simTauJetPt));
      double rno = rndm->Rndm();
      if( rno < mistag) oriJetCSVS = 1.0;
      //      cout<<mistag<<" "<<rno<<endl;      

      vector<TLorentzVector> combNJetVec;
      vector<double> combJetsBtag;

      bool includeTauJet = false;
      for(unsigned int ij=0; ij<cleanJetVec.size(); ij++){
         if( tauJetLVec.Pt() > cleanJetVec.at(ij).Pt() && !includeTauJet ){
            combNJetVec.push_back(tauJetLVec); combJetsBtag.push_back(oriJetCSVS);
            includeTauJet = true;
         }
         combNJetVec.push_back(cleanJetVec.at(ij)); combJetsBtag.push_back(cleanJetBtag.at(ij));
      }
      // it's possible that the tau jet is the least energetic jet so that it's not added into the combNJetVec during the loop
      if( !includeTauJet ){ combNJetVec.push_back(tauJetLVec); combJetsBtag.push_back(oriJetCSVS); }

      // Taking into account the simulated tau jet, recompute                                                                                 
      // HT, MHT, and N(jets)                                                                                                                 
      double simHt = cleanHt;
      TLorentzVector simMhtLVec;

      // If simulted tau-jet meets same criteria as as HT jets,                                                                               
      // recompute HT and MH
     
      if( tauJetLVec.Pt() > htJetPtMin() && fabs(tauJetLVec.Eta()) < htJetEtaMax()) {
        simHt += tauJetLVec.Pt();
      }

      if( tauJetLVec.Pt() > mhtJetPtMin() && fabs(tauJetLVec.Eta()) < mhtJetEtaMax()) {
        simMhtLVec.SetVectM( (selMhtLVec-tauJetLVec).Vect(), 0);
      }
      //recompute met                                                                                                                          
      TLorentzVector simmetLVec; simmetLVec.SetVectM( (selmetLVec - tauJetLVec).Vect(), 0);

      const double simMht = simMhtLVec.Pt();

      const double simmet = simmetLVec.Pt();
      const double simmetPhi = simmetLVec.Phi();

      //recompute jetVec

      int combNJetPt30Eta24 = AnaFunctions::countJets(combNJetVec, AnaConsts::pt30Eta24Arr);
      int combNJetPt50Eta24 = AnaFunctions::countJets(combNJetVec, AnaConsts::pt50Eta24Arr);

      //recompute deltaphi

      std::vector<double> * deltaPhiVec = new std::vector<double>();
      (*deltaPhiVec) = AnaFunctions::calcDPhi(combNJetVec, simmetPhi, 3, AnaConsts::dphiArr);
      bool passdeltaPhi = true;
      if( deltaPhiVec->at(0) < AnaConsts::dPhi0_CUT || deltaPhiVec->at(1) < AnaConsts::dPhi1_CUT || deltaPhiVec->at(2) < AnaConsts::dPhi2_CUT){
	  passdeltaPhi = false;
      }

      //recompute bjet

      int cnt1CSVS = AnaFunctions::countCSVS(combNJetVec, combJetsBtag, AnaConsts::cutCSVS, AnaConsts::bTagArr);
      bool passbJets = true;
      if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cnt1CSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cnt1CSVS < AnaConsts::high_nJetsSelBtagged ) ) ){
	 passbJets = false;
      }
	
      //top tagger input
      int comb30_pre = AnaFunctions::countJets(combNJetVec, AnaConsts::pt30Arr);
      std::vector<TLorentzVector> *jetsLVec_forTagger_pre = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger_pre = new std::vector<double>();
      AnaFunctions::prepareJetsForTagger(combNJetVec, combJetsBtag, (*jetsLVec_forTagger_pre), (*recoJetsBtag_forTagger_pre));
        int nTopCandSortedCnt_pre = -1;
	double MT2_pre = -1;
	double mTcomb_pre = -1;
      //Apply baseline cut

      if(combNJetPt30Eta24<AnaConsts::nJetsSelPt30Eta24) continue;
      if(combNJetPt50Eta24<AnaConsts::nJetsSelPt50Eta24) continue;
      if(simmet<AnaConsts::defaultMETcut) continue;
      if(!passdeltaPhi) continue;
      if(!passbJets) continue;

      //Apply Top tagger
      if(comb30_pre >= AnaConsts::nJetsSel ){
	 type3Ptr->processEvent((*jetsLVec_forTagger_pre), (*recoJetsBtag_forTagger_pre), simmetLVec);
         nTopCandSortedCnt_pre = type3Ptr->nTopCandSortedCnt;
         MT2_pre = type3Ptr->best_had_brJet_MT2;
         mTcomb_pre = type3Ptr->best_had_brJet_mTcomb;
      }

      bool passTopTagger = type3Ptr->passNewTaggerReq();
      //New tagger cuts
      bool passNewCuts = type3Ptr->passNewCuts();
      if(!passTopTagger) continue;
      
      //mtw count
      cnt_nomtw++;
      if(pass_mtw)cnt_mtw++;

      //Activity variable calculation:
      double muact = AnaFunctions::getMuonActivity(muLVec, jetsLVec, tr.getVec<double>("recoJetschargedHadronEnergyFraction"), tr.getVec<double>("recoJetschargedEmEnergyFraction"),AnaConsts::muonsAct);

      //correction factor:
      const double corrBRWToTauHad = 0.65;  // Correction for the BR of hadronic tau decays                          
      const double corrBRTauToMu = Efficiency::taumucor(Efficiency::Ptbin1(muLVec.Pt()));//correction from tauonic mu contamination
      const double corrMuAcc = 1./Efficiency::acc(Efficiency::Njetbin(combNJetPt30Eta24)); // Correction for muon acceptance
      const double corrMuRecoEff = 1./Efficiency::reco(Efficiency::Ptbin(muLVec.Pt()), Efficiency::Actbin(muact)); // Correction for muon reconstruction efficiency             
      const double corrMuIsoEff = 1./Efficiency::iso(Efficiency::Ptbin(muLVec.Pt()), Efficiency::Actbin(muact)); // Correction for muon isolation efficiency 

      const double corrmtWEff = 1./0.8897;
      //The overall correction factor                                                                                                          
//      const double corr = corrBRTauToMu * corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff;
// For MC, no need of applying the corrBRTauToMu as you know if this event is from W->tau->mu or not
// However, we need know the fraction so that we can apply it to data (can be got from the printout)
  const double corr = corrBRWToTauHad * corrMuAcc * corrMuRecoEff * corrMuIsoEff * corrmtWEff;

// iSR: this should be determined by search region requirement
// Now we have only one bin which is the baseline
      int iSR = 0;

      predVec[iSR] += corr;
      if( istaumu_genRecoMatch ) pred_from_taumuVec[iSR] += corr;

      // Fill the prediction
      if( !istaumu_genRecoMatch && pass_mtw){
         myBaseHistgram.hPredHt->Fill(simHt,corr);
         myBaseHistgram.hPredmet->Fill(simmet,corr);
         myBaseHistgram.hPredNJets->Fill(combNJetPt30Eta24,corr);
         myBaseHistgram.hPredNbJets->Fill(cnt1CSVS,corr);
         myBaseHistgram.hPredNTops->Fill(nTopCandSortedCnt_pre,corr);
         myBaseHistgram.hPredMT2->Fill(MT2_pre,corr);
         myBaseHistgram.hPredmTcomb->Fill(mTcomb_pre,corr);
      }
    }//control sample loop
  }
  // --- Save the Histograms to File -----------------------------------
  std::cout<<"lastbin:"<<myBaseHistgram.hPredmet->GetBinContent(myBaseHistgram.hPredmet->GetXaxis()->GetNbins())<<std::endl;
  std::cout<<"overflow:"<<myBaseHistgram.hPredmet->GetBinContent(myBaseHistgram.hPredmet->GetXaxis()->GetNbins() + 1)<<std::endl;

  drawOverFlowBin(myBaseHistgram.hPredmet);
  drawOverFlowBin(myBaseHistgram.hTruemet);
  std::cout<<"newlastbin:"<<myBaseHistgram.hPredmet->GetBinContent(myBaseHistgram.hPredmet->GetXaxis()->GetNbins())<<std::endl;
  drawOverFlowBin(myBaseHistgram.hPredMT2);
  drawOverFlowBin(myBaseHistgram.hTrueMT2);
  drawOverFlowBin(myBaseHistgram.hPredNbJets);
  drawOverFlowBin(myBaseHistgram.hTrueNbJets);
  (myBaseHistgram.oFile)->Write();

// This print out can be used to extract the corrBRTauToMu ratio
  std::cout<<"\nPrediction in numbers ..."<<std::endl;
  for(unsigned int ib=0; ib<nSR; ib++){
     std::cout<<"ib : "<<ib<<"  true : "<<trueVec[ib]<<"  pred : "<<predVec[ib]<<"  from_taumu : "<<pred_from_taumuVec[ib]<<std::endl;
  }
  std::cout<<std::endl;

  std::cout<<cnt_nomtw<<"  "<<cnt_mtw<<std::endl;
  return 0;
}
