#ifndef ANABASELINE_BASELINEDEF_H
#define ANABASELINE_BASELINEDEF_H

#include "NTupleReader.h"
#include "customize.h"

#include "Math/VectorUtil.h"

#include <sstream>
#include <iostream>
#include <fstream>

class BaselineVessel
{
private:
    const std::string spec;

public:
    BaselineVessel(const std::string specialization = "") : spec(specialization) { }

    void passBaseline(NTupleReader &tr)
    {
        bool debug = false;
        bool doIsoTrksVeto = true;
        bool doMuonVeto = true;
        bool incZEROtop = false;

        bool passBaseline = true;
        bool passBaselineNoTag = true;

        std::string jetVecLabel = "jetsLVec";
        std::string CSVVecLabel = "recoJetsBtag_0";
        std::string METLabel    = "met";
        std::string METPhiLabel = "metphi";
        if( spec.compare("noIsoTrksVeto") == 0)
        {
           doIsoTrksVeto = false;
        }
        if( spec.compare("incZEROtop") == 0)
        {
           incZEROtop = true;
        }
        if( spec.compare("hadtau") == 0)
        {
           doMuonVeto = false;
        }
        if( spec.compare("lostlept") == 0)
        {
           doMuonVeto = false;
        }
        if(spec.compare("Zinv") == 0) 
        {
            jetVecLabel = "cleanJetpt30ArrVec";//"jetsLVec";//"prodJetsNoMu_jetsLVec";
            CSVVecLabel = "cleanJetpt30ArrBTag";//"recoJetsBtag_0";
            METLabel    = "cleanMetPt";
            METPhiLabel = "cleanMetPhi";
            doMuonVeto  = false;
        }

        // Form TLorentzVector of MET
        TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr.getVar<double>(METLabel), 0, tr.getVar<double>(METPhiLabel), 0);

        // Calculate number of leptons
        int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsMiniIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsMiniIsoArr);
        int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesMiniIso"), tr.getVec<double>("elesMtw"), tr.getVec<unsigned int>("elesisEB"), AnaConsts::elesMiniIsoArr);
        int nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<double>("loose_isoTrks_iso"), tr.getVec<double>("loose_isoTrks_mtw"), tr.getVec<int>("loose_isoTrks_pdgId"));

        // Calculate number of jets and b-tagged jets
        int cntCSVS = AnaFunctions::countCSVS(tr.getVec<TLorentzVector>(jetVecLabel), tr.getVec<double>(CSVVecLabel), AnaConsts::cutCSVS, AnaConsts::bTagArr);
        int cntNJetsPt50Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>(jetVecLabel), AnaConsts::pt50Eta24Arr);
        int cntNJetsPt30Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>(jetVecLabel), AnaConsts::pt30Eta24Arr);
        int cntNJetsPt30      = AnaFunctions::countJets(tr.getVec<TLorentzVector>(jetVecLabel), AnaConsts::pt30Arr);

        // Calculate deltaPhi
        std::vector<double> * dPhiVec = new std::vector<double>();
        (*dPhiVec) = AnaFunctions::calcDPhi(tr.getVec<TLorentzVector>(jetVecLabel), metLVec.Phi(), 3, AnaConsts::dphiArr);

        // Prepare jets and b-tag working points for top tagger
        std::vector<TLorentzVector> *jetsLVec_forTagger = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger = new std::vector<double>();
        AnaFunctions::prepareJetsForTagger(tr.getVec<TLorentzVector>(jetVecLabel), tr.getVec<double>(CSVVecLabel), (*jetsLVec_forTagger), (*recoJetsBtag_forTagger));
        if( debug ) std::cout<<"\njetsLVec_forTagger->size : "<<jetsLVec_forTagger->size()<<"  recoJetsBtag_forTagger->size : "<<recoJetsBtag_forTagger->size()<<"  passBaseline : "<<passBaseline<<std::endl;

        // Pass lepton veto?
        bool passLeptVeto = true, passMuonVeto = true, passEleVeto = true, passIsoTrkVeto = true;
        if( doMuonVeto && nMuons != AnaConsts::nMuonsSel ){ passBaseline = false; passBaselineNoTag = false; passLeptVeto = false; passMuonVeto = false; }
        if( nElectrons != AnaConsts::nElectronsSel ){ passBaseline = false; passBaselineNoTag = false; passLeptVeto = false; passEleVeto = false; }
        // Isolated track veto is disabled for now
        if( doIsoTrksVeto && nIsoTrks != AnaConsts::nIsoTrksSel ){ passBaseline = false; passBaselineNoTag = false; passLeptVeto = false; passIsoTrkVeto = false; }
        if( debug ) std::cout<<"nMuons : "<<nMuons<<"  nElectrons : "<<nElectrons<<"  nIsoTrks : "<<nIsoTrks<<"  passBaseline : "<<passBaseline<<std::endl;

        // Pass number of jets?
        bool passnJets = true;
        if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){ passBaseline = false; passBaselineNoTag = false; passnJets = false; }
        if( cntNJetsPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 ){ passBaseline = false; passBaselineNoTag = false; passnJets = false; }
        if( debug ) std::cout<<"cntNJetsPt50Eta24 : "<<cntNJetsPt50Eta24<<"  cntNJetsPt30Eta24 : "<<cntNJetsPt30Eta24<<"  cntNJetsPt30 : "<<cntNJetsPt30<<"  passBaseline : "<<passBaseline<<std::endl;

        // Pass deltaPhi?
        bool passdPhis = true;
        if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ){ passBaseline = false; passBaselineNoTag = false; passdPhis = false; }
        if( debug ) std::cout<<"dPhi0 : "<<dPhiVec->at(0)<<"  dPhi1 : "<<dPhiVec->at(1)<<"  dPhi2 : "<<dPhiVec->at(2)<<"  passBaseline : "<<passBaseline<<std::endl;

        // Pass number of b-tagged jets?
        bool passBJets = true;
        if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ){ passBaseline = false; passBJets = false; }
        if( debug ) std::cout<<"cntCSVS : "<<cntCSVS<<"  passBaseline : "<<passBaseline<<std::endl;

        // Pass the baseline MET requirement?
        bool passMET = true;
        if( metLVec.Pt() < AnaConsts::defaultMETcut ){ passBaseline = false; passBaselineNoTag = false; passMET = false; }
        if( debug ) std::cout<<"met : "<<tr.getVar<double>("met")<<"  defaultMETcut : "<<AnaConsts::defaultMETcut<<"  passBaseline : "<<passBaseline<<std::endl;

        // Calculate top tagger related variables. 
        // Note that to save speed, only do the calculation after previous base line requirements.
        int nTopCandSortedCnt = -1;

        if( passnJets && cntNJetsPt30 >= AnaConsts::nJetsSel ){
            type3Ptr->processEvent((*jetsLVec_forTagger), (*recoJetsBtag_forTagger), metLVec);
            nTopCandSortedCnt = type3Ptr->nTopCandSortedCnt;
        }

        // Pass top tagger requirement?
        bool passTagger = type3Ptr->passNewTaggerReq() && (incZEROtop || nTopCandSortedCnt >= AnaConsts::low_nTopCandSortedSel);

        if( !passTagger ) passBaseline = false;

//        bool passNewCuts = type3Ptr->passNewCuts();

        // Register all the calculated variables
        tr.registerDerivedVar("nMuons_CUT" + spec, nMuons);
        tr.registerDerivedVar("nElectrons_CUT" + spec, nElectrons);
        tr.registerDerivedVar("nIsoTrks_CUT" + spec, nIsoTrks);

        tr.registerDerivedVar("cntNJetsPt50Eta24" + spec, cntNJetsPt50Eta24);
        tr.registerDerivedVar("cntNJetsPt30Eta24" + spec, cntNJetsPt30Eta24);

        tr.registerDerivedVec("dPhiVec" + spec, dPhiVec);

        tr.registerDerivedVar("cntCSVS" + spec, cntCSVS);

        tr.registerDerivedVec("jetsLVec_forTagger" + spec, jetsLVec_forTagger);
        tr.registerDerivedVec("recoJetsBtag_forTagger" + spec, recoJetsBtag_forTagger);

        tr.registerDerivedVar("cntNJetsPt30" + spec, cntNJetsPt30);

        tr.registerDerivedVar("passLeptVeto" + spec, passLeptVeto);
        tr.registerDerivedVar("passMuonVeto" + spec, passMuonVeto);
        tr.registerDerivedVar("passEleVeto" + spec, passEleVeto);
        tr.registerDerivedVar("passIsoTrkVeto" + spec, passIsoTrkVeto);
        tr.registerDerivedVar("passnJets" + spec, passnJets);
        tr.registerDerivedVar("passdPhis" + spec, passdPhis);
        tr.registerDerivedVar("passBJets" + spec, passBJets);
        tr.registerDerivedVar("passMET" + spec, passMET);
        tr.registerDerivedVar("passTagger" + spec, passTagger);
        tr.registerDerivedVar("passBaseline" + spec, passBaseline);
        tr.registerDerivedVar("passBaselineNoTag" + spec, passBaselineNoTag);
//        tr.registerDerivedVar("passNewCuts" + spec, passNewCuts);

        tr.registerDerivedVar("nTopCandSortedCnt" + spec, nTopCandSortedCnt);

        tr.registerDerivedVar("best_lept_brJet_MT" + spec, type3Ptr->best_lept_brJet_MT);
        tr.registerDerivedVar("best_had_brJet_MT" + spec, type3Ptr->best_had_brJet_MT);
        tr.registerDerivedVar("best_had_brJet_mTcomb" + spec, type3Ptr->best_had_brJet_mTcomb);
        tr.registerDerivedVar("best_had_brJet_MT2" + spec, type3Ptr->best_had_brJet_MT2);

        double HT = AnaFunctions::calcHT(tr.getVec<TLorentzVector>(jetVecLabel), AnaConsts::pt50Eta24Arr);
        tr.registerDerivedVar("HT" + spec, HT);

        if( debug ) std::cout<<"passBaseline : "<<passBaseline<<"  passBaseline : "<<passBaseline<<std::endl;

    } 

    void operator()(NTupleReader &tr)
    {
        passBaseline(tr);
    }
} blv;

void passBaselineFunc(NTupleReader &tr)
{
    blv(tr);
}

namespace stopFunctions
{
class CleanJets
{
 public:
  void operator()(NTupleReader& tr) {internalCleanJets(tr);}
  void setMuonIso(const std::string muIsoFlag)
  {
    if(muIsoFlag.compare("mini") == 0)
      {
	muIsoStr_ = "muonsMiniIso";
	muIsoReq_ = AnaConsts::muonsMiniIsoArr;
      }
    else if(muIsoFlag.compare("rel") == 0)
      {
	muIsoStr_ = "muonsRelIso";
	muIsoReq_ = AnaConsts::muonsArr;
      }
    else
      {
	std::cout << "cleanJets(...): muon iso mode not recognized!!! Using \"rel iso\" settings." << std::endl;
	muIsoStr_ = "muonsRelIso";
	muIsoReq_ = AnaConsts::muonsArr;
      }
  }
  void setElecIso(const std::string elecIsoFlag)
  {
    if(elecIsoFlag.compare("mini") == 0)
      {
	std::cout << "cleanJets(...): electron mini iso mode not implemented yet!!! Using \"rel iso\" settings." << std::endl;
	//elecIsoStr = "elesMiniIso";
	//elecIsoReq = AnaConsts::elessMiniIsoArr;
	elecIsoStr_ = "elesRelIso";
	elecIsoReq_ = AnaConsts::elesArr;
      }
    else if(elecIsoFlag.compare("rel") == 0)
      {
	elecIsoStr_ = "elesRelIso";
	elecIsoReq_ = AnaConsts::elesArr;
      }
    else
      {
	std::cout << "cleanJets(...): muon iso mode not recognized!!! Using \"rel iso\" settings." << std::endl;
	elecIsoStr_ = "elesRelIso";
	elecIsoReq_ = AnaConsts::elesArr;
      }
  }
  CleanJets()
    {
      setMuonIso("rel");
      setElecIso("rel");
    }
 private:
  std::string muIsoStr_, elecIsoStr_;
  AnaConsts::IsoAccRec muIsoReq_;
  AnaConsts::ElecIsoAccRec elecIsoReq_;
  void internalCleanJets(NTupleReader& tr)
  {
    const std::vector<TLorentzVector>& jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
    const std::vector<TLorentzVector>& elesLVec = tr.getVec<TLorentzVector>("elesLVec");
    const std::vector<TLorentzVector>& muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
    const std::vector<double>& elesIso = tr.getVec<double>(elecIsoStr_);
    const std::vector<double>& muonsIso = tr.getVec<double>(muIsoStr_);
    const std::vector<double>& recoJetsBtag_0 = tr.getVec<double>("recoJetsBtag_0");
    const std::vector<int>& muMatchedJetIdx = tr.getVec<int>("muMatchedJetIdx");
    const std::vector<int>& eleMatchedJetIdx = tr.getVec<int>("eleMatchedJetIdx");
    const std::vector<unsigned int>& elesisEB = tr.getVec<unsigned int>("elesisEB");
    if(elesLVec.size() != elesIso.size()
       || elesLVec.size() != eleMatchedJetIdx.size()
       || elesLVec.size() != elesisEB.size()
       || muonsLVec.size() != muonsIso.size()
       || muonsLVec.size() != muMatchedJetIdx.size()
       || jetsLVec.size() != recoJetsBtag_0.size())
      {
	std::cout << "MISMATCH IN VECTOR SIZE!!!!! Aborting jet cleaning algorithm!!!!!!" << std::endl;
	return;
      }
    std::vector<TLorentzVector>* cleanJetVec = new std::vector<TLorentzVector>();
    std::vector<double>* cleanJetBTag = new std::vector<double>;
    std::vector<TLorentzVector>* cleanJetpt30ArrVec = new std::vector<TLorentzVector>();
    std::vector<double>* cleanJetpt30ArrBTag = new std::vector<double>;
    std::vector<int>* rejectJetIdx_formuVec = new std::vector<int>();
    std::vector<int>* rejectJetIdx_foreleVec = new std::vector<int>();
    const double jldRMax = 0.15;
    const double HT_jetPtMin = 50;
    const double HT_jetEtaMax = 2.4;
    const double MTH_jetPtMin = 30.0;
    double HT = 0.0, HTNoIso = 0.0;
    TLorentzVector MHT;
    std::vector<bool> keepJetPFCandMatch(jetsLVec.size(), true);
    for(int iM = 0; iM < muonsLVec.size() && iM < muonsIso.size() && iM < muMatchedJetIdx.size(); ++iM)
      {
	if(!AnaFunctions::passMuon(muonsLVec[iM], muonsIso[iM], 0.0, muIsoReq_)){ rejectJetIdx_formuVec->push_back(-1); continue; }
	if(muMatchedJetIdx[iM] >= 0){ keepJetPFCandMatch[muMatchedJetIdx[iM]] = false; rejectJetIdx_formuVec->push_back(muMatchedJetIdx[iM]); }
	else
	  {
	    //If muon matching to PF candidate has failed, use dR matching as fallback
	    int match = AnaFunctions::jetLepdRMatch(muonsLVec[iM], jetsLVec, jldRMax);
	    if(match >= 0){ keepJetPFCandMatch[match] = false; rejectJetIdx_formuVec->push_back(match); }
	    else rejectJetIdx_formuVec->push_back(-1);
	  }
      }
    for(int iE = 0; iE < elesLVec.size() && iE < elesIso.size() && iE < eleMatchedJetIdx.size(); ++iE)
      {
	if(!AnaFunctions::passElectron(elesLVec[iE], elesIso[iE], 0.0, elesisEB[iE], elecIsoReq_)){ rejectJetIdx_foreleVec->push_back(-1); continue; }
	if(eleMatchedJetIdx[iE] >= 0){ keepJetPFCandMatch[eleMatchedJetIdx[iE]] = false; rejectJetIdx_foreleVec->push_back(eleMatchedJetIdx[iE]); }
	else
	  {
	    //If electron matching to PF candidate has failed, use dR matching as fallback
	    int match = AnaFunctions::jetLepdRMatch(elesLVec[iE], jetsLVec, jldRMax);
	    if(match >= 0){ keepJetPFCandMatch[match] = false; rejectJetIdx_foreleVec->push_back(match); }
	    else rejectJetIdx_foreleVec->push_back(-1);
	  }
      }
    int jetsKept = 0;
    for(int iJet = 0; iJet < jetsLVec.size(); ++iJet)
      {
	if(keepJetPFCandMatch[iJet])
	  {
	    ++jetsKept;
	    cleanJetVec->push_back(jetsLVec[iJet]);
	    cleanJetBTag->push_back(recoJetsBtag_0[iJet]);
	    if(AnaFunctions::jetPassCuts(jetsLVec[iJet], AnaConsts::pt30Arr))
	      {
		cleanJetpt30ArrVec->push_back(jetsLVec[iJet]);
		cleanJetpt30ArrBTag->push_back(recoJetsBtag_0[iJet]);
	      }
	    if(jetsLVec[iJet].Pt() > HT_jetPtMin && fabs(jetsLVec[iJet].Eta()) < HT_jetEtaMax) HT += jetsLVec[iJet].Pt();
	    if(jetsLVec[iJet].Pt() > MTH_jetPtMin) MHT += jetsLVec[iJet];
	  }
      }
    tr.registerDerivedVar("nJetsRemoved", static_cast<int>(jetsLVec.size() - jetsKept));
    tr.registerDerivedVar("cleanHt", HT);
    tr.registerDerivedVar("cleanMHt", MHT.Pt());
    tr.registerDerivedVar("cleanMHtPhi", MHT.Phi());
    tr.registerDerivedVec("cleanJetVec", cleanJetVec);
    tr.registerDerivedVec("cleanJetBTag", cleanJetBTag);
    tr.registerDerivedVec("cleanJetpt30ArrVec", cleanJetpt30ArrVec);
    tr.registerDerivedVec("cleanJetpt30ArrBTag", cleanJetpt30ArrBTag);
    tr.registerDerivedVec("rejectJetIdx_formuVec", rejectJetIdx_formuVec);
    tr.registerDerivedVec("rejectJetIdx_foreleVec", rejectJetIdx_foreleVec);
  }
} cjh;

    void cleanJets(NTupleReader& tr)
    {
        cjh(tr);
    }
}

#endif
