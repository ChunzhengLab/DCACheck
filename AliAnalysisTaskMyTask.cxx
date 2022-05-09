/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TMath.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskMyTask.h"
#include "AliMultSelection.h"

class AliAnalysisTaskMyTask;    // your analysis class
ClassImp(AliAnalysisTaskMyTask) // classimp: necessary for root

AliAnalysisTaskMyTask::AliAnalysisTaskMyTask() : AliAnalysisTaskSE(),
  fPeriod("LHC10h"),
  fTrigger("kMB"),
  fDcaXYCut(2.4),
  fDcaZCut(3.2),
  fAOD(nullptr),
  fVtx(nullptr),
  fPIDResponse(nullptr),
  fOutputList(nullptr)
{
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 10; j++) {
      fHistDcaXY[i][j]               = nullptr;
      fHistDcaZ[i][j]                = nullptr;
      fHistPt[i][j]                  = nullptr;
      fHistEta[i][j]                 = nullptr;
      fHistPhi[i][j]                 = nullptr;
      fHistChi2[i][j]                = nullptr;
      fHistNCl[i][j]                 = nullptr;
      fHistSigmaPionTPCTOF[i][j]     = nullptr;
      fHistSigmaKaonTPCTOF[i][j]     = nullptr;
      fHistSigmaProtonTPCTOF[i][j]   = nullptr;
      fHistSigmaDeuteronTPCTOF[i][j] = nullptr;
      fHistSigmaTritonTPCTOF[i][j]   = nullptr;
      fHistSigmaHe3TPCTOF[i][j]      = nullptr;
      fHistPtVsDCAXY[i][j]           = nullptr;
      fHistPtVsDCAZ[i][j]            = nullptr;
    }
  }
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::AliAnalysisTaskMyTask(const char* name) : AliAnalysisTaskSE(name),
  fPeriod("LHC10h"),
  fTrigger("kMB"),
  fDcaXYCut(2.4),
  fDcaZCut(3.2),
  fAOD(nullptr),
  fVtx(nullptr),
  fPIDResponse(nullptr),
  fOutputList(nullptr)
{
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 10; j++) {
      fHistDcaXY[i][j]               = nullptr;
      fHistDcaZ[i][j]                = nullptr;
      fHistPt[i][j]                  = nullptr;
      fHistEta[i][j]                 = nullptr;
      fHistPhi[i][j]                 = nullptr;
      fHistChi2[i][j]                = nullptr;
      fHistNCl[i][j]                 = nullptr;
      fHistSigmaPionTPCTOF[i][j]     = nullptr;
      fHistSigmaKaonTPCTOF[i][j]     = nullptr;
      fHistSigmaProtonTPCTOF[i][j]   = nullptr;
      fHistSigmaDeuteronTPCTOF[i][j] = nullptr;
      fHistSigmaTritonTPCTOF[i][j]   = nullptr;
      fHistSigmaHe3TPCTOF[i][j]      = nullptr;
      fHistPtVsDCAXY[i][j]           = nullptr;
      fHistPtVsDCAZ[i][j]            = nullptr;
    }
  }
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::~AliAnalysisTaskMyTask()
{
    // destructor
    if(fOutputList) delete fOutputList;
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  //Histogram
  char charDCACut[20];
  char charFB[20];
  for (int i = 0; i < 2; i++) {
  if (i==0) sprintf(charDCACut, "bfDCACut");
  if (i==1) sprintf(charDCACut, "afDCACut");
    for (int jFB = 0; jFB < 10; jFB++) {
      if(jFB == 0) sprintf(charFB, "FB1");
      if(jFB == 1) sprintf(charFB, "FB16");
      if(jFB == 2) sprintf(charFB, "FB32");
      if(jFB == 3) sprintf(charFB, "FB64");
      if(jFB == 4) sprintf(charFB, "FB96");
      if(jFB == 5) sprintf(charFB, "FB128");
      if(jFB == 6) sprintf(charFB, "FB256");
      if(jFB == 7) sprintf(charFB, "FB272");
      if(jFB == 8) sprintf(charFB, "FB512");
      if(jFB == 9) sprintf(charFB, "FB768");

      fHistDcaXY[i][jFB]               = new TH1F(Form("fHistDcaXY_%s_%s", charFB, charDCACut), Form("DcaXY_%s_%s", charFB, charDCACut), 1000, -10, 10);
      fHistDcaZ[i][jFB]                = new TH1F(Form("fHistDcaZ_%s_%s", charFB, charDCACut), Form("DcaZ_%s_%s", charFB, charDCACut), 1000, -10, 10);
      fHistPt[i][jFB]                  = new TH1F(Form("fHistPt_%s_%s", charFB, charDCACut), Form("pT_%s_%s", charFB, charDCACut), 200, 0, 10);
      fHistEta[i][jFB]                 = new TH1F(Form("fHistEta_%s_%s", charFB, charDCACut), Form("Eta_%s_%s", charFB, charDCACut), 100, -3, 3);
      fHistPhi[i][jFB]                 = new TH1F(Form("fHistPhi_%s_%s", charFB, charDCACut), Form("Phi_%s_%s", charFB, charDCACut), 360, 0, TMath::TwoPi());
      fHistChi2[i][jFB]                = new TH1F(Form("fHistChi2_%s_%s", charFB, charDCACut), Form("Chi2_%s_%s", charFB, charDCACut), 200, -1,10);
      fHistNCl[i][jFB]                 = new TH1F(Form("fHistNCl_%s_%s", charFB, charDCACut), Form("NCl_%s_%s", charFB, charDCACut), 200, 0, 200.);
      fHistSigmaPionTPCTOF[i][jFB]     = new TH2F(Form("fHist2DnSigmaPionTPCTOF_%s_%s", charFB, charDCACut), Form("nSigmaPionTPC vs. TOF_%s_%s", charFB, charDCACut), 500,-50,50,500,-50,50);
      fHistSigmaKaonTPCTOF[i][jFB]     = new TH2F(Form("fHist2DnSigmaKaonTPCTOF_%s_%s", charFB, charDCACut), Form("nSigmaKaonTPC vs. TOF_%s_%s", charFB, charDCACut), 500,-50,50,500,-50,50);
      fHistSigmaProtonTPCTOF[i][jFB]   = new TH2F(Form("fHist2DnSigmaProtonTPCTOF_%s_%s", charFB, charDCACut), Form("nSigmaProton vs. TOF_%s_%s", charFB, charDCACut), 500,-50,50,500,-50,50);
      fHistSigmaDeuteronTPCTOF[i][jFB] = new TH2F(Form("fHist2DnSigmaDeuteronTPCTOF_%s_%s", charFB, charDCACut), Form("nSigmaDeuteronTPC vs. TOF_%s_%s", charFB, charDCACut), 500,-50,50,500,-50,50);
      fHistSigmaTritonTPCTOF[i][jFB]   = new TH2F(Form("fHist2DnSigmaTritonTPCTOF_%s_%s", charFB, charDCACut), Form("nSigmaTritonTPC vs. TOF_%s_%s", charFB, charDCACut), 500,-50,50,500,-50,50);
      fHistSigmaHe3TPCTOF[i][jFB]      = new TH2F(Form("fHist2DnSigmaHe3TPCTOF_%s_%s", charFB, charDCACut), Form("nSigmaHe3Proton TPC vs. TOF_%s_%s", charFB, charDCACut), 500,-50,50,500,-50,50);
      fHistPtVsDCAXY[i][jFB]           = new TH2F(Form("fHist2DPtVsDCAXY_%s_%s", charFB, charDCACut), Form("pT Vs DCAXY_%s_%s", charFB, charDCACut), 200,0,10,200,-10,10);
      fHistPtVsDCAZ[i][jFB]            = new TH2F(Form("fHist2DPtVsDCAZ_%s_%s", charFB, charDCACut), Form("pT Vs DCAZ_%s_%s", charFB, charDCACut), 200,0,10,200,-10,10);
      fOutputList -> Add(fHistDcaXY[i][jFB]);
      fOutputList -> Add(fHistDcaZ[i][jFB]);
      fOutputList -> Add(fHistPt[i][jFB]);
      fOutputList -> Add(fHistEta[i][jFB]);
      fOutputList -> Add(fHistPhi[i][jFB]);
      fOutputList -> Add(fHistChi2[i][jFB]);
      fOutputList -> Add(fHistNCl[i][jFB]);
      fOutputList -> Add(fHistSigmaPionTPCTOF[i][jFB]);
      fOutputList -> Add(fHistSigmaKaonTPCTOF[i][jFB]);
      fOutputList -> Add(fHistSigmaProtonTPCTOF[i][jFB]);
      fOutputList -> Add(fHistSigmaDeuteronTPCTOF[i][jFB]);
      fOutputList -> Add(fHistSigmaTritonTPCTOF[i][jFB]);
      fOutputList -> Add(fHistSigmaHe3TPCTOF[i][jFB]);
      fOutputList -> Add(fHistPtVsDCAXY[i][jFB]);
      fOutputList -> Add(fHistPtVsDCAZ[i][jFB]);
    }
  }
  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserExec(Option_t *)
{
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) return;
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  AliAnalysisManager *fMan = AliAnalysisManager::GetAnalysisManager();
  if (!fMan) return;
  AliInputEventHandler* fHandler = (AliInputEventHandler*)(fMan->GetInputEventHandler());
  if (!fHandler) return;
  fPIDResponse = fHandler->GetPIDResponse();

  //----------------------------
  // Trigger
  //----------------------------
  UInt_t mask = fHandler->IsEventSelected();
  bool isTrigselected = false;
  if (fTrigger.EqualTo("kMB"))
  isTrigselected = mask & AliVEvent::kMB;
  else if (fTrigger.EqualTo("kINT7"))
  isTrigselected = mask & AliVEvent::kINT7;
  else if (fTrigger.EqualTo("kINT7+kCentral+kSemiCentral"))
  isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kCentral + AliVEvent::kSemiCentral);
  if (isTrigselected == false) return;

  //----------------------------
  // Vertex
  //----------------------------
  double vtx[3] = {0};
  fVtx -> GetXYZ(vtx);
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  if (fabs(vtx[0]) < 1e-6 || fabs(vtx[1]) < 1e-6 || fabs(vtx[2]) < 1e-6) return;
  double dz = vtx[2] - fAOD->GetPrimaryVertexSPD()->GetZ();
  if (fabs(vtx[2]) > 10) return;
  if (!fVtx || fVtx -> GetNContributors() < 2 || vtSPD -> GetNContributors() < 1) return;
  if (fPeriod.EqualTo("LHC10h")) if (fabs(dz) > 0.5) return;
  if (fPeriod.EqualTo("LHC15o")) {
      double covTrc[6], covSPD[6];
      fVtx -> GetCovarianceMatrix(covTrc);
      fAOD -> GetPrimaryVertexSPD() -> GetCovarianceMatrix(covSPD);
      double errTot = TMath::Sqrt(covTrc[5] + covSPD[5]);
      double errTrc = TMath::Sqrt(covTrc[5]);
      double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
      if (fabs(dz) > 0.2 || nsigTot > 10 || nsigTrc > 20) return;
  }
  
  //----------------------------
  // Centrality
  //----------------------------
  double centV0M = -999;
  double centCL1 = -999;
  if (fPeriod.EqualTo("LHC10h")) {
    centV0M  = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
    centCL1  = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
  } else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    centV0M  = fMultSel->GetMultiplicityPercentile("V0M");
    centCL1 = fMultSel->GetMultiplicityPercentile("CL1");
  }
  if (fabs(centV0M - centCL1) > 7.5) return;
  if (centV0M < 0 || centV0M > 100) return;

  //----------------------------
  // Track loop
  //----------------------------
  int nTracks(fAOD->GetNumberOfTracks());
  for(int iTrack(0); iTrack < nTracks; iTrack++) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track) continue;

    double mag = fAOD->GetMagneticField();
    double dca[2] = {-999,-999};
    double cov[3] = {-999,-999,-999};
    if (!track->PropagateToDCA(fVtx, mag, 100., dca, cov)) continue;
    double pt             = track->Pt();
    double eta            = track->Eta();
    double phi            = track->Phi();
    short  charge         = track->Charge();
    double dedx           = track->GetTPCsignal();
    double chi2           = track->Chi2perNDF();
    unsigned short nhits  = track->GetTPCNcls();


    float nSigmaTPCPion     = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kPion);
    float nSigmaTPCKaon     = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kKaon);
    float nSigmaTPCProton   = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kProton);
    float nSigmaTPCDeuteron = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kDeuteron);
    float nSigmaTPCTriton   = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kTriton);
    float nSigmaTPCHe3      = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kHe3);

    float nSigmaTOFPion     = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kPion);
    float nSigmaTOFKaon     = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kKaon);
    float nSigmaTOFProton   = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kProton);
    float nSigmaTOFDeuteron = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kDeuteron);
    float nSigmaTOFTriton   = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kTriton);
    float nSigmaTOFHe3      = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kHe3);

    int nFBBits = 10;
    TBits FBBits(nFBBits);
    FBBits.SetBitNumber(0, track->TestFilterBit(1));
    FBBits.SetBitNumber(1, track->TestFilterBit(16));
    FBBits.SetBitNumber(2, track->TestFilterBit(32));
    FBBits.SetBitNumber(3, track->TestFilterBit(64));
    FBBits.SetBitNumber(4, track->TestFilterBit(96));
    FBBits.SetBitNumber(5, track->TestFilterBit(128));
    FBBits.SetBitNumber(6, track->TestFilterBit(256));
    FBBits.SetBitNumber(7, track->TestFilterBit(272));
    FBBits.SetBitNumber(8, track->TestFilterBit(512));
    FBBits.SetBitNumber(9, track->TestFilterBit(768));




    for (int iFBBits = 0; iFBBits < nFBBits; iFBBits++) {
      if (FBBits.TestBitNumber(iFBBits)) {
        fHistDcaXY[0][iFBBits]              -> Fill(dca[0]);
        fHistDcaZ[0][iFBBits]               -> Fill(dca[1]);
        fHistPt[0][iFBBits]                 -> Fill(pt);
        fHistEta[0][iFBBits]                -> Fill(eta);
        fHistPhi[0][iFBBits]                -> Fill(phi);
        fHistChi2[0][iFBBits]               -> Fill(chi2);
        fHistNCl[0][iFBBits]                -> Fill(nhits);
        fHistSigmaPionTPCTOF[0][iFBBits]    -> Fill(nSigmaTPCPion, nSigmaTOFPion);
        fHistSigmaKaonTPCTOF[0][iFBBits]    -> Fill(nSigmaTPCKaon, nSigmaTOFKaon);
        fHistSigmaProtonTPCTOF[0][iFBBits]  -> Fill(nSigmaTPCProton, nSigmaTOFProton);   
        fHistSigmaDeuteronTPCTOF[0][iFBBits]-> Fill(nSigmaTPCDeuteron, nSigmaTOFDeuteron);
        fHistSigmaTritonTPCTOF[0][iFBBits]  -> Fill(nSigmaTPCTriton, nSigmaTOFTriton);
        fHistSigmaHe3TPCTOF[0][iFBBits]     -> Fill(nSigmaTPCHe3, nSigmaTOFHe3);
        fHistPtVsDCAXY[0][iFBBits]          -> Fill(pt, dca[0]);
        fHistPtVsDCAZ[0][iFBBits]           -> Fill(pt, dca[1]);
      }
    }

    if (fabs(dca[0]) > fDcaXYCut) continue;
    if (fabs(dca[1]) > fDcaZCut)  continue;

    for (int iFBBits = 0; iFBBits < nFBBits; iFBBits++) {
      if (FBBits.TestBitNumber(iFBBits)) {
        fHistDcaXY[1][iFBBits]              -> Fill(dca[0]);
        fHistDcaZ[1][iFBBits]               -> Fill(dca[1]);
        fHistPt[1][iFBBits]                 -> Fill(pt);
        fHistEta[1][iFBBits]                -> Fill(eta);
        fHistPhi[1][iFBBits]                -> Fill(phi);
        fHistChi2[1][iFBBits]               -> Fill(chi2);
        fHistNCl[1][iFBBits]                -> Fill(nhits);
        fHistSigmaPionTPCTOF[1][iFBBits]    -> Fill(nSigmaTPCPion, nSigmaTOFPion);
        fHistSigmaKaonTPCTOF[1][iFBBits]    -> Fill(nSigmaTPCKaon, nSigmaTOFKaon);
        fHistSigmaProtonTPCTOF[1][iFBBits]  -> Fill(nSigmaTPCProton, nSigmaTOFProton);   
        fHistSigmaDeuteronTPCTOF[1][iFBBits]-> Fill(nSigmaTPCDeuteron, nSigmaTOFDeuteron);
        fHistSigmaTritonTPCTOF[1][iFBBits]  -> Fill(nSigmaTPCTriton, nSigmaTOFTriton);
        fHistSigmaHe3TPCTOF[1][iFBBits]     -> Fill(nSigmaTPCHe3, nSigmaTOFHe3);
        fHistPtVsDCAXY[1][iFBBits]          -> Fill(pt, dca[0]);
        fHistPtVsDCAZ[1][iFBBits]           -> Fill(pt, dca[1]);
      }
    }
  }

  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________