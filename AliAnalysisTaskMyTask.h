/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskMyTask_H
#define AliAnalysisTaskMyTask_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMyTask : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskMyTask();
                                AliAnalysisTaskMyTask(const char *name);
        virtual                 ~AliAnalysisTaskMyTask();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        TString                 fPeriod;
        TString                 fTrigger;
        double                  fDcaXYCut;
        double                  fDcaZCut;

        AliAODEvent*            fAOD;
        AliAODVertex*           fVtx;
        AliPIDResponse*         fPIDResponse;
        TList*                  fOutputList;

        TH1F*                   fHistDcaXY[2][10];
        TH1F*                   fHistDcaZ[2][10];
        TH1F*                   fHistPt[2][10];
        TH1F*                   fHistEta[2][10];
        TH1F*                   fHistPhi[2][10];
        TH1F*                   fHistChi2[2][10];
        TH1F*                   fHistNCl[2][10];
        TH2F*                   fHistSigmaPionTPCTOF[2][10];
        TH2F*                   fHistSigmaKaonTPCTOF[2][10];
        TH2F*                   fHistSigmaProtonTPCTOF[2][10];
        TH2F*                   fHistSigmaDeuteronTPCTOF[2][10];
        TH2F*                   fHistSigmaTritonTPCTOF[2][10];
        TH2F*                   fHistSigmaHe3TPCTOF[2][10];
        TH2F*                   fHistPtVsDCAXY[2][10];
        TH2F*                   fHistPtVsDCAZ[2][10];

	//\ dowang
        TH2F*                   fHistTOFMassVsPt[2][10];
        TH2F*                   fHistdEdxVsPt[2][10];

        TH2F*                   fHistSigmaPionTPCVsPt[2][10];
        TH2F*                   fHistSigmaKaonTPCVsPt[2][10];
	TH2F*                   fHistSigmaProtonTPCVsPt[2][10];
        TH2F*                   fHistSigmaDeuteronTPCVsPt[2][10];
        TH2F*                   fHistSigmaTritonTPCVsPt[2][10];
        TH2F*                   fHistSigmaHe3TPCVsPt[2][10];

        TH2F*                   fHistSigmaPionTOFVsPt[2][10];
        TH2F*                   fHistSigmaKaonTOFVsPt[2][10];
	TH2F*                   fHistSigmaProtonTOFVsPt[2][10];
        TH2F*                   fHistSigmaDeuteronTOFVsPt[2][10];
        TH2F*                   fHistSigmaTritonTOFVsPt[2][10];
        TH2F*                   fHistSigmaHe3TOFVsPt[2][10];

        AliAnalysisTaskMyTask(const AliAnalysisTaskMyTask&); // not implemented
        AliAnalysisTaskMyTask& operator=(const AliAnalysisTaskMyTask&); // not implemented

        ClassDef(AliAnalysisTaskMyTask, 1);
};

#endif
