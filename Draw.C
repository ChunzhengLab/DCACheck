#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectoryFile.h"
#include "TCanvas.h"
#include "TPad.h"

template<typename TH> 
void DrawHist(TH *hist, TCanvas* c, int bLogy, bool bColz);

void Draw() {
  TFile* inputFile = TFile::Open("AnalysisResults_10h_139310.root");
  TDirectoryFile *inputFolder = (TDirectoryFile*)inputFile->Get("MyTask");

  TList* inputList = nullptr;
  inputFolder->GetObject("MyOutputContainer", inputList);

  TH1F* fHistDcaXY[2][10];
  TH1F* fHistDcaZ[2][10];
  TH1F* fHistPt[2][10];
  TH1F* fHistEta[2][10];
  TH1F* fHistPhi[2][10];
  TH1F* fHistChi2[2][10];
  TH1F* fHistNCl[2][10];
  TH2F* fHist2DnSigmaPionTPCTOF[2][10];
  TH2F* fHist2DnSigmaKaonTPCTOF[2][10];
  TH2F* fHist2DnSigmaProtonTPCTOF[2][10];
  TH2F* fHist2DnSigmaDeuteronTPCTOF[2][10];
  TH2F* fHist2DnSigmaTritonTPCTOF[2][10];
  TH2F* fHist2DnSigmaHe3TPCTOF[2][10];
  TH2F* fHist2DPtVsDCAXY[2][10];
  TH2F* fHist2DPtVsDCAZ[2][10];

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

      fHistDcaXY[i][jFB] = (TH1F*)inputList->FindObject(Form("fHistDcaXY_%s_%s", charFB, charDCACut));
      fHistDcaZ[i][jFB] = (TH1F*)inputList->FindObject(Form("fHistDcaZ_%s_%s", charFB, charDCACut));
      fHistPt[i][jFB] = (TH1F*)inputList->FindObject(Form("fHistPt_%s_%s", charFB, charDCACut));
      fHistEta[i][jFB] = (TH1F*)inputList->FindObject(Form("fHistEta_%s_%s", charFB, charDCACut));
      fHistPhi[i][jFB] = (TH1F*)inputList->FindObject(Form("fHistPhi_%s_%s", charFB, charDCACut));
      fHistChi2[i][jFB] = (TH1F*)inputList->FindObject(Form("fHistChi2_%s_%s", charFB, charDCACut));
      fHistNCl[i][jFB] = (TH1F*)inputList->FindObject(Form("fHistNCl_%s_%s", charFB, charDCACut));
      fHist2DnSigmaPionTPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaPionTPCTOF_%s_%s", charFB, charDCACut));
      fHist2DnSigmaKaonTPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaKaonTPCTOF_%s_%s", charFB, charDCACut));
      fHist2DnSigmaProtonTPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaProtonTPCTOF_%s_%s", charFB, charDCACut));
      fHist2DnSigmaDeuteronTPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaDeuteronTPCTOF_%s_%s", charFB, charDCACut));
      fHist2DnSigmaTritonTPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaTritonTPCTOF_%s_%s", charFB, charDCACut));
      fHist2DnSigmaHe3TPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaHe3TPCTOF_%s_%s", charFB, charDCACut));
      fHist2DPtVsDCAXY[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DPtVsDCAXY_%s_%s", charFB, charDCACut));
      fHist2DPtVsDCAZ[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DPtVsDCAZ_%s_%s", charFB, charDCACut));

      if(i == 1) {
        fHistDcaXY[i][jFB] ->SetLineColor(kRed);
        fHistDcaZ[i][jFB] ->SetLineColor(kRed);
        fHistPt[i][jFB] ->SetLineColor(kRed);
        fHistEta[i][jFB] ->SetLineColor(kRed);
        fHistPhi[i][jFB] ->SetLineColor(kRed);
        fHistChi2[i][jFB] ->SetLineColor(kRed);
        fHistNCl[i][jFB] ->SetLineColor(kRed);
      }

      fHist2DnSigmaPionTPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaPionTPCTOF_%s_%s", charFB, charDCACut));
      fHist2DnSigmaKaonTPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaKaonTPCTOF_%s_%s", charFB, charDCACut));
      fHist2DnSigmaProtonTPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaProtonTPCTOF_%s_%s", charFB, charDCACut));
      fHist2DnSigmaDeuteronTPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaDeuteronTPCTOF_%s_%s", charFB, charDCACut));
      fHist2DnSigmaTritonTPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaTritonTPCTOF_%s_%s", charFB, charDCACut));
      fHist2DnSigmaHe3TPCTOF[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DnSigmaHe3TPCTOF_%s_%s", charFB, charDCACut));
      fHist2DPtVsDCAXY[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DPtVsDCAXY_%s_%s", charFB, charDCACut));
      fHist2DPtVsDCAZ[i][jFB] = (TH2F*)inputList->FindObject(Form("fHist2DPtVsDCAZ_%s_%s", charFB, charDCACut));

      fHistDcaXY[i][jFB] ->GetXaxis()->SetRangeUser(-5, 5);
      fHistDcaZ[i][jFB] ->GetXaxis()->SetRangeUser(-5, 5);
      fHistPt[i][jFB] ->GetXaxis()->SetRangeUser(0, 5);
      fHistPhi[i][jFB] -> SetMinimum(0);
      fHistPhi[i][jFB] -> SetMaximum(3e5);

      //这里可以写入DCA的Cut;
    }
  }

  TCanvas* cHistDcaXY                  = new TCanvas("HistDcaXY"                 ,"HistDcaXY"                 ,3000,300); 
  TCanvas* cHistDcaZ                   = new TCanvas("HistDcaZ"                  ,"HistDcaZ"                  ,3000,300);
  TCanvas* cHistPt                     = new TCanvas("HistPt"                    ,"HistPt"                    ,3000,300);
  TCanvas* cHistEta                    = new TCanvas("HistEta"                   ,"HistEta"                   ,3000,300);
  TCanvas* cHistPhi                    = new TCanvas("HistPhi"                   ,"HistPhi"                   ,3000,300);
  TCanvas* cHistChi2                   = new TCanvas("HistChi2"                  ,"HistChi2"                  ,3000,300);
  TCanvas* cHistNCl                    = new TCanvas("HistNCl"                   ,"HistNCl"                   ,3000,300);
  // TCanvas* cHist2DnSigmaPionTPCTOF     = new TCanvas("Hist2DnSigmaPionTPCTOF"    ,"Hist2DnSigmaPionTPCTOF"    ,3000,300);
  // TCanvas* cHist2DnSigmaKaonTPCTOF     = new TCanvas("Hist2DnSigmaKaonTPCTOF"    ,"Hist2DnSigmaKaonTPCTOF"    ,3000,300);
  // TCanvas* cHist2DnSigmaProtonTPCTOF   = new TCanvas("Hist2DnSigmaProtonTPCTOF"  ,"Hist2DnSigmaProtonTPCTOF"  ,3000,300);
  // TCanvas* cHist2DnSigmaDeuteronTPCTOF = new TCanvas("Hist2DnSigmaDeuteronTPCTOF","Hist2DnSigmaDeuteronTPCTOF",3000,300);
  // TCanvas* cHist2DnSigmaTritonTPCTOF   = new TCanvas("Hist2DnSigmaTritonTPCTOF"  ,"Hist2DnSigmaTritonTPCTOF"  ,3000,300);
  // TCanvas* cHist2DnSigmaHe3TPCTOF      = new TCanvas("Hist2DnSigmaHe3TPCTOF"     ,"Hist2DnSigmaHe3TPCTOF"     ,3000,300);
  // TCanvas* cHist2DPtVsDCAXY            = new TCanvas("Hist2DPtVsDCAXY"           ,"Hist2DPtVsDCAXY"           ,3000,300);
  // TCanvas* cHist2DPtVsDCAZ             = new TCanvas("Hist2DPtVsDCAZ"            ,"Hist2DPtVsDCAZ"            ,3000,300);
  DrawHist(fHistDcaXY, cHistDcaXY, 1, 0);
  DrawHist(fHistDcaZ, cHistDcaZ, 1, 0);
  DrawHist(fHistPt , cHistPt, 0, 0);
  DrawHist(fHistEta, cHistEta, 0, 0);
  DrawHist(fHistPhi, cHistPhi, 0, 0);
  DrawHist(fHistChi2, cHistChi2, 0, 0);
  DrawHist(fHistNCl , cHistNCl, 0, 0);
  // DrawHist(fHist2DnSigmaPionTPCTOF, cHist2DnSigmaPionTPCTOF, 0, 1);
  // DrawHist(fHist2DnSigmaKaonTPCTOF, cHist2DnSigmaKaonTPCTOF, 0, 1);
  // DrawHist(fHist2DnSigmaProtonTPCTOF, cHist2DnSigmaProtonTPCTOF, 0, 1);
  // DrawHist(fHist2DnSigmaDeuteronTPCTOF, cHist2DnSigmaDeuteronTPCTOF, 0, 1);
  // DrawHist(fHist2DnSigmaTritonTPCTOF, cHist2DnSigmaTritonTPCTOF, 0, 1);
  // DrawHist(fHist2DnSigmaHe3TPCTOF, cHist2DnSigmaHe3TPCTOF, 0, 1);
  // DrawHist(fHist2DPtVsDCAXY, cHist2DPtVsDCAXY, 0, 1);
  // DrawHist(fHist2DPtVsDCAZ, cHist2DPtVsDCAZ, 0, 1);
}


template<typename TH> 
void DrawHist(TH* hist, TCanvas* c, int bLogy, bool bColz) {
  c->Divide(10,2);
  int nPad = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 10; j++) {
      nPad++;
      if (bLogy) gPad->SetLogy();
      c->cd(nPad);
      if (bColz) hist[i][j]->Draw("colz");
      else hist[i][j]->Draw();
    }
  }
  TString name = c->GetName();
  c->SaveAs(name+".pdf");
}

