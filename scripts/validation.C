#include "TROOT.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChainElement.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"

#include <iostream>
#include <string>
#include <vector>

#include "/uscms/home/aperloff/Scripts/tdrstyle_mod14.C"

using namespace std;

//______________________________________________________________________________
void setChainElementNames(TChain* c, string dirTreeName) {
  TObjArray *fileElements=c->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
    chEl->SetName(dirTreeName.c_str());
  }
}

//______________________________________________________________________________
TCanvas* makeMatchCanvas(TH1F* match, TH1F* strictMatch) {
   TH1D* frame = new TH1D("frame", "frame;N_{match}(sim, emu);Events", 11, 0, 11);
   frame->GetYaxis()->SetRangeUser(0, strictMatch->GetMaximum()*1.1);
   TCanvas* canvas = tdrCanvas("NMatch", frame, 13, 11, true);
   match->SetMarkerStyle(kFullCircle);
   match->SetMarkerColor(kRed);
   match->SetLineColor(kRed);
   match->Draw("E1 same");
   strictMatch->SetMarkerStyle(kOpenCircle);
   strictMatch->SetMarkerColor(kRed-2);
   strictMatch->SetLineColor(kRed-2);
   strictMatch->Draw("E1 same");

   TLegend* l = tdrLeg(0.58,0.9,0.86,0.7);
   l->AddEntry((TObject*)0, "FastHisto","");
   l->AddEntry((TObject*)0, "#epsilon=0.10cm","");
   l->AddEntry((TObject*)0, "#delta=0.10cm","");
   l->AddEntry(match, "Total Matches", "pl");
   l->AddEntry(strictMatch, "Strict Ordering", "pl");
   return canvas;
}

//______________________________________________________________________________
TCanvas* makeSumPtCanvas(TH1F* reco, TH1F* emulation) {
   TH1D* frame = new TH1D("frame", "frame;#Sigma_{tracks #in vtx} p_{T};Fraction of Events / 5 GeV", 100, 0, 500);
   frame->GetYaxis()->SetRangeUser(0, 0.07);
   TCanvas* canvas = tdrCanvas("SumPt", frame, 13, 11, true);
   reco->SetMarkerStyle(kFullCircle);
   //emulation->SetMarkerStyle(kOpenCircle);
   reco->SetMarkerColor(kRed);
   //emulation->SetMarkerColor(kRed);
   reco->SetLineColor(kRed);
   //emulation->SetLineColor(kRed);
   reco->Scale(1.0/reco->Integral());
   //emulation->Scale(1.0/emulation->Integral());
   reco->Draw("E1 same");
   //emulation->Draw("E1 same");

   TLegend* l = tdrLeg(0.66,0.9,0.86,0.7);
   l->AddEntry((TObject*)0, "FastHisto","");
   l->AddEntry((TObject*)0, "#epsilon=0.10cm","");
   l->AddEntry((TObject*)0, "#delta=0.10cm","");
   l->AddEntry(reco, "Simulation", "pl");
   //l->AddEntry(emulation, "Emulation", "pl");
   return canvas;
}

//______________________________________________________________________________
TCanvas* makeNTrackCanvas(TH1F* reco, TH1F* emulation) {
   TH1D* frame = new TH1D("frame", "frame;N_{tracks};Entries", 80, 0, 80);
   frame->GetYaxis()->SetRangeUser(0, 0.1);
   TCanvas* canvas = tdrCanvas("NTracks", frame, 13, 11, true);
   reco->SetMarkerStyle(kFullCircle);
   //emulation->SetMarkerStyle(kOpenCircle);
   reco->SetMarkerColor(kRed);
   //emulation->SetMarkerColor(kRed);
   reco->SetLineColor(kRed);
   //emulation->SetLineColor(kRed);
   reco->Scale(1.0/reco->Integral());
   //emulation->Scale(1.0/emulation->Integral());
   reco->Draw("E1 same");
   //emulation->Draw("E1 same");

   TLegend* l = tdrLeg(0.66,0.9,0.86,0.7);
   l->AddEntry((TObject*)0, "FastHisto","");
   l->AddEntry((TObject*)0, "#epsilon=0.10cm","");
   l->AddEntry((TObject*)0, "#delta=0.10cm","");
   l->AddEntry(reco, "Simulation", "pl");
   //l->AddEntry(emulation, "Emulation", "pl");
   return canvas;
}

//______________________________________________________________________________
TCanvas* makeResidualCanvas(TH1F* reco, TH1F* emulation) {
   TH1D* frame = new TH1D("frame", "frame;z_{L1}-z_{gen} (cm);Fraction of Events / 0.025 cm", 80, -1, 1);
   frame->GetYaxis()->SetRangeUser(0, 0.16);
   TCanvas* canvas = tdrCanvas("VertexingResidual", frame, 13, 11, true);
   reco->SetMarkerStyle(kFullCircle);
   reco->SetMarkerColor(kRed);
   reco->SetLineColor(kRed);
   reco->Scale(1.0/reco->Integral());
   reco->Draw("E1 sames");
   emulation->SetMarkerStyle(kOpenCircle);
   emulation->SetMarkerColor(kRed-2);
   emulation->SetLineColor(kRed-2);
   emulation->Scale(1.0/emulation->Integral());
   emulation->Draw("E1 sames");

   TF1* reco_fit = new TF1("reco_residual_fit", "gaus", -1, 1);
   reco_fit->SetLineStyle(kSolid);
   reco->Fit(reco_fit);
   TF1* emulation_fit = new TF1("emulation_residual_fit", "gaus", -1, 1);
   emulation_fit->SetLineStyle(kDashed);
   emulation->Fit(emulation_fit);

   TLine* reco_line = new TLine(reco_fit->GetParameter(1),0,reco_fit->GetParameter(1),0.16);
   reco_line->SetLineColor(kRed);
   reco_line->SetLineStyle(kSolid);
   reco_line->Draw("same");
   TLine* emulation_line = new TLine(emulation_fit->GetParameter(1),0,emulation_fit->GetParameter(1),0.16);
   emulation_line->SetLineColor(kRed-2);
   emulation_line->SetLineStyle(kDashed);
   emulation_line->Draw("same");

   gPad->Update();

   TLegend* l = tdrLeg(0.66,0.9,0.86,0.7);
   l->AddEntry((TObject*)0, "FastHisto","");
   l->AddEntry((TObject*)0, "#epsilon=0.10cm","");
   l->AddEntry((TObject*)0, "#delta=0.10cm","");
   l->AddEntry(reco, "Simulation", "pl");
   l->AddEntry(emulation, "Emulation", "pl");

   TLegend* stats = tdrLeg(0.55,0.68,0.86,0.53);
   stats->SetTextSize(0.02);
   stats->AddEntry((TObject*)0, Form("Reco\\space Mean \\hfill %0.3g\\pm%0.3g", reco_fit->GetParameter(1), reco_fit->GetParError(1)), "");
   stats->AddEntry((TObject*)0, Form("Reco\\space Sigma \\hfill %0.3f\\pm%0.3f", reco_fit->GetParameter(2), reco_fit->GetParError(2)), "");
   stats->AddEntry((TObject*)0, Form("Reco\\space\\chi^{2}/ndf \\hfill %0.2f/%i", reco_fit->GetChisquare(), reco_fit->GetNDF()), "");
   stats->AddEntry((TObject*)0, Form("Emulation\\space Mean \\hfill %0.3g\\pm%0.3g", emulation_fit->GetParameter(1), emulation_fit->GetParError(1)), "");
   stats->AddEntry((TObject*)0, Form("Emulation\\space Sigma \\hfill %0.3f\\pm%0.3f", emulation_fit->GetParameter(2), emulation_fit->GetParError(2)), "");
   stats->AddEntry((TObject*)0, Form("Emulation\\space\\chi^{2}/ndf \\hfill %0.2f/%i", emulation_fit->GetChisquare(), emulation_fit->GetNDF()), "");

   return canvas;
}

//______________________________________________________________________________
TCanvas* makeResolutionCanvas(TProfile* reco, TProfile* emulation) {
   TH1D* frame = new TH1D("frame", "frame;z_{gen} (cm);|z_{L1}-z_{gen}| (cm)", 60, -15, 15);
   frame->GetYaxis()->SetRangeUser(0, 0.1);
   TCanvas* canvas = tdrCanvas("VertexResolution", frame, 13, 11, true);
   reco->SetMarkerStyle(kFullCircle);
   emulation->SetMarkerStyle(kOpenCircle);
   reco->SetMarkerColor(kRed);
   emulation->SetMarkerColor(kRed-2);
   reco->SetLineColor(kRed);
   emulation->SetLineColor(kRed-2);
   reco->Draw("same");
   emulation->Draw("same");

   TF1* reco_fit = new TF1("reco_resolution_fit", "pol0", -15, 15);
   TF1* emulation_fit = new TF1("emulation_resolution_fit", "pol0", -15, 15);
   reco_fit->SetLineStyle(kSolid);
   emulation_fit->SetLineStyle(kDashed);
   reco->Fit(reco_fit);
   emulation->Fit(emulation_fit);

   TLegend* l = tdrLeg(0.4,0.35,0.6,0.15);
   l->AddEntry((TObject*)0, "FastHisto","");
   l->AddEntry((TObject*)0, "#epsilon=0.10cm","");
   l->AddEntry((TObject*)0, "#delta=0.10cm","");
   l->AddEntry(reco, "Simulation", "pl");
   l->AddEntry(emulation, "Emulation", "pl");
   return canvas;
}

//______________________________________________________________________________
TCanvas* makeEfficiencyCanvas(TEfficiency* reco, TEfficiency* emulation) {
   TH1D* frame = new TH1D("frame", "frame;z_{gen} (cm);Vertex Reconstruction Efficiency", 30, -15, 15);
   frame->GetYaxis()->SetRangeUser(0, 1.1);
   TCanvas* canvas = tdrCanvas("VertexReconstructionEfficiency", frame, 13, 11, true);
   reco->SetMarkerStyle(kFullCircle);
   emulation->SetMarkerStyle(kOpenCircle);
   reco->SetMarkerColor(kRed);
   emulation->SetMarkerColor(kRed-2);
   reco->SetLineColor(kRed);
   emulation->SetLineColor(kRed-2);
   reco->Draw("same");
   emulation->Draw("same");
   TLegend* l = tdrLeg(0.4,0.35,0.6,0.15);
   l->AddEntry((TObject*)0, "FastHisto","");
   l->AddEntry((TObject*)0, "#epsilon=0.10cm","");
   l->AddEntry((TObject*)0, "#delta=0.10cm","");
   l->AddEntry(reco, "Simulation", "pl");
   l->AddEntry(emulation, "Emulation", "pl");
   return canvas;
}

//______________________________________________________________________________
// Run with: `root -n -b validation.C+`
void validation() {
   gStyle->SetOptStat(0);
   gEnv->SetValue("TFile.AsyncPrefetching", 1);

   //TFile* fin = TFile::Open("output_numEvent10000.root","READ");
   //TTree* t = (TTree*)fin->Get("L1TVertexNTupler/l1VertexReco");
   TChain* t = new TChain();
   std::string url_string = "root://cmseos.fnal.gov/";
   std::string input = "/store/user/aperloff/2021_04_17-VertexNtuples/Phase2HLTTDRSummer20ReRECOMiniAOD/TT_TuneCP5/";
   t->Add((url_string+input+"14TeV-powheg-pythia8_*.root").c_str());
   setChainElementNames(t,"L1TVertexNTupler/l1VertexReco");
   
   vector<float> *tp_pt = new vector<float>;
   vector<float> *tp_z0 = new vector<float>;
   vector<float> *tracks_pt = new vector<float>;
   vector<float> *tracks_z0 = new vector<float>;
   vector<float> *fh_z0 = new vector<float>;
   vector<float> *fh_sumpt = new vector<float>;
   vector<unsigned int> *fh_ntracks = new vector<unsigned int>;
   vector<float> *fh_em_z0 = new vector<float>;
   vector<float> *fh_em_sumpt = new vector<float>;
   float gen_vtx_z0;

   TBranch *b_tp_pt, *b_tp_z0,
           *b_tracks_pt, *b_tracks_z0,
           *b_fh_z0, *b_fh_sumpt, *b_fh_ntracks,
           *b_fh_em_z0, *b_fh_em_sumpt,
           *b_genvtx_z0;
   t->SetBranchAddress("trueTracks_pt",&tp_pt,&b_tp_pt);
   t->SetBranchAddress("trueTracks_z0",&tp_z0,&b_tp_z0);
   t->SetBranchAddress("recoTracks_hybrid_pt",&tracks_pt,&b_tracks_pt);
   t->SetBranchAddress("recoTracks_hybrid_z0",&tracks_z0,&b_tracks_z0);
   t->SetBranchAddress("recoVertices_FastHisto_z0",&fh_z0,&b_fh_z0);
   t->SetBranchAddress("recoVertices_FastHisto_sumPt",&fh_sumpt,&b_fh_sumpt);
   t->SetBranchAddress("recoVertices_FastHisto_numTracks",&fh_ntracks,&b_fh_ntracks);
   t->SetBranchAddress("emulationVertices_FastHistoEmulation_z0",&fh_em_z0,&b_fh_em_z0);
   t->SetBranchAddress("emulationVertices_FastHistoEmulation_sumPt",&fh_em_sumpt,&b_fh_em_sumpt);
   t->SetBranchAddress("genVertex_z0",&gen_vtx_z0,&b_genvtx_z0);

   TFile* fout = TFile::Open("validation.root","RECREATE");

   TH1F* h_tp_pt = new TH1F("tp_pt","tp_pt;p_{T}^{TP};entries",30,0,150); h_tp_pt->Sumw2();
   TH1F* h_tp_z0 = new TH1F("tp_z0","tp_z0;z0^{TP} [cm];entries;",300,-15,15); h_tp_z0->Sumw2();
   TH1F* h_tracks_pt = new TH1F("tracks_pt","tracks_pt;p_{T}^{tracks};entries",30,0,150); h_tracks_pt->Sumw2();
   TH1F* h_tracks_z0 = new TH1F("tracks_z0","tracks_z0;z0^{tracks} [cm];entries",300,-15,15); h_tracks_z0->Sumw2();
   TH1F* h_fh_z0 = new TH1F("fh_z0","fh_z0;z0^{FH};events",300,-15,15); h_tracks_z0->Sumw2();
   TH1F* h_fh_sumpt = new TH1F("fh_sumpt","fh_sumpt;#Sigma_{tracks #in vtx} p_{T};events",100,0,500); h_fh_sumpt->Sumw2();
   TH1F* h_fh_ntracks = new TH1F("fh_ntracks","fh_ntracks;N_{tracks};events",80,0,80); h_fh_ntracks->Sumw2();
   TH1F* h_reco_gen_z0 = new TH1F("reco_gen_z0","reco_gen_z0;z_{reco}-z_{gen} (cm);Fraction of Events / 0.025 cm",80,-1,1); h_reco_gen_z0->Sumw2();
   TH1F* h_emulation_gen_z0 = new TH1F("emulation_gen_z0","emulation_gen_z0;z_{emu}-z_{gen} (cm);Fraction of Events / 0.025 cm",80,-1,1); h_emulation_gen_z0->Sumw2();
   TH1F* h_vertexMatch = new TH1F("h_vertexMatch","h_vertexMatch;N_{match}(sim, emu);Events",11,0,11); h_vertexMatch->Sumw2();
   TH1F* h_vertexMatch_ordered = new TH1F("h_vertexMatch_ordered","h_vertexMatch_ordered;N_{match}^{ordered}(sim, emu);Events",11,0,11); h_vertexMatch_ordered->Sumw2();
   TProfile* p_reco_res_gen_z0 = new TProfile("reco_res_gen_z0", "reco_res_gen_z0;z_{gen} (cm);|z_{reco}-z_{gen}| (cm)", 60, -15., 15., 0, 0.1);
   TProfile* p_emulation_res_gen_z0 = new TProfile("emulation_res_gen_z0", "emulation_res_gen_z0;z_{gen} (cm);|z_{emu}-z_{gen}| (cm)", 60, -15., 15., 0, 0.1);
   TEfficiency* e_vertexRecoEff_reco_gen_z0 = new TEfficiency("e_vertexRecoEff_reco_gen_z0","e_vertexRecoEff_reco_gen_z0;z_{gen} (cm);Vertex Reconstruction Efficiency",30,-15.,15.);
   TEfficiency* e_vertexRecoEff_emulation_gen_z0 = new TEfficiency("e_vertexRecoEff_emulation_gen_z0","e_vertexRecoEff_emulation_gen_z0;z_{gen} (cm);Vertex Reconstruction Efficiency",30,-15.,15.);


   size_t nEntries = t->GetEntries();
   for (unsigned int ientry = 0; ientry < nEntries; ientry++) {
      if (ientry % 500 == 0) cout << "Doing entry " << ientry << " ... " << endl;
      t->GetEntry(ientry);
      for(unsigned int itp = 0; itp < tp_pt->size(); itp++) {
         h_tp_pt->Fill(tp_pt->at(itp));
         h_tp_z0->Fill(tp_z0->at(itp));
      }
      for(unsigned int itk = 0; itk < tracks_pt->size(); itk++) {
         h_tracks_pt->Fill(tracks_pt->at(itk));
         h_tracks_z0->Fill(tracks_z0->at(itk));
      }
      int nmatch = 0;
      int nmatch_before_fail = 0;
      bool good_match = true;
      for(unsigned int ivtx = 0; ivtx < fh_z0->size(); ivtx++) {
         if (std::abs(fh_z0->at(ivtx)-fh_em_z0->at(ivtx)) < 0.1) {
            nmatch++;
            if (good_match) {
               nmatch_before_fail++;
            }
         }
         else {
            good_match = false;
         }
      }
      h_fh_z0->Fill(fh_z0->at(0));
      h_fh_sumpt->Fill(fh_sumpt->at(0));
      h_fh_ntracks->Fill(fh_ntracks->at(0));
      h_reco_gen_z0->Fill(fh_z0->at(0)-gen_vtx_z0);
      h_emulation_gen_z0->Fill(fh_em_z0->at(0)-gen_vtx_z0);
      h_vertexMatch->Fill(nmatch);
      h_vertexMatch_ordered->Fill(nmatch_before_fail);
      p_reco_res_gen_z0->Fill(gen_vtx_z0, std::abs(fh_z0->at(0)-gen_vtx_z0));
      p_emulation_res_gen_z0->Fill(gen_vtx_z0, std::abs(fh_em_z0->at(0)-gen_vtx_z0));
      e_vertexRecoEff_reco_gen_z0->Fill((std::abs(fh_z0->at(0)-gen_vtx_z0) < 0.1), gen_vtx_z0);
      e_vertexRecoEff_emulation_gen_z0->Fill((std::abs(fh_em_z0->at(0)-gen_vtx_z0) < 0.1), gen_vtx_z0);
   }
   
   h_tp_pt->Write();
   h_tp_z0->Write();
   h_tracks_pt->Write();
   h_tracks_z0->Write();
   h_fh_z0->Write();
   h_fh_sumpt->Write();
   h_fh_ntracks->Write();
   h_reco_gen_z0->Write();
   h_emulation_gen_z0->Write();
   h_vertexMatch->Write();
   h_vertexMatch_ordered->Write();
   p_reco_res_gen_z0->Write();
   p_emulation_res_gen_z0->Write();
   e_vertexRecoEff_reco_gen_z0->Write();
   e_vertexRecoEff_emulation_gen_z0->Write();
   
   TCanvas* cEff = makeEfficiencyCanvas(e_vertexRecoEff_reco_gen_z0, e_vertexRecoEff_emulation_gen_z0);
   cEff->SaveAs("efficiency.png");
   cEff->Write();
   TCanvas* cResolution = makeResolutionCanvas(p_reco_res_gen_z0, p_emulation_res_gen_z0);
   cResolution->SaveAs("resolution.png");
   cResolution->Write();
   TCanvas* cResidual = makeResidualCanvas(h_reco_gen_z0, h_emulation_gen_z0);
   cResidual->SaveAs("residual.png");
   cResidual->Write();
   TCanvas* cNTrack = makeNTrackCanvas(h_fh_ntracks, nullptr);
   cNTrack->SaveAs("ntracks.png");
   cNTrack->Write();
   TCanvas* cSumPt = makeSumPtCanvas(h_fh_sumpt, nullptr);
   cSumPt->SaveAs("sumpt.png");
   cSumPt->Write();
   TCanvas* cMatch = makeMatchCanvas(h_vertexMatch, h_vertexMatch_ordered);
   cMatch->SaveAs("n_matched_vertices.png");
   cMatch->Write();

   fout->Close();   
}
