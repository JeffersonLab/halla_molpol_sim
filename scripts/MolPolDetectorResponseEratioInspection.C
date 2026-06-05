///////////////////////////////////////////////////////////////
// MolPolDetectorResponseEratioInspection.C
// Standalone investigation of the E_detected / E_actual ratio
// distribution bumps in the detector response output.
//
// Reads a *_detresponse.root file and analyzes events by their
// energy ratio, dumping detailed per-hit and per-throw info
// for events in configurable ratio windows.
//
// Usage:
//   root -l 'MolPolDetectorResponseEratioInspection.C("input_detresponse.root")'
//   root -l 'MolPolDetectorResponseEratioInspection.C("input_detresponse.root", 0.15, 0.35)'
//   root -l 'MolPolDetectorResponseEratioInspection.C("input_detresponse.root", 0.45, 0.55, 50)'
//
///////////////////////////////////////////////////////////////

#include "MolPolDetectorResponse.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"

#include <cstdio>
#include <cmath>

void MolPolDetectorResponseEratioInspection(const char *filename) {

    TFile *f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        printf("Error: cannot open %s\n", filename);
        return;
    }

    TTree *tree = (TTree *)f->Get("TDetResp");
    if (!tree) {
        printf("Error: TDetResp tree not found in %s\n", filename);
        f->Close();
        return;
    }

    SetupDetRespBranches(tree);

    Long64_t nEntries = tree->GetEntries();
    printf("=== Energy Ratio Investigation ===\n");
    printf("File: %s\n", filename);
    printf("Entries: %lld\n\n", nEntries);

    // --- Histograms: ratio split by hit multiplicity ---
    TH1F *hRatioAll  = new TH1F("hRatioAll",  "All events;E_{det}/E_{actual};Events",    300, 0.0, 1.5);
    TH1F *hRatio1hit = new TH1F("hRatio1hit", "1 hit;E_{det}/E_{actual};Events",         300, 0.0, 1.5);
    TH1F *hRatio2hit = new TH1F("hRatio2hit", "2 hits;E_{det}/E_{actual};Events",        300, 0.0, 1.5);
    TH1F *hRatio3plus = new TH1F("hRatio3plus","3+ hits;E_{det}/E_{actual};Events",      300, 0.0, 1.5);

    hRatioAll->SetLineColor(kBlack);    hRatioAll->SetLineWidth(2);
    hRatio1hit->SetLineColor(kRed);     hRatio1hit->SetLineWidth(2);
    hRatio2hit->SetLineColor(kBlue);    hRatio2hit->SetLineWidth(2);
    hRatio3plus->SetLineColor(kGreen+2);hRatio3plus->SetLineWidth(2);

    // 2D: ratio vs hitZ (throw 0) — maps out the incomplete gamma function
    TH2F *hRatioVsHitZ = new TH2F("hRatioVsHitZ",
        "E_{det}/E_{actual} vs Shower Start Depth (per-hit);hitZ [cm];E_{det}/E_{actual}",
        300, 0.0, 30.0, 300, 0.0, 1.5);

    // 2D: ratio vs hitE — shows threshold discontinuity (single-hit only)
    TH2F *hRatioVsHitE = new TH2F("hRatioVsHitE",
        "E_{det}/E_{actual} vs Hit Energy (per-hit);E_{hit} [GeV];E_{det}/E_{actual}",
        300, 0.0, 6.0, 300, 0.0, 1.5);

    // 2D: ratio vs sumHitE — per-event capture (all events)
    TH2F *hRatioVsSumE = new TH2F("hRatioVsSumE",
        "E_{det}/E_{actual} vs Event Energy (per-event);E_{actual} [GeV];E_{det}/E_{actual}",
        300, 0.0, 12.0, 300, 0.0, 1.5);

    Int_t nDumped = 0;
    Long64_t nInWindow = 0;
    Long64_t nWithHits = 0;

    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        tree->GetEntry(iEntry);

        if (dr_nHitsDet9 == 0) continue;
        nWithHits++;

        // Compute actual and detected energies (throw 0)
        Double_t sumHitE = 0.0;
        for (Int_t h = 0; h < dr_nHitsDet9; h++)
            sumHitE += dr_hitE[h];

        if (sumHitE <= 0.0) continue;

        Double_t detTotal = 0.0;
        for (Int_t s = 0; s < 8; s++)
            detTotal += dr_pmtE[s];

        Double_t ratio = detTotal / sumHitE;

        // Fill multiplicity-separated histograms
        hRatioAll->Fill(ratio);
        if (dr_nHitsDet9 == 1)      hRatio1hit->Fill(ratio);
        else if (dr_nHitsDet9 == 2) hRatio2hit->Fill(ratio);
        else                        hRatio3plus->Fill(ratio);

        // Ratio vs start depth (throw 0)
        hRatioVsHitZ->Fill(dr_hitZ[0], ratio);

        // Ratio vs hit energy (single-hit events only for clean picture)
        if (dr_nHitsDet9 == 1)
            hRatioVsHitE->Fill(dr_hitE[0], ratio);

        // Ratio vs total event energy (all events)
        hRatioVsSumE->Fill(sumHitE, ratio);
    }

    printf("Events with hits: %lld\n", nWithHits);

    // --- Plot ---
    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);

    TCanvas *cInspect = new TCanvas("cInspect",
        "Energy Ratio Investigation", 1200, 1000);
    cInspect->Divide(2, 2);

    // Panel 1: ratio by multiplicity
    cInspect->cd(1);
    gPad->SetLogy();
    hRatioAll->SetMinimum(0.9);
    hRatioAll->SetTitle("E_{det}/E_{actual} by Hit Multiplicity");
    hRatioAll->Draw("HISTE");
    hRatio1hit->Draw("HISTE same");
    hRatio2hit->Draw("HISTE same");
    hRatio3plus->Draw("HISTE same");

    // Shade the dump window
    TLine *lLow = new TLine(ratioLow, hRatioAll->GetMinimum(),
                              ratioLow, hRatioAll->GetMaximum());
    lLow->SetLineColor(kMagenta);
    lLow->SetLineStyle(2);
    lLow->SetLineWidth(2);
    lLow->Draw("same");
    TLine *lHigh = new TLine(ratioHigh, hRatioAll->GetMinimum(),
                               ratioHigh, hRatioAll->GetMaximum());
    lHigh->SetLineColor(kMagenta);
    lHigh->SetLineStyle(2);
    lHigh->SetLineWidth(2);
    lHigh->Draw("same");

    TLegend *leg1 = new TLegend(0.15, 0.65, 0.45, 0.88);
    leg1->AddEntry(hRatioAll, "All", "l");
    leg1->AddEntry(hRatio1hit, "1 hit", "l");
    leg1->AddEntry(hRatio2hit, "2 hits", "l");
    leg1->AddEntry(hRatio3plus, "3+ hits", "l");
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->Draw();

    // Panel 2: Ratio vs hit energy — threshold + lateral leakage (per-hit)
    cInspect->cd(2);
    gPad->SetRightMargin(0.15);
    gPad->SetLogx();
    gPad->SetLogz();
    hRatioVsHitE->SetStats(kFALSE);
    hRatioVsHitE->Draw("COLZ");
    TLine *lThresh = new TLine(MIN_SHOWER_E_GEV, 0.0,
                                MIN_SHOWER_E_GEV, 1.2);
    lThresh->SetLineColor(kRed);
    lThresh->SetLineWidth(2);
    lThresh->SetLineStyle(2);
    lThresh->Draw("same");

    // Panel 3: Ratio vs hitZ — incomplete gamma mapping (per-hit)
    cInspect->cd(3);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    hRatioVsHitZ->SetStats(kFALSE);
    hRatioVsHitZ->Draw("COLZ");

    // Panel 4: Ratio vs total event energy (per-event)
    cInspect->cd(4);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    hRatioVsSumE->SetStats(kFALSE);
    hRatioVsSumE->Draw("COLZ");

    cInspect->Update();
}