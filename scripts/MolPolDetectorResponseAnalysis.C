///////////////////////////////////////////////////////////////
// MolPolDetectorResponseAnalysis.C
// Combined analysis and plotting for the MolPol detector
// response model output. Computes analyzing power with
// configurable active PMT rows and energy threshold, and
// produces PMT energy spectra, L/R correlations, and A_zz
// histograms.
//
// Also provides standalone interactive functions:
//   PlotEventShower(tShower, eventNum)
//   PlotShowerParams(tShower, eventNum)
//
// Usage:
//   root 'MolPolDetectorResponseAnalysis.C("input_detresponse.root")'
//   root 'MolPolDetectorResponseAnalysis.C("input_detresponse.root", 0.5)'
//   root 'MolPolDetectorResponseAnalysis.C("input_detresponse.root", 0.5, "23")'
//
// Interactive shower visualization (requires writeShowerDiag):
//   root [0] .L MolPolDetectorResponseAnalysis.C
//   root [1] TChain *t = new TChain("TShower");
//   root [2] t->Add("input_detresponse.root");
//   root [3] PlotEventShower(t, 1);    // 3D shower for event 1
//   root [4] PlotShowerParams(t, 1);   // RC, RT, p, dE/dt for event 1
//
///////////////////////////////////////////////////////////////

#ifndef MOLPOL_DETECTOR_RESPONSE_ANALYSIS_C
#define MOLPOL_DETECTOR_RESPONSE_ANALYSIS_C

#include "MolPolDetectorResponseAnalysis.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TView.h"
#include "TAxis.h"
#include "TColor.h"
#include "TPad.h"
#include "TGaxis.h"

#include <cstdio>
#include <vector>
#include <cmath>

// PlotEventShower()
// 3D visualization of the modeled shower for a specific event
// in the TShower diagnostic tree.
//
// eventNum is 1-based: PlotEventShower(tShower, 1) plots the
// first unique event, PlotEventShower(tShower, 10) the tenth.
//
// Each depth slice is drawn as a single TPolyLine3D ring at
// the 90% energy containment radius. Ring color indicates
// deposited energy (blue=low, red=high). Line width scales
// with energy fraction.
//
// Axes: X = depth into calorimeter [cm] (beam enters from left)
//       Y = horizontal position on face [cm] (hitLx)
//       Z = vertical position on face [cm] (hitLy)
///////////////////////////////////////////////////////////////

void PlotEventShower(TTree *tShower, Int_t eventNum = 1) {

    if (!tShower) {
        printf("PlotEventShower: TShower tree is null. "
               "Was writeShowerDiag enabled?\n");
        return;
    }

    Long64_t nEntries = tShower->GetEntries();
    if (nEntries == 0) {
        printf("PlotEventShower: TShower tree is empty.\n");
        return;
    }

    // Wire up branches for reading
    Int_t   sh_eventIdx, sh_hitIdx, sh_nDepthSteps;
    Float_t sh_hitX, sh_hitY, sh_hitE;
    Float_t sh_depthT[MAX_DEPTH_STEPS];
    Float_t sh_sliceE[MAX_DEPTH_STEPS];
    Float_t sh_RC[MAX_DEPTH_STEPS];
    Float_t sh_RT[MAX_DEPTH_STEPS];
    Float_t sh_p[MAX_DEPTH_STEPS];

    tShower->SetBranchAddress("eventIdx",    &sh_eventIdx);
    tShower->SetBranchAddress("hitIdx",      &sh_hitIdx);
    tShower->SetBranchAddress("hitX",        &sh_hitX);
    tShower->SetBranchAddress("hitY",        &sh_hitY);
    tShower->SetBranchAddress("hitE",        &sh_hitE);
    tShower->SetBranchAddress("nDepthSteps", &sh_nDepthSteps);
    tShower->SetBranchAddress("depthT",      sh_depthT);
    tShower->SetBranchAddress("sliceE",      sh_sliceE);
    tShower->SetBranchAddress("RC",          sh_RC);
    tShower->SetBranchAddress("RT",          sh_RT);
    tShower->SetBranchAddress("p",           sh_p);

    // --- Build ordered list of unique eventIdx values ---
    std::vector<Int_t> uniqueList;
    for (Long64_t i = 0; i < nEntries; i++) {
        tShower->GetEntry(i);
        Bool_t found = false;
        for (size_t u = 0; u < uniqueList.size(); u++) {
            if (uniqueList[u] == sh_eventIdx) { found = true; break; }
        }
        if (!found) uniqueList.push_back(sh_eventIdx);
    }

    Int_t nUniqueEvents = (Int_t)uniqueList.size();
    printf("PlotEventShower: %d unique events in TShower\n", nUniqueEvents);

    if (eventNum < 1 || eventNum > nUniqueEvents) {
        printf("PlotEventShower: eventNum %d out of range [1, %d]\n",
               eventNum, nUniqueEvents);
        return;
    }

    Int_t targetEventIdx = uniqueList[eventNum - 1];
    printf("PlotEventShower: plotting event #%d (eventIdx = %d)\n",
           eventNum, targetEventIdx);

    // --- Collect all hits for this event ---
    std::vector<HitSlice> slices;

    Float_t maxEnergy = 0.0;
    Float_t maxDepth_cm = 0.0;
    Int_t   nHitsEvent = 0;
    Int_t   lastHitIdx = -1;
    Float_t totalHitE = 0.0;

    // Store hit entry points for drawing markers
    std::vector<Float_t> hitEntryX, hitEntryY, hitEntryE;
    std::vector<Int_t>   hitEntryIdx;

    for (Long64_t i = 0; i < nEntries; i++) {
        tShower->GetEntry(i);
        if (sh_eventIdx != targetEventIdx) continue;

        if (sh_hitIdx != lastHitIdx) {
            nHitsEvent++;
            totalHitE += sh_hitE;
            lastHitIdx = sh_hitIdx;
            hitEntryX.push_back(sh_hitX);
            hitEntryY.push_back(sh_hitY);
            hitEntryE.push_back(sh_hitE);
            hitEntryIdx.push_back(sh_hitIdx);
            printf("  hit[%d]: E=%.4f GeV, x=%.2f cm, y=%.2f cm, %d depth steps\n",
                   sh_hitIdx, sh_hitE, sh_hitX, sh_hitY, sh_nDepthSteps);
        }

        for (Int_t s = 0; s < sh_nDepthSteps; s++) {
            if (sh_sliceE[s] <= 0.0f) continue;

            Double_t r90_rm = Solve90PctRadius(sh_p[s], sh_RC[s], sh_RT[s]);
            Double_t r90_cm = r90_rm * PbConst::RM;
            Double_t depth_cm = sh_depthT[s] * PbConst::X0;

            HitSlice hs;
            hs.x        = sh_hitX;
            hs.y        = sh_hitY;
            hs.depth_cm = (Float_t)depth_cm;
            hs.energy   = sh_sliceE[s];
            hs.r90      = (Float_t)r90_cm;
            hs.rc       = sh_RC[s];
            hs.rt       = sh_RT[s];
            hs.p        = sh_p[s];
            hs.hitIdx   = sh_hitIdx;
            slices.push_back(hs);

            if (sh_sliceE[s] > maxEnergy) maxEnergy = sh_sliceE[s];
            if (depth_cm > maxDepth_cm)   maxDepth_cm = (Float_t)depth_cm;
        }
    }

    printf("  %d hit(s), %.3f GeV total, %d depth slices\n",
           nHitsEvent, totalHitE, (Int_t)slices.size());

    if (slices.empty()) {
        printf("PlotEventShower: no shower slices found for this event.\n");
        return;
    }

    // --- Debug: print first 20 slices per hit for large-ring diagnosis ---
    printf("  --- First 20 slices per hit (for large-ring diagnosis) ---\n");
    Int_t prevHit = -1;
    Int_t sliceCount = 0;
    for (size_t si = 0; si < slices.size(); si++) {
        HitSlice &hs = slices[si];
        if (hs.hitIdx != prevHit) {
            prevHit = hs.hitIdx;
            sliceCount = 0;
            printf("  hit[%d]:\n", hs.hitIdx);
        }
        if (sliceCount < 20) {
            printf("    slice %d: depth=%.2f cm [%.1f X0], E=%.4f GeV, "
                   "r90=%.2f cm, RC=%.4f RT=%.4f p=%.4f\n",
                   sliceCount, hs.depth_cm, hs.depth_cm / PbConst::X0,
                   hs.energy,
                   hs.r90, hs.rc, hs.rt, hs.p);
        }
        sliceCount++;
    }


    // --- Create canvas with space for colorbar ---
    TString cName  = TString::Format("cShower_ev%d", eventNum);
    TString cTitle = TString::Format("Event #%d Shower", eventNum);

    TCanvas *c = new TCanvas(cName.Data(), cTitle.Data(), 1100, 700);

    // Set palette 55 (kRainBow): violet->blue->cyan->green->yellow->orange->red
    gStyle->SetPalette(55);

    // Main pad for the 3D shower (leave room for colorbar)
    TPad *padMain = new TPad("padMain", "", 0.0, 0.0, 0.88, 1.0);
    padMain->Draw();
    padMain->cd();

    // --- Establish 3D axes with a dummy TGraph2D ---
    // X = depth [cm] (0 to caloDepth)
    // Y = hitX [cm] (-9 to +9, detector half-width)
    // Z = hitY [cm] (-15 to +15, detector half-height)
    TGraph2D *gAxis = new TGraph2D(2);
    gAxis->SetPoint(0, 0.0,  -9.0, -15.0);
    gAxis->SetPoint(1, DetGeom::caloDepth, 9.0, 15.0);
    gAxis->SetTitle(TString::Format(
        "Event #%d Shower;Depth [cm];X [cm];Y [cm]", eventNum).Data());
    gAxis->SetMarkerStyle(1);
    gAxis->SetMarkerSize(0);
    gAxis->Draw("p0");
    padMain->Update();

    // --- Draw shower discs as 3 TPolyLine3D rings per slice ---
    // Three rings span the depth of each 1-X0 slice, giving a
    // disc-like appearance. Ring at the slice midpoint +/- X0/3.
    const Int_t nPtsRing = 48;
    const Int_t nRingsPerSlice = 3;
    Double_t dz = PbConst::X0 / 3.0;  // ring spacing in cm

    for (size_t si = 0; si < slices.size(); si++) {
        HitSlice &hs = slices[si];

        Double_t eFrac = hs.energy / maxEnergy;
        if (eFrac < 0.01) continue;

        // Color from palette (matches colorbar)
        Int_t nColors = gStyle->GetNumberOfColors();
        Int_t palIdx = (Int_t)(eFrac * (nColors - 1));
        if (palIdx >= nColors) palIdx = nColors - 1;
        Int_t color = gStyle->GetColorPalette(palIdx);

        // Line width scales with energy: 1 (low) to 3 (high)
        Int_t lineW = 1 + (Int_t)(2.0 * eFrac);

        // Draw 3 rings at different depths within the slice,
        // clipped to the physical block boundaries.
        Double_t clipY = DetGeom::blockHalfX;  // ±9 cm in Y (hitX dir)
        Double_t clipZ = DetGeom::blockHalfY;  // ±15 cm in Z (hitY dir)

        for (Int_t ir = -1; ir <= 1; ir++) {
            Double_t ringDepth = hs.depth_cm + ir * dz;
            if (ringDepth < 0.0) continue;

            // Generate ring points and check bounds
            Double_t yPts[nPtsRing + 1];
            Double_t zPts[nPtsRing + 1];
            Bool_t   inside[nPtsRing + 1];

            for (Int_t ip = 0; ip <= nPtsRing; ip++) {
                Double_t theta = 2.0 * TMath::Pi() * ip / nPtsRing;
                yPts[ip] = hs.x + hs.r90 * cos(theta);
                zPts[ip] = hs.y + hs.r90 * sin(theta);
                inside[ip] = (yPts[ip] >= -clipY && yPts[ip] <= clipY &&
                              zPts[ip] >= -clipZ && zPts[ip] <= clipZ);
            }

            // Draw runs of consecutive inside points as separate polylines
            Int_t runStart = -1;
            for (Int_t ip = 0; ip <= nPtsRing; ip++) {
                if (inside[ip]) {
                    if (runStart < 0) runStart = ip;
                } else {
                    // End of an inside run
                    if (runStart >= 0) {
                        Int_t nRunPts = ip - runStart;
                        TPolyLine3D *seg = new TPolyLine3D(nRunPts);
                        for (Int_t k = 0; k < nRunPts; k++)
                            seg->SetPoint(k, ringDepth,
                                yPts[runStart + k], zPts[runStart + k]);
                        seg->SetLineColor(color);
                        seg->SetLineWidth(lineW);
                        seg->Draw("same");
                        runStart = -1;
                    }
                }
            }
            // Flush final run if ring ends inside
            if (runStart >= 0) {
                Int_t nRunPts = nPtsRing + 1 - runStart;
                TPolyLine3D *seg = new TPolyLine3D(nRunPts);
                for (Int_t k = 0; k < nRunPts; k++)
                    seg->SetPoint(k, ringDepth,
                        yPts[runStart + k], zPts[runStart + k]);
                seg->SetLineColor(color);
                seg->SetLineWidth(lineW);
                seg->Draw("same");
            }
        }
    }

    // --- Draw hit entry points on the detector face (depth = 0) ---
    Int_t hitMarkerColors[6] = { kRed, kBlue, kGreen+2, kMagenta, kCyan+1, kOrange+1 };
    for (size_t ih = 0; ih < hitEntryX.size(); ih++) {
        TGraph2D *gHit = new TGraph2D(1);
        gHit->SetPoint(0, 0.0, hitEntryX[ih], hitEntryY[ih]);
        gHit->SetMarkerStyle(29);  // star
        gHit->SetMarkerSize(2.0);
        gHit->SetMarkerColor(hitMarkerColors[ih % 6]);
        gHit->Draw("p0 same");
    }

    // --- Draw PMT segment grid at depth = 0 (detector face) ---
    Double_t halfW = 9.0;
    for (Int_t iy = 0; iy < 5; iy++) {
        TPolyLine3D *line = new TPolyLine3D(2);
        line->SetPoint(0, 0.0, -halfW, DetGeom::rowYBound[iy]);
        line->SetPoint(1, 0.0,  halfW, DetGeom::rowYBound[iy]);
        line->SetLineColor(kGray + 2);
        line->SetLineWidth(1);
        line->Draw("same");
    }
    {
        TPolyLine3D *line = new TPolyLine3D(2);
        line->SetPoint(0, 0.0, 0.0, -DetGeom::blockHalfY);
        line->SetPoint(1, 0.0, 0.0,  DetGeom::blockHalfY);
        line->SetLineColor(kGray + 2);
        line->SetLineWidth(1);
        line->Draw("same");
    }
    {
        TPolyLine3D *box = new TPolyLine3D(5);
        box->SetPoint(0, 0.0, -halfW, -DetGeom::blockHalfY);
        box->SetPoint(1, 0.0,  halfW, -DetGeom::blockHalfY);
        box->SetPoint(2, 0.0,  halfW,  DetGeom::blockHalfY);
        box->SetPoint(3, 0.0, -halfW,  DetGeom::blockHalfY);
        box->SetPoint(4, 0.0, -halfW, -DetGeom::blockHalfY);
        box->SetLineColor(kBlack);
        box->SetLineWidth(2);
        box->Draw("same");
    }

    padMain->Update();

    // --- Colorbar on the right side ---
    c->cd();
    TPad *padCbar = new TPad("padCbar", "", 0.88, 0.10, 1.0, 0.90);
    padCbar->SetLeftMargin(0.02);
    padCbar->SetRightMargin(0.55);
    padCbar->SetTopMargin(0.02);
    padCbar->SetBottomMargin(0.02);
    padCbar->Draw();
    padCbar->cd();

    gStyle->SetOptStat(0);

    TH2F *hCbar = new TH2F("hCbar", "",
        1, 0.0, 1.0,
        100, 0.0, maxEnergy);
    for (Int_t ib = 1; ib <= 100; ib++)
        hCbar->SetBinContent(1, ib, hCbar->GetYaxis()->GetBinCenter(ib));

    // Hide the default axes (we draw our own on the right)
    hCbar->GetXaxis()->SetLabelSize(0);
    hCbar->GetXaxis()->SetTickLength(0);
    hCbar->GetYaxis()->SetLabelSize(0);
    hCbar->GetYaxis()->SetTickLength(0);
    hCbar->GetYaxis()->SetTitle("");
    hCbar->Draw("COL");

    padCbar->Update();

    // Draw axis with labels on the right side
    TGaxis *rightAxis = new TGaxis(
        1.0, 0.0, 1.0, maxEnergy,   // x1,y1 -> x2,y2 (user coords)
        0.0, maxEnergy,               // axis value range
        510, "+L");                    // 510 = ~10 primary divisions, +L = labels on right
    rightAxis->SetTitle("Slice Energy Deposited [GeV]");
    rightAxis->SetTitleOffset(1.8);
    rightAxis->SetTitleSize(0.14);
    rightAxis->SetTitleFont(42);
    rightAxis->SetLabelSize(0.13);
    rightAxis->SetLabelFont(42);
    rightAxis->Draw();

    padCbar->Modified();
    padCbar->Update();

    gStyle->SetOptStat(1111);

    c->Update();
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// PlotShowerParams()
// Four-panel plot of the Grindhammer-Peters shower parameters
// for a specific event:
//   Panel 1: RC vs tau (core radius)
//   Panel 2: RT vs tau (tail radius)
//   Panel 3: p vs tau (core weight)
//   Panel 4: <1/E dE/dr> vs r (integrated radial profile)
//
// Panels 1-3 reproduce the style of Figs. 9-11 in
// arXiv:hep-ex/0001020. Panel 4 reproduces Fig. 13 (upper).
// All hits overlaid, distinguished by color.
// Uses the same 1-based event indexing as PlotEventShower().
///////////////////////////////////////////////////////////////

void PlotShowerParams(TTree *tShower, Int_t eventNum = 1) {

    if (!tShower) {
        printf("PlotShowerParams: TShower tree is null.\n");
        return;
    }

    Long64_t nEntries = tShower->GetEntries();
    if (nEntries == 0) {
        printf("PlotShowerParams: TShower tree is empty.\n");
        return;
    }

    // Wire up branches
    Int_t   sh_eventIdx, sh_hitIdx, sh_nDepthSteps;
    Float_t sh_hitX, sh_hitY, sh_hitE;
    Float_t sh_depthT[MAX_DEPTH_STEPS];
    Float_t sh_sliceE[MAX_DEPTH_STEPS];
    Float_t sh_RC[MAX_DEPTH_STEPS];
    Float_t sh_RT[MAX_DEPTH_STEPS];
    Float_t sh_p[MAX_DEPTH_STEPS];

    tShower->SetBranchAddress("eventIdx",    &sh_eventIdx);
    tShower->SetBranchAddress("hitIdx",      &sh_hitIdx);
    tShower->SetBranchAddress("hitX",        &sh_hitX);
    tShower->SetBranchAddress("hitY",        &sh_hitY);
    tShower->SetBranchAddress("hitE",        &sh_hitE);
    tShower->SetBranchAddress("nDepthSteps", &sh_nDepthSteps);
    tShower->SetBranchAddress("depthT",      sh_depthT);
    tShower->SetBranchAddress("sliceE",      sh_sliceE);
    tShower->SetBranchAddress("RC",          sh_RC);
    tShower->SetBranchAddress("RT",          sh_RT);
    tShower->SetBranchAddress("p",           sh_p);

    // --- Build unique event list ---
    std::vector<Int_t> uniqueList;
    for (Long64_t i = 0; i < nEntries; i++) {
        tShower->GetEntry(i);
        Bool_t found = false;
        for (size_t u = 0; u < uniqueList.size(); u++) {
            if (uniqueList[u] == sh_eventIdx) { found = true; break; }
        }
        if (!found) uniqueList.push_back(sh_eventIdx);
    }

    Int_t nUniqueEvents = (Int_t)uniqueList.size();
    if (eventNum < 1 || eventNum > nUniqueEvents) {
        printf("PlotShowerParams: eventNum %d out of range [1, %d]\n",
               eventNum, nUniqueEvents);
        return;
    }

    Int_t targetEventIdx = uniqueList[eventNum - 1];
    printf("PlotShowerParams: event #%d (eventIdx = %d)\n",
           eventNum, targetEventIdx);

    // --- Collect per-hit data ---
    Int_t hitColors[6] = { kRed, kBlue, kGreen+2, kMagenta, kCyan+1, kOrange+1 };

    std::vector<Int_t> hitIdxList;
    std::vector<Float_t> hitEList;
    std::vector<Double_t> hitTList;  // shower max T in X0
    std::vector< std::vector<Double_t> > tauVecs;
    std::vector< std::vector<Double_t> > depthVecs;  // depth in X0
    std::vector< std::vector<Double_t> > rcVecs;
    std::vector< std::vector<Double_t> > rtVecs;
    std::vector< std::vector<Double_t> > pVecs;
    std::vector< std::vector<Double_t> > sliceWtVecs;  // sliceE / hitE

    Int_t lastHitIdx = -1;
    Int_t currentHitSlot = -1;

    for (Long64_t i = 0; i < nEntries; i++) {
        tShower->GetEntry(i);
        if (sh_eventIdx != targetEventIdx) continue;

        if (sh_hitIdx != lastHitIdx) {
            lastHitIdx = sh_hitIdx;
            currentHitSlot++;
            hitIdxList.push_back(sh_hitIdx);
            hitEList.push_back(sh_hitE);
            tauVecs.push_back(std::vector<Double_t>());
            depthVecs.push_back(std::vector<Double_t>());
            rcVecs.push_back(std::vector<Double_t>());
            rtVecs.push_back(std::vector<Double_t>());
            pVecs.push_back(std::vector<Double_t>());
            sliceWtVecs.push_back(std::vector<Double_t>());
        }

        Double_t hitEMeV = sh_hitE * 1000.0;
        Double_t y = hitEMeV / PbConst::Ec;
        Double_t T = GPShower::Thom(TMath::Log(y));

        // Store T on first encounter of this hit
        if ((Int_t)hitTList.size() <= currentHitSlot)
            hitTList.push_back(T);

        for (Int_t s = 0; s < sh_nDepthSteps; s++) {
            // Require positive, finite values in all fields.
            // NaN fails all comparisons, so use positive-range checks
            // (NaN < X is false, so !(NaN < X) rejects it).
            if (!(sh_depthT[s] > 0.0f && sh_depthT[s] < 100.0f)) continue;
            if (!(sh_sliceE[s] > 0.0f)) continue;
            if (!(sh_RC[s] >= 0.0f && sh_RC[s] < 10.0f)) continue;
            if (!(sh_RT[s] >= 0.0f && sh_RT[s] < 10.0f)) continue;
            if (!(sh_p[s] >= 0.0f && sh_p[s] <= 1.0f)) continue;
            Double_t tau = sh_depthT[s] / T;
            tauVecs[currentHitSlot].push_back(tau);
            depthVecs[currentHitSlot].push_back((Double_t)sh_depthT[s]);
            rcVecs[currentHitSlot].push_back((Double_t)sh_RC[s]);
            rtVecs[currentHitSlot].push_back((Double_t)sh_RT[s]);
            pVecs[currentHitSlot].push_back((Double_t)sh_p[s]);
            Double_t wt = (sh_hitE > 0.0f) ? sh_sliceE[s] / sh_hitE : 0.0;
            sliceWtVecs[currentHitSlot].push_back(wt);
        }
    }

    Int_t nHits = (Int_t)hitIdxList.size();
    printf("  %d hit(s)\n", nHits);

    if (nHits == 0) {
        printf("PlotShowerParams: no hits found.\n");
        return;
    }

    // --- Find axis ranges ---
    Double_t tauMax  = 0.0;
    Double_t rcMax   = 0.0;
    Double_t rtMax   = 0.0;
    Double_t depthMax = 0.0;
    for (Int_t ih = 0; ih < nHits; ih++) {
        for (size_t s = 0; s < tauVecs[ih].size(); s++) {
            if (tauVecs[ih][s]  > tauMax)  tauMax  = tauVecs[ih][s];
            if (rcVecs[ih][s]   > rcMax)   rcMax   = rcVecs[ih][s];
            if (rtVecs[ih][s]   > rtMax)   rtMax   = rtVecs[ih][s];
            if (depthVecs[ih][s] > depthMax) depthMax = depthVecs[ih][s];
        }
    }

    // Print hit info
    for (Int_t ih = 0; ih < nHits; ih++) {
        Double_t hitEMeV = hitEList[ih] * 1000.0;
        Double_t y = hitEMeV / PbConst::Ec;
        Double_t T = GPShower::Thom(TMath::Log(y));
        printf("  hit[%d]: E=%.4f GeV, T=%.2f X0, %d points\n",
               hitIdxList[ih], hitEList[ih], T, (Int_t)tauVecs[ih].size());
    }

    // --- Create 2x2 canvas ---
    TString cName = TString::Format("cShowerParams_ev%d", eventNum);
    TCanvas *cParams = new TCanvas(cName.Data(),
        TString::Format("Event #%d Shower Parameters", eventNum).Data(),
        900, 700);
    cParams->Divide(2, 2, 0.01, 0.01);

    // --- Panel 1: RC vs tau ---
    cParams->cd(1);
    TH1F *hFrameRC = gPad->DrawFrame(0.0, 0.0,
        tauMax * 1.1, rcMax * 1.3);
    hFrameRC->SetTitle("R_{C} vs #tau;#tau = t/T;R_{C} [R_{M}]");

    for (Int_t ih = 0; ih < nHits; ih++) {
        Int_t nPts = (Int_t)tauVecs[ih].size();
        if (nPts == 0) continue;
        TGraph *g = new TGraph(nPts,
            &tauVecs[ih][0], &rcVecs[ih][0]);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.6);
        g->SetMarkerColor(hitColors[ih % 6]);
        g->SetLineColor(hitColors[ih % 6]);
        g->Draw("PL same");
    }

    // --- Panel 2: RT vs tau ---
    cParams->cd(2);
    Double_t rtDisplayMax = std::min(rtMax * 1.3, 5.0);
    TH1F *hFrameRT = gPad->DrawFrame(0.0, 0.0,
        tauMax * 1.1, rtDisplayMax);
    hFrameRT->SetTitle("R_{T} vs #tau;#tau = t/T;R_{T} [R_{M}]");

    for (Int_t ih = 0; ih < nHits; ih++) {
        Int_t nPts = (Int_t)tauVecs[ih].size();
        if (nPts == 0) continue;
        TGraph *g = new TGraph(nPts,
            &tauVecs[ih][0], &rtVecs[ih][0]);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.6);
        g->SetMarkerColor(hitColors[ih % 6]);
        g->SetLineColor(hitColors[ih % 6]);
        g->Draw("PL same");
    }

    // --- Panel 3: p vs tau ---
    cParams->cd(3);
    TH1F *hFrameP = gPad->DrawFrame(0.0, 0.0,
        tauMax * 1.1, 1.1);
    hFrameP->SetTitle("p vs #tau;#tau = t/T;p (core weight)");

    for (Int_t ih = 0; ih < nHits; ih++) {
        Int_t nPts = (Int_t)tauVecs[ih].size();
        if (nPts == 0) continue;
        TGraph *g = new TGraph(nPts,
            &tauVecs[ih][0], &pVecs[ih][0]);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.6);
        g->SetMarkerColor(hitColors[ih % 6]);
        g->SetLineColor(hitColors[ih % 6]);
        g->Draw("PL same");
    }

    // --- Panel 4: Longitudinal profile 1/E dE/dt vs t ---
    // Energy deposition per radiation length, matching Fig. 4
    // of arXiv:hep-ex/0001020.
    //
    // Data points: sliceWt[s] = sliceE/hitE is the fraction per
    // 1 X0 slice, so 1/E dE/dt ≈ sliceWt (since dt = 1 X0).
    // Analytic curve: gamma distribution f(t) from Eq. 2.
    cParams->cd(4);

    // Find y-range from data
    Double_t longMax = 0.0;
    for (Int_t ih = 0; ih < nHits; ih++) {
        for (size_t s = 0; s < sliceWtVecs[ih].size(); s++) {
            if (sliceWtVecs[ih][s] > longMax) longMax = sliceWtVecs[ih][s];
        }
    }

    TH1F *hFrameLong = gPad->DrawFrame(0.0, 0.0,
        depthMax * 1.1, longMax * 1.3);
    hFrameLong->SetTitle(
        "1/E  dE/dt vs t;t [X_{0}];"
        "1/E  dE/dt [X_{0}^{-1}]");

    for (Int_t ih = 0; ih < nHits; ih++) {
        Int_t nPts = (Int_t)depthVecs[ih].size();
        if (nPts == 0) continue;

        // Data points (discrete slice weights)
        TGraph *gData = new TGraph(nPts,
            &depthVecs[ih][0], &sliceWtVecs[ih][0]);
        gData->SetMarkerStyle(20);
        gData->SetMarkerSize(0.6);
        gData->SetMarkerColor(hitColors[ih % 6]);
        gData->SetLineColor(hitColors[ih % 6]);
        gData->Draw("P same");

        // Analytic gamma distribution curve (Eq. 2)
        Double_t T = hitTList[ih];
        Double_t hitEMeV = hitEList[ih] * 1000.0;
        Double_t y = hitEMeV / PbConst::Ec;
        Double_t alpha = GPShower::AlphaHom(TMath::Log(y), PbConst::Z);
        if (alpha <= 1.0) alpha = 1.01;
        Double_t beta = (alpha - 1.0) / T;

        const Int_t nCurvePts = 200;
        Double_t tCurve[nCurvePts], fCurve[nCurvePts];
        for (Int_t it = 0; it < nCurvePts; it++) {
            tCurve[it] = (it + 0.5) * depthMax * 1.1 / nCurvePts;
            fCurve[it] = GPShower::GammaPDF(tCurve[it], alpha, beta);
        }

        TGraph *gCurve = new TGraph(nCurvePts, tCurve, fCurve);
        gCurve->SetLineColor(hitColors[ih % 6]);
        gCurve->SetLineWidth(2);
        gCurve->SetLineStyle(2);  // dashed for analytic
        gCurve->Draw("L same");
    }

    cParams->Update();
}

///////////////////////////////////////////////////////////////
// Main macro function
///////////////////////////////////////////////////////////////

void MolPolDetectorResponseAnalysis(const char *fileList,
                                     Float_t energyCut = 0.0,
                                     const char *activeRowStr = "1234") {

    printf("=== MolPolDetectorResponseAnalysis ===\n");
    printf("File: %s\n", fileList);
    printf("Gate: PMT Response\n");
    printf("Energy cut: %.4f GeV\n", energyCut);

    // --- Parse active PMT rows ---
    Bool_t activeRow[4];
    ParseActiveRows(activeRowStr, activeRow);
    printf("Active PMT rows: ");
    for (Int_t r = 0; r < 4; r++) {
        if (activeRow[r]) printf("%d ", r + 1);
    }
    printf("\n");
    printf("Active segments: ");
    for (Int_t r = 0; r < 4; r++) {
        if (activeRow[r])
            printf("%s/%s ", DetGeom::PmtSegName[r], DetGeom::PmtSegName[4 + r]);
    }
    printf("\n");

    // --- Load data ---
    TChain *tree = SetupDetRespChain(fileList);
    if (!tree) return;
    SetupDetRespBranches(tree);

    Long64_t nEntries = tree->GetEntries();
    printf("Entries: %lld\n", nEntries);

    // --- Initialize accumulators ---
    DetRespAsymAccum accum[kDetRespAsymN];
    for (Int_t ch = 0; ch < kDetRespAsymN; ch++)
        DetRespAccumInit(accum[ch]);

    Double_t weightTotalThrown = 0.0;
    Double_t totalMollerEvXs = 0.0;

    gStyle->SetOptStat(111111);

    // --- First pass: find energy range for axis scaling ---
    Double_t maxPmtETotal_active = 0.0;
    Double_t maxPmtETotal_all    = 0.0;
    Double_t maxPmtESeg          = 0.0;
    Double_t maxPmtELeft         = 0.0;
    Double_t maxPmtERight        = 0.0;
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        Double_t eL_active = 0.0, eR_active = 0.0;
        Double_t eL_all = 0.0, eR_all = 0.0;
        for (Int_t r = 0; r < 4; r++) {
            eL_all += dr_pmtE[r];
            eR_all += dr_pmtE[4 + r];
            if (activeRow[r]) {
                eL_active += dr_pmtE[r];
                eR_active += dr_pmtE[4 + r];
            }
        }

        // Gate: at least one arm above threshold
        if (eL_active + eR_active > 0.0) {
            if (eL_active + eR_active > maxPmtETotal_active)
                maxPmtETotal_active = eL_active + eR_active;
            if (eL_all + eR_all > maxPmtETotal_all)
                maxPmtETotal_all = eL_all + eR_all;
            if (eL_active > maxPmtELeft)  maxPmtELeft  = eL_active;
            if (eR_active > maxPmtERight) maxPmtERight = eR_active;
            for (Int_t s = 0; s < 8; s++) {
                if (dr_pmtE[s] > maxPmtESeg) maxPmtESeg = dr_pmtE[s];
            }
        }
    }

    Double_t eMaxSeg = ceil(maxPmtESeg * 2.0) / 2.0;
    Double_t eMaxSum = ceil(maxPmtETotal_all * 2.0) / 2.0;
    Double_t eMaxLR  = ceil(std::max(maxPmtELeft, maxPmtERight) * 2.0) / 2.0;
    printf("Max pmtETotal (active): %.3f GeV\n", maxPmtETotal_active);
    printf("Max pmtETotal (all):    %.3f GeV\n", maxPmtETotal_all);
    printf("Max arm (L/R):          %.3f / %.3f GeV\n", maxPmtELeft, maxPmtERight);

    // --- Create histograms ---
    const Int_t nBins = 200;

    // Per-segment energy spectra
    TH1F *hPmtSeg[8];
    TH1F *hPmtSegAbove[8];
    TH1F *hPmtSegBelow[8];
    for (Int_t s = 0; s < 8; s++) {
        hPmtSeg[s] = new TH1F(
            TString::Format("hPmt%s", DetGeom::PmtSegName[s]).Data(),
            TString::Format("%s Energy;E [GeV];Events", DetGeom::PmtSegName[s]).Data(),
            nBins, 0.0, eMaxSeg);
        hPmtSeg[s]->SetLineWidth(2);

        hPmtSegAbove[s] = new TH1F(
            TString::Format("hPmt%sAbove", DetGeom::PmtSegName[s]).Data(),
            TString::Format("%s Energy - %s Sum Cut: %.2f GeV;E [GeV];Events",
                            DetGeom::PmtSegName[s],
                            (s < 4) ? "Left" : "Right", energyCut).Data(),
            nBins, 0.0, eMaxSeg);
        hPmtSegAbove[s]->SetLineColor(kBlue);
        hPmtSegAbove[s]->SetLineWidth(2);

        hPmtSegBelow[s] = new TH1F(
            TString::Format("hPmt%sBelow", DetGeom::PmtSegName[s]).Data(),
            TString::Format("%s Energy - %s Sum Cut: %.2f GeV;E [GeV];Events",
                            DetGeom::PmtSegName[s],
                            (s < 4) ? "Left" : "Right", energyCut).Data(),
            nBins, 0.0, eMaxSeg);
        hPmtSegBelow[s]->SetLineColor(kMagenta);
        hPmtSegBelow[s]->SetLineWidth(2);
        hPmtSegBelow[s]->SetLineStyle(2);
    }

    // Per-segment: coincidence-gated above/below
    TH1F *hPmtSegAboveCoinc[8];
    TH1F *hPmtSegBelowCoinc[8];
    for (Int_t s = 0; s < 8; s++) {
        hPmtSegAboveCoinc[s] = new TH1F(
            TString::Format("hPmt%sAboveCoinc", DetGeom::PmtSegName[s]).Data(),
            TString::Format("%s Energy - Threshold & Coinc Cut: %.2f GeV;E [GeV];Events",
                            DetGeom::PmtSegName[s], energyCut).Data(),
            nBins, 0.0, eMaxSeg);
        hPmtSegAboveCoinc[s]->SetLineColor(kCyan+1);
        hPmtSegAboveCoinc[s]->SetLineWidth(2);

        hPmtSegBelowCoinc[s] = new TH1F(
            TString::Format("hPmt%sBelowCoinc", DetGeom::PmtSegName[s]).Data(),
            TString::Format("%s Energy - Threshold & Coinc Cut: %.2f GeV;E [GeV];Events",
                            DetGeom::PmtSegName[s], energyCut).Data(),
            nBins, 0.0, eMaxSeg);
        hPmtSegBelowCoinc[s]->SetLineColor(kMagenta);
        hPmtSegBelowCoinc[s]->SetLineWidth(2);
        hPmtSegBelowCoinc[s]->SetLineStyle(2);
    }

    // Sum histograms: active rows and all rows
    const char *sumLabels[3] = { "Left", "Right", "Total" };
    TH1F *hSumActive[3];
    TH1F *hSumAll[3];
    TH1F *hSumActiveCut[3];
    TH1F *hSumAllCut[3];
    for (Int_t i = 0; i < 3; i++) {
        hSumActive[i] = new TH1F(
            TString::Format("hSum%sActive", sumLabels[i]).Data(),
            TString::Format("%s Energy;E [GeV];Events", sumLabels[i]).Data(),
            nBins, 0.0, eMaxSum);
        hSumActive[i]->SetLineColor(kBlue);
        hSumActive[i]->SetLineWidth(2);

        hSumAll[i] = new TH1F(
            TString::Format("hSum%sAll", sumLabels[i]).Data(),
            TString::Format("%s Energy;E [GeV];Events", sumLabels[i]).Data(),
            nBins, 0.0, eMaxSum);
        hSumAll[i]->SetLineColor(kRed);
        hSumAll[i]->SetLineWidth(2);
        hSumAll[i]->SetLineStyle(2);

        hSumActiveCut[i] = new TH1F(
            TString::Format("hSum%sActiveCut", sumLabels[i]).Data(),
            TString::Format("%s Energy - Threshold & Coincidence Cut: %.2f GeV;E [GeV];Events",
                            sumLabels[i], energyCut).Data(),
            nBins, 0.0, eMaxSum);
        hSumActiveCut[i]->SetLineColor(kCyan+1);
        hSumActiveCut[i]->SetLineWidth(2);

        hSumAllCut[i] = new TH1F(
            TString::Format("hSum%sAllCut", sumLabels[i]).Data(),
            TString::Format("%s Energy - Threshold & Coincidence Cut: %.2f GeV;E [GeV];Events",
                            sumLabels[i], energyCut).Data(),
            nBins, 0.0, eMaxSum);
        hSumAllCut[i]->SetLineColor(kMagenta);
        hSumAllCut[i]->SetLineWidth(2);
        hSumAllCut[i]->SetLineStyle(2);
    }

    // Sum histograms: single-arm threshold cut
    TH1F *hSumActiveSingle[3];
    TH1F *hSumAllSingle[3];
    for (Int_t i = 0; i < 3; i++) {
        hSumActiveSingle[i] = new TH1F(
            TString::Format("hSum%sActiveSingle", sumLabels[i]).Data(),
            TString::Format("%s Energy - %s Arm Cut: %.2f GeV;E [GeV];Events",
                            sumLabels[i],
                            (i < 2) ? ((i == 0) ? "Left" : "Right") : "Either",
                            energyCut).Data(),
            nBins, 0.0, eMaxSum);
        hSumActiveSingle[i]->SetLineColor(kBlue);
        hSumActiveSingle[i]->SetLineWidth(2);

        hSumAllSingle[i] = new TH1F(
            TString::Format("hSum%sAllSingle", sumLabels[i]).Data(),
            TString::Format("%s Energy - %s Arm Cut: %.2f GeV;E [GeV];Events",
                            sumLabels[i],
                            (i < 2) ? ((i == 0) ? "Left" : "Right") : "Either",
                            energyCut).Data(),
            nBins, 0.0, eMaxSum);
        hSumAllSingle[i]->SetLineColor(kRed);
        hSumAllSingle[i]->SetLineWidth(2);
        hSumAllSingle[i]->SetLineStyle(2);
    }

    // L-vs-R heatmap (active rows, square aspect)
    TH2F *hLR = new TH2F("hLR",
        "Left vs Right (Active PMT);pmtELeft [GeV];pmtERight [GeV]",
        nBins, 0.0, eMaxLR, nBins, 0.0, eMaxLR);

    // A_zz histograms (coincidence gets double bins for finer resolution)
    const Double_t azzMax = 7.0 / 9.0;
    TH1F *hAzz[kDetRespAsymN];
    hAzz[kDetRespAsymCoinc] = new TH1F("hAzzCoinc",
        "A_{zz} Coincidence;A_{zz};Weighted counts", 400, 0.0, azzMax);
    hAzz[kDetRespAsymLeft] = new TH1F("hAzzLeft",
        "A_{zz} Left singles;A_{zz};Weighted counts", 200, 0.0, azzMax);
    hAzz[kDetRespAsymRight] = new TH1F("hAzzRight",
        "A_{zz} Right singles;A_{zz};Weighted counts", 200, 0.0, azzMax);
    for (Int_t ch = 0; ch < kDetRespAsymN; ch++)
        hAzz[ch]->SetLineWidth(2);

    // Energy ratio diagnostic histograms
    TH1F *hRatioAll  = new TH1F("hRatioAll",
        "All events;E_{det}/E_{actual};Events", 300, 0.0, 1.5);
    TH1F *hRatio1hit = new TH1F("hRatio1hit",
        "1 hit;E_{det}/E_{actual};Events", 300, 0.0, 1.5);
    TH1F *hRatio2hit = new TH1F("hRatio2hit",
        "2 hits;E_{det}/E_{actual};Events", 300, 0.0, 1.5);
    TH1F *hRatio3plus = new TH1F("hRatio3plus",
        "3+ hits;E_{det}/E_{actual};Events", 300, 0.0, 1.5);
    hRatioAll->SetLineColor(kBlack);     hRatioAll->SetLineWidth(2);
    hRatio1hit->SetLineColor(kRed);      hRatio1hit->SetLineWidth(2);
    hRatio2hit->SetLineColor(kBlue);     hRatio2hit->SetLineWidth(2);
    hRatio3plus->SetLineColor(kGreen+2); hRatio3plus->SetLineWidth(2);

    TH2F *hRatioVsHitZ = new TH2F("hRatioVsHitZ",
        "E_{det}/E_{actual} vs Shower Start Depth (per-hit);hitZ [cm];E_{det}/E_{actual}",
        240, 0.0, 12.0, 200, 0.0, 1.5);
    TH2F *hRatioVsHitE = new TH2F("hRatioVsHitE",
        "E_{det}/E_{actual} vs Hit Energy (per-hit);E_{hit} [GeV];E_{det}/E_{actual}",
        240, 0.0, 12.0, 200, 0.0, 1.5);
    TH2F *hRatioVsSumE = new TH2F("hRatioVsSumE",
        "E_{det}/E_{actual} vs Event Energy (per-event);E_{actual} [GeV];E_{det}/E_{actual}",
        240, 0.0, 12.0, 200, 0.0, 1.5);

    // --- Main event loop ---
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        tree->GetEntry(iEntry);

        weightTotalThrown += dr_evUnpolWght;
        totalMollerEvXs   += dr_evXs;

        // Compute energy sums
        Double_t pmtELeft_active  = 0.0, pmtERight_active  = 0.0;
        Double_t pmtELeft_all     = 0.0, pmtERight_all     = 0.0;
        for (Int_t r = 0; r < 4; r++) {
            pmtELeft_all  += dr_pmtE[r];
            pmtERight_all += dr_pmtE[4 + r];
            if (activeRow[r]) {
                pmtELeft_active  += dr_pmtE[r];
                pmtERight_active += dr_pmtE[4 + r];
            }
        }

        // Gate: pass if at least one arm (active) has energy above threshold
        Bool_t hasActivity = (pmtELeft_active + pmtERight_active > 0.0);

        // Channel gates (using active-row sums)
        Bool_t passCoinc = (pmtELeft_active > energyCut &&
                            pmtERight_active > energyCut);
        Bool_t passLeft  = (pmtELeft_active > energyCut);
        Bool_t passRight = (pmtERight_active > energyCut);

        // Fill PMT plots for events with any activity
        if (hasActivity) {
            // Per-segment spectra
            for (Int_t s = 0; s < 8; s++) {
                if (dr_pmtE[s] > 0.0) {
                    hPmtSeg[s]->Fill(dr_pmtE[s]);
                    // Per-arm gating: left segs (0-3) by passLeft,
                    // right segs (4-7) by passRight
                    Bool_t passArm = (s < 4) ? passLeft : passRight;
                    if (passArm)
                        hPmtSegAbove[s]->Fill(dr_pmtE[s]);
                    else
                        hPmtSegBelow[s]->Fill(dr_pmtE[s]);
                    // Coincidence gating
                    if (passCoinc)
                        hPmtSegAboveCoinc[s]->Fill(dr_pmtE[s]);
                    else
                        hPmtSegBelowCoinc[s]->Fill(dr_pmtE[s]);
                }
            }

            // Sum histograms
            if (pmtELeft_active  > 0.0) hSumActive[0]->Fill(pmtELeft_active);
            if (pmtERight_active > 0.0) hSumActive[1]->Fill(pmtERight_active);
            if (pmtELeft_active + pmtERight_active > 0.0)
                hSumActive[2]->Fill(pmtELeft_active + pmtERight_active);

            if (pmtELeft_all  > 0.0) hSumAll[0]->Fill(pmtELeft_all);
            if (pmtERight_all > 0.0) hSumAll[1]->Fill(pmtERight_all);
            if (pmtELeft_all + pmtERight_all > 0.0)
                hSumAll[2]->Fill(pmtELeft_all + pmtERight_all);

            // L-vs-R heatmap (coincidence gate: both arms above energyCut)
        }
        if (passCoinc) {
            hLR->Fill(pmtELeft_active, pmtERight_active);
        }

        // Fill threshold-cut sum histograms
        if (passLeft) {
            hSumActiveCut[0]->Fill(pmtELeft_active);
            if (pmtELeft_all > 0.0) hSumAllCut[0]->Fill(pmtELeft_all);
        }
        if (passRight) {
            hSumActiveCut[1]->Fill(pmtERight_active);
            if (pmtERight_all > 0.0) hSumAllCut[1]->Fill(pmtERight_all);
        }
        if (passCoinc) {
            hSumActiveCut[2]->Fill(pmtELeft_active + pmtERight_active);
            if (pmtELeft_all + pmtERight_all > 0.0)
                hSumAllCut[2]->Fill(pmtELeft_all + pmtERight_all);
        }

        // Single-arm threshold-cut sum histograms
        if (passLeft) {
            hSumActiveSingle[0]->Fill(pmtELeft_active);
            if (pmtELeft_all > 0.0) hSumAllSingle[0]->Fill(pmtELeft_all);
        }
        if (passRight) {
            hSumActiveSingle[1]->Fill(pmtERight_active);
            if (pmtERight_all > 0.0) hSumAllSingle[1]->Fill(pmtERight_all);
        }
        if (passLeft || passRight) {
            hSumActiveSingle[2]->Fill(pmtELeft_active + pmtERight_active);
            if (pmtELeft_all + pmtERight_all > 0.0)
                hSumAllSingle[2]->Fill(pmtELeft_all + pmtERight_all);
        }

        // Energy ratio diagnostics (throw 0)
        if (dr_nHitsDet9 > 0) {
            Double_t sumHitE = 0.0;
            for (Int_t h = 0; h < dr_nHitsDet9; h++)
                sumHitE += dr_hitE[h];

            if (sumHitE > 0.0) {
                Double_t detTotal = 0.0;
                for (Int_t s = 0; s < 8; s++) detTotal += dr_pmtE[s];
                Double_t ratio = detTotal / sumHitE;

                hRatioAll->Fill(ratio);
                if (dr_nHitsDet9 == 1)      hRatio1hit->Fill(ratio);
                else if (dr_nHitsDet9 == 2) hRatio2hit->Fill(ratio);
                else                        hRatio3plus->Fill(ratio);

                hRatioVsHitZ->Fill(dr_hitZ[0], ratio);
                hRatioVsSumE->Fill(sumHitE, ratio);
                if (dr_nHitsDet9 == 1)
                    hRatioVsHitE->Fill(dr_hitE[0], ratio);
            }
        }

        // A_zz accumulation (all three channels independently)
        Double_t eventAzz = DetRespEventAsymmetry();

        if (passCoinc) {
            DetRespAccumulate(accum[kDetRespAsymCoinc], eventAzz);
            if (eventAzz > 0.0)
                hAzz[kDetRespAsymCoinc]->Fill(eventAzz,
                    dr_evPolPlusWghtZ + dr_evPolMinusWghtZ);
        }
        if (passLeft) {
            DetRespAccumulate(accum[kDetRespAsymLeft], eventAzz);
            if (eventAzz > 0.0)
                hAzz[kDetRespAsymLeft]->Fill(eventAzz,
                    dr_evPolPlusWghtZ + dr_evPolMinusWghtZ);
        }
        if (passRight) {
            DetRespAccumulate(accum[kDetRespAsymRight], eventAzz);
            if (eventAzz > 0.0)
                hAzz[kDetRespAsymRight]->Fill(eventAzz,
                    dr_evPolPlusWghtZ + dr_evPolMinusWghtZ);
        }
    }

    // --- Summary printout ---
    printf("\n========= Summary =========\n");
    printf("Entries: %lld\n", nEntries);

    const Double_t targetPolarization = 0.08014;

    for (Int_t ch = 0; ch < kDetRespAsymN; ch++) {
        printf("--- %s ---\n", DetRespAsymChannelName(ch));

        Double_t acceptedEventFrac = 0.0;
        if (nEntries > 0)
            acceptedEventFrac = (Double_t)accum[ch].nPass / (Double_t)nEntries;
        Double_t weightedEventFrac = 0.0;
        if (weightTotalThrown > 0.0)
            weightedEventFrac = accum[ch].weightAccepted / weightTotalThrown;
        Double_t mollerWtdEventFrac = 0.0;
        if (totalMollerEvXs > 0.0)
            mollerWtdEventFrac = accum[ch].sumMollerEvXs / totalMollerEvXs;

        printf("      Raw event fraction: %.5f\n", acceptedEventFrac);
        printf(" Moller-wtd event fraction: %.5f\n", mollerWtdEventFrac);
        printf("   Weighted event fraction: %.5f\n", weightedEventFrac);

        if (accum[ch].nPass == 0) {
            printf("  (no passing events — asymmetry omitted)\n");
            continue;
        }

        Double_t totalSignalWeight = accum[ch].sumPolWtPos + accum[ch].sumPolWtNeg;
        if (totalSignalWeight <= 0.0) {
            printf("  (zero polarized weight — asymmetry omitted)\n");
            continue;
        }

        Double_t meanAzz  = accum[ch].summedAsym / totalSignalWeight;
        Double_t meanAzz2 = accum[ch].summedAsym2 / totalSignalWeight;

        Double_t polarizationFactor = 0.0;
        Double_t polarizationFacErr = 0.0;
        Double_t meanAnalyzingPower = 0.0;
        Double_t meanAnalyzingPowerErr = 0.0;

        if (meanAzz != 0.0) {
            polarizationFactor = 1.0 / meanAzz;
            Double_t varAzz = meanAzz2 - meanAzz * meanAzz;
            polarizationFacErr = sqrt(fabs(varAzz) / (Double_t)accum[ch].nPass)
                                 / fabs(meanAzz) * fabs(polarizationFactor);
            meanAnalyzingPower = 1.0 / polarizationFactor / targetPolarization;
            meanAnalyzingPowerErr = polarizationFacErr / polarizationFactor
                                   / polarizationFactor / targetPolarization;
        }

        printf("    Polarization factor: %.5f +/- %.5f\n",
               polarizationFactor, polarizationFacErr);
        printf("   Mean analyzing power: %.7f +/- %.7f\n",
               meanAnalyzingPower, meanAnalyzingPowerErr);
    }

    // =================================================================
    // PLOTS
    // =================================================================

    // --- Canvas 1: PMT segment energy spectra (4R x 2C) ---
    TString segCanvasTitle = TString::Format(
        "PMT Segment Energy - Per-Arm Threshold Cut: %.2f GeV", energyCut);

    TCanvas *cPmtSeg = new TCanvas("cPmtSeg", segCanvasTitle.Data(), 800, 900);
    cPmtSeg->Divide(2, 4, 0.01, 0.01);

    Int_t padMap[8] = { 1, 3, 5, 7, 2, 4, 6, 8 };

    for (Int_t s = 0; s < 8; s++) {
        cPmtSeg->cd(padMap[s]);
        Int_t row = s % 4;

        Double_t yMaxSeg = std::max(hPmtSegAbove[s]->GetMaximum(),
                                     hPmtSegBelow[s]->GetMaximum()) * 1.15;
        if (yMaxSeg <= 0.0) yMaxSeg = 1.0;

        if (!activeRow[row]) {
            gPad->SetFillColor(TColor::GetColor(1.0f, 0.9f, 0.9f));
            hPmtSegBelow[s]->SetMaximum(yMaxSeg);
            hPmtSegBelow[s]->SetMinimum(0.0);
            hPmtSegBelow[s]->Draw("HISTE");
            hPmtSegAbove[s]->Draw("HISTE same");
            TPaveText *ptInactive = new TPaveText(0.3, 0.4, 0.7, 0.6, "NDC");
            ptInactive->SetBorderSize(1);
            ptInactive->SetFillStyle(1001);
            ptInactive->SetFillColor(kWhite);
            ptInactive->SetTextColor(TColor::GetColor(139, 0, 0));
            ptInactive->AddText("Inactive PMT");
            ptInactive->Draw();
        } else {
            gPad->SetFillColor(kWhite);
            hPmtSegBelow[s]->SetMaximum(yMaxSeg);
            hPmtSegBelow[s]->SetMinimum(0.0);
            hPmtSegBelow[s]->Draw("HISTE");
            hPmtSegAbove[s]->Draw("HISTE same");
        }

        // Legend on first pad only
        if (s == 0) {
            TLegend *legSeg = new TLegend(0.45, 0.7, 0.88, 0.88);
            legSeg->AddEntry(hPmtSegAbove[s], "Above cut", "l");
            legSeg->AddEntry(hPmtSegBelow[s], "Below cut", "l");
            legSeg->SetBorderSize(0);
            legSeg->SetFillStyle(0);
            legSeg->Draw();
        }
    }
    cPmtSeg->Update();

    // --- Canvas 1b: PMT segment energy (log) ---
    TCanvas *cPmtSegLog = new TCanvas("cPmtSegLog",
        (segCanvasTitle + " (log)").Data(), 800, 900);
    cPmtSegLog->Divide(2, 4, 0.01, 0.01);

    for (Int_t s = 0; s < 8; s++) {
        cPmtSegLog->cd(padMap[s]);
        Int_t row = s % 4;
        gPad->SetLogy();

        if (!activeRow[row]) {
            gPad->SetFillColor(TColor::GetColor(1.0f, 0.9f, 0.9f));
        } else {
            gPad->SetFillColor(kWhite);
        }

        TH1F *hBelowLog = (TH1F*)hPmtSegBelow[s]->Clone(
            TString::Format("hPmtSegBelowLog_%d", s).Data());
        TH1F *hAboveLog = (TH1F*)hPmtSegAbove[s]->Clone(
            TString::Format("hPmtSegAboveLog_%d", s).Data());
        hBelowLog->SetMinimum(0.9);
        hBelowLog->Draw("HISTE");
        hAboveLog->Draw("HISTE same");

        if (!activeRow[row]) {
            TPaveText *ptInactive = new TPaveText(0.3, 0.4, 0.7, 0.6, "NDC");
            ptInactive->SetBorderSize(1);
            ptInactive->SetFillStyle(1001);
            ptInactive->SetFillColor(kWhite);
            ptInactive->SetTextColor(TColor::GetColor(139, 0, 0));
            ptInactive->AddText("Inactive PMT");
            ptInactive->Draw();
        }

        if (s == 0) {
            TLegend *legSegLog = new TLegend(0.45, 0.7, 0.88, 0.88);
            legSegLog->AddEntry(hAboveLog, "Above cut", "l");
            legSegLog->AddEntry(hBelowLog, "Below cut", "l");
            legSegLog->SetBorderSize(0);
            legSegLog->SetFillStyle(0);
            legSegLog->Draw();
        }
    }
    cPmtSegLog->Update();

    // --- Canvas 1c: PMT segment energy - Coincidence cut (linear) ---
    TString segCoincTitle = TString::Format(
        "PMT Segment Energy - Threshold & Coincidence Cut: %.2f GeV", energyCut);

    TCanvas *cPmtSegCoinc = new TCanvas("cPmtSegCoinc",
        segCoincTitle.Data(), 800, 900);
    cPmtSegCoinc->Divide(2, 4, 0.01, 0.01);

    for (Int_t s = 0; s < 8; s++) {
        cPmtSegCoinc->cd(padMap[s]);
        Int_t row = s % 4;

        Double_t yMaxSC = std::max(hPmtSegAboveCoinc[s]->GetMaximum(),
                                    hPmtSegBelowCoinc[s]->GetMaximum()) * 1.15;
        if (yMaxSC <= 0.0) yMaxSC = 1.0;

        if (!activeRow[row]) {
            gPad->SetFillColor(TColor::GetColor(1.0f, 0.9f, 0.9f));
        } else {
            gPad->SetFillColor(kWhite);
        }

        hPmtSegBelowCoinc[s]->SetMaximum(yMaxSC);
        hPmtSegBelowCoinc[s]->SetMinimum(0.0);
        hPmtSegBelowCoinc[s]->Draw("HISTE");
        hPmtSegAboveCoinc[s]->Draw("HISTE same");

        if (!activeRow[row]) {
            TPaveText *ptI = new TPaveText(0.3, 0.4, 0.7, 0.6, "NDC");
            ptI->SetBorderSize(1);
            ptI->SetFillStyle(1001);
            ptI->SetFillColor(kWhite);
            ptI->SetTextColor(TColor::GetColor(139, 0, 0));
            ptI->AddText("Inactive PMT");
            ptI->Draw();
        }

        if (s == 0) {
            TLegend *legSC = new TLegend(0.45, 0.7, 0.88, 0.88);
            legSC->AddEntry(hPmtSegAboveCoinc[s], "Above cut", "l");
            legSC->AddEntry(hPmtSegBelowCoinc[s], "Below cut", "l");
            legSC->SetBorderSize(0);
            legSC->SetFillStyle(0);
            legSC->Draw();
        }
    }
    cPmtSegCoinc->Update();

    // --- Canvas 1d: PMT segment energy - Coincidence cut (log) ---
    TCanvas *cPmtSegCoincLog = new TCanvas("cPmtSegCoincLog",
        (segCoincTitle + " (log)").Data(), 800, 900);
    cPmtSegCoincLog->Divide(2, 4, 0.01, 0.01);

    for (Int_t s = 0; s < 8; s++) {
        cPmtSegCoincLog->cd(padMap[s]);
        Int_t row = s % 4;
        gPad->SetLogy();

        if (!activeRow[row]) {
            gPad->SetFillColor(TColor::GetColor(1.0f, 0.9f, 0.9f));
        } else {
            gPad->SetFillColor(kWhite);
        }

        TH1F *hBCLog = (TH1F*)hPmtSegBelowCoinc[s]->Clone(
            TString::Format("hPmtSegBelowCoincLog_%d", s).Data());
        TH1F *hACLog = (TH1F*)hPmtSegAboveCoinc[s]->Clone(
            TString::Format("hPmtSegAboveCoincLog_%d", s).Data());
        hBCLog->SetMinimum(0.9);
        hBCLog->Draw("HISTE");
        hACLog->Draw("HISTE same");

        if (!activeRow[row]) {
            TPaveText *ptI = new TPaveText(0.3, 0.4, 0.7, 0.6, "NDC");
            ptI->SetBorderSize(1);
            ptI->SetFillStyle(1001);
            ptI->SetFillColor(kWhite);
            ptI->SetTextColor(TColor::GetColor(139, 0, 0));
            ptI->AddText("Inactive PMT");
            ptI->Draw();
        }

        if (s == 0) {
            TLegend *legSCL = new TLegend(0.45, 0.7, 0.88, 0.88);
            legSCL->AddEntry(hACLog, "Above cut", "l");
            legSCL->AddEntry(hBCLog, "Below cut", "l");
            legSCL->SetBorderSize(0);
            legSCL->SetFillStyle(0);
            legSCL->Draw();
        }
    }
    cPmtSegCoincLog->Update();

    // --- Canvas 2: PMT sum energy (3R x 2C: linear left, log right) ---
    TCanvas *cPmtSums = new TCanvas("cPmtSums", "PMT Energy Sums", 1000, 900);
    cPmtSums->Divide(2, 3, 0.01, 0.01);

    for (Int_t i = 0; i < 3; i++) {
        // Left column: linear
        cPmtSums->cd(i * 2 + 1);
        Double_t yMax = std::max(hSumAll[i]->GetMaximum(),
                                 hSumActive[i]->GetMaximum()) * 1.15;
        hSumAll[i]->SetMaximum(yMax);
        hSumAll[i]->SetMinimum(0.0);
        hSumAll[i]->Draw("HISTE");
        hSumActive[i]->Draw("HISTE same");

        TLegend *legLin = new TLegend(0.15, 0.75, 0.45, 0.88);
        legLin->AddEntry(hSumActive[i], "Active PMT", "l");
        legLin->AddEntry(hSumAll[i], "All PMT", "l");
        legLin->SetBorderSize(0);
        legLin->SetFillStyle(0);
        legLin->Draw();

        // Right column: log
        cPmtSums->cd(i * 2 + 2);
        gPad->SetLogy();
        // Clone for independent y-axis range in log pad
        TH1F *hAllLog = (TH1F*)hSumAll[i]->Clone(
            TString::Format("hSumAllLog_%d", i).Data());
        TH1F *hActiveLog = (TH1F*)hSumActive[i]->Clone(
            TString::Format("hSumActiveLog_%d", i).Data());
        hAllLog->SetMinimum(0.9);
        hAllLog->SetMaximum(yMax);
        hAllLog->Draw("HISTE");
        hActiveLog->Draw("HISTE same");

        TLegend *legLog = new TLegend(0.15, 0.75, 0.45, 0.88);
        legLog->AddEntry(hActiveLog, "Active PMT", "l");
        legLog->AddEntry(hAllLog, "All PMT", "l");
        legLog->SetBorderSize(0);
        legLog->SetFillStyle(0);
        legLog->Draw();
    }
    cPmtSums->Update();

    // --- Canvas 2b: PMT sum energy - Threshold cut (3R x 2C: linear left, log right) ---
    TString sumCutTitle = TString::Format(
        "PMT Energy Sums - Threshold & Coincidence Cut: %.2f GeV", energyCut);
    TCanvas *cPmtSumsCut = new TCanvas("cPmtSumsCut",
        sumCutTitle.Data(), 1000, 900);
    cPmtSumsCut->Divide(2, 3, 0.01, 0.01);

    for (Int_t i = 0; i < 3; i++) {
        // Left column: linear
        cPmtSumsCut->cd(i * 2 + 1);
        Double_t yMaxCut = std::max(hSumAllCut[i]->GetMaximum(),
                                    hSumActiveCut[i]->GetMaximum()) * 1.15;
        if (yMaxCut <= 0.0) yMaxCut = 1.0;
        hSumAllCut[i]->SetMaximum(yMaxCut);
        hSumAllCut[i]->SetMinimum(0.0);
        hSumAllCut[i]->Draw("HISTE");
        hSumActiveCut[i]->Draw("HISTE same");

        TLegend *legCutLin = new TLegend(0.15, 0.75, 0.45, 0.88);
        legCutLin->AddEntry(hSumActiveCut[i], "Active PMT", "l");
        legCutLin->AddEntry(hSumAllCut[i], "All PMT", "l");
        legCutLin->SetBorderSize(0);
        legCutLin->SetFillStyle(0);
        legCutLin->Draw();

        // Right column: log
        cPmtSumsCut->cd(i * 2 + 2);
        gPad->SetLogy();
        TH1F *hAllCutLog = (TH1F*)hSumAllCut[i]->Clone(
            TString::Format("hSumAllCutLog_%d", i).Data());
        TH1F *hActiveCutLog = (TH1F*)hSumActiveCut[i]->Clone(
            TString::Format("hSumActiveCutLog_%d", i).Data());
        hAllCutLog->SetMinimum(0.9);
        hAllCutLog->SetMaximum(yMaxCut);
        hAllCutLog->Draw("HISTE");
        hActiveCutLog->Draw("HISTE same");

        TLegend *legCutLog = new TLegend(0.15, 0.75, 0.45, 0.88);
        legCutLog->AddEntry(hActiveCutLog, "Active PMT", "l");
        legCutLog->AddEntry(hAllCutLog, "All PMT", "l");
        legCutLog->SetBorderSize(0);
        legCutLog->SetFillStyle(0);
        legCutLog->Draw();
    }
    cPmtSumsCut->Update();

    // --- Canvas 2c: PMT sum energy - Single arm cut (3R x 2C) ---
    TString sumSingleTitle = TString::Format(
        "PMT Energy Sums - Single Arm Cut: %.2f GeV", energyCut);
    TCanvas *cPmtSumsSingle = new TCanvas("cPmtSumsSingle",
        sumSingleTitle.Data(), 1000, 900);
    cPmtSumsSingle->Divide(2, 3, 0.01, 0.01);

    for (Int_t i = 0; i < 3; i++) {
        // Left column: linear
        cPmtSumsSingle->cd(i * 2 + 1);
        Double_t yMaxSng = std::max(hSumAllSingle[i]->GetMaximum(),
                                    hSumActiveSingle[i]->GetMaximum()) * 1.15;
        if (yMaxSng <= 0.0) yMaxSng = 1.0;
        hSumAllSingle[i]->SetMaximum(yMaxSng);
        hSumAllSingle[i]->SetMinimum(0.0);
        hSumAllSingle[i]->Draw("HISTE");
        hSumActiveSingle[i]->Draw("HISTE same");

        TLegend *legSngLin = new TLegend(0.15, 0.75, 0.45, 0.88);
        legSngLin->AddEntry(hSumActiveSingle[i], "Active PMT", "l");
        legSngLin->AddEntry(hSumAllSingle[i], "All PMT", "l");
        legSngLin->SetBorderSize(0);
        legSngLin->SetFillStyle(0);
        legSngLin->Draw();

        // Right column: log
        cPmtSumsSingle->cd(i * 2 + 2);
        gPad->SetLogy();
        TH1F *hAllSngLog = (TH1F*)hSumAllSingle[i]->Clone(
            TString::Format("hSumAllSingleLog_%d", i).Data());
        TH1F *hActiveSngLog = (TH1F*)hSumActiveSingle[i]->Clone(
            TString::Format("hSumActiveSingleLog_%d", i).Data());
        hAllSngLog->SetMinimum(0.9);
        hAllSngLog->SetMaximum(yMaxSng);
        hAllSngLog->Draw("HISTE");
        hActiveSngLog->Draw("HISTE same");

        TLegend *legSngLog = new TLegend(0.15, 0.75, 0.45, 0.88);
        legSngLog->AddEntry(hActiveSngLog, "Active PMT", "l");
        legSngLog->AddEntry(hAllSngLog, "All PMT", "l");
        legSngLog->SetBorderSize(0);
        legSngLog->SetFillStyle(0);
        legSngLog->Draw();
    }
    cPmtSumsSingle->Update();

    // --- Canvas 3: L-vs-R heatmap (linear) ---
    gStyle->SetPalette(55);

    TCanvas *cPmtLR = new TCanvas("cPmtLR", "Left vs Right", 700, 600);
    cPmtLR->SetRightMargin(0.15);
    hLR->SetStats(kFALSE);
    hLR->Draw("COLZ");
    cPmtLR->Update();

    // --- Canvas 4: L-vs-R heatmap (log-z) ---
    TCanvas *cPmtLRLog = new TCanvas("cPmtLRLog", "Left vs Right (log)", 700, 600);
    cPmtLRLog->SetRightMargin(0.15);
    cPmtLRLog->SetLogz();
    hLR->Draw("COLZ");
    cPmtLRLog->Update();

    // --- Canvas 5: A_zz histograms ---
    for (Int_t ch = 0; ch < kDetRespAsymN; ch++) {
        Int_t firstBin = 1;
        for (Int_t b = 1; b <= hAzz[ch]->GetNbinsX(); b++) {
            if (hAzz[ch]->GetBinContent(b) > 0.0) {
                firstBin = b;
                break;
            }
        }
        Int_t lowBin = firstBin - 10;
        if (lowBin < 1) lowBin = 1;
        Double_t xMin = hAzz[ch]->GetBinLowEdge(lowBin);
        hAzz[ch]->GetXaxis()->SetRangeUser(xMin, azzMax);
    }

    TCanvas *cAzz = new TCanvas("cAzz", "Analyzing Power", 1400, 400);
    cAzz->Divide(3, 1);

    for (Int_t ch = 0; ch < kDetRespAsymN; ch++) {
        cAzz->cd(ch + 1);
        hAzz[ch]->SetMinimum(0.0);
        hAzz[ch]->Draw("HISTE");
    }
    cAzz->Update();

    // --- Canvas: Energy Ratio Diagnostics (2x2) ---
    TCanvas *cEratio = new TCanvas("cEratio",
        "Energy Ratio Diagnostics", 1200, 1000);
    cEratio->Divide(2, 2);

    // Panel 1: ratio by multiplicity
    cEratio->cd(1);
    gPad->SetLogy();
    hRatioAll->SetMinimum(0.9);
    hRatioAll->SetTitle("E_{det}/E_{actual} by Hit Multiplicity");
    hRatioAll->Draw("HISTE");
    hRatio1hit->Draw("HISTE same");
    hRatio2hit->Draw("HISTE same");
    hRatio3plus->Draw("HISTE same");
    TLegend *legR = new TLegend(0.15, 0.65, 0.45, 0.88);
    legR->AddEntry(hRatioAll, "All", "l");
    legR->AddEntry(hRatio1hit, "1 hit", "l");
    legR->AddEntry(hRatio2hit, "2 hits", "l");
    legR->AddEntry(hRatio3plus, "3+ hits", "l");
    legR->SetBorderSize(0);
    legR->SetFillStyle(0);
    legR->Draw();

    // Panel 2: ratio vs hitE (per-hit, single-hit only)
    cEratio->cd(2);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    hRatioVsHitE->SetStats(kFALSE);
    hRatioVsHitE->Draw("COLZ");

    // Panel 3: ratio vs hitZ (per-hit)
    cEratio->cd(3);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    hRatioVsHitZ->SetStats(kFALSE);
    hRatioVsHitZ->Draw("COLZ");

    // Panel 4: ratio vs sumHitE (per-event)
    cEratio->cd(4);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    hRatioVsSumE->SetStats(kFALSE);
    hRatioVsSumE->Draw("COLZ");

    cEratio->Update();
}

#endif // MOLPOL_DETECTOR_RESPONSE_ANALYSIS_C