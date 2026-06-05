///////////////////////////////////////////////////////////////
// MolPolDetectorResponseShowerInspection.C
// Eric King 06/02/2026
//
// Interactive shower inspection for MolPol detector response
// output files containing a TShower diagnostic tree.
//
// Usage:
//   root -l 'MolPolDetectorResponseShowerInspection.C("input_detresponse.root")'
//
// Once loaded, call interactively at the ROOT prompt:
//   ListShowerEvents(N, M, energyCut) -- list N events starting from M
//   PlotEventShower(eventNum)         -- 3D shower visualization
//   PlotShowerParams(eventNum)        -- Grindhammer-Peters parameter plots
//
// eventNum is 1-based. Defaults to 1 if omitted.
///////////////////////////////////////////////////////////////

#ifndef MOLPOL_DETECTOR_RESPONSE_SHOWER_INSPECTION_C
#define MOLPOL_DETECTOR_RESPONSE_SHOWER_INSPECTION_C

#include "MolPolDetectorResponse.h"

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
#include <map>

///////////////////////////////////////////////////////////////
// Globals — set by main macro, used by all interactive
// functions.
///////////////////////////////////////////////////////////////

TTree *gTShower  = nullptr;
TTree *gTDetResp = nullptr;

///////////////////////////////////////////////////////////////
// EventIndexEntry — per-event arm energies cached at load
// time for fast lookup in ListShowerEvents and future use.
// Keyed by eventIdx (== TDetResp entry number).
///////////////////////////////////////////////////////////////

struct EventIndexEntry {
    Double_t pmtELeft;
    Double_t pmtERight;
    Double_t pmtE[8];  // individual segment energies (post-smearing)
};

// Ordered list of eventIdx values present in TShower (1-based
// position matches the eventNum used by PlotEventShower etc.).
// gEventIndex maps eventIdx -> arm energies from TDetResp.
std::vector<Int_t>               gShowerEventList;
std::map<Int_t, EventIndexEntry> gEventIndex;


///////////////////////////////////////////////////////////////
// Solve90PctRadius()
// Bisection solver for the radius containing 90% of the
// shower energy given the Grindhammer-Peters core/tail params.
///////////////////////////////////////////////////////////////

Double_t Solve90PctRadius(Double_t p, Double_t Rc, Double_t Rt) {
    Double_t rLo = 0.0;
    Double_t rHi = 20.0;  // well beyond any physical radius in RM
    for (Int_t iter = 0; iter < 60; iter++) {
        Double_t rMid = 0.5 * (rLo + rHi);
        Double_t cdf = GPShower::RadialCDFTotal(rMid, p, Rc, Rt);
        if (cdf < 0.9)
            rLo = rMid;
        else
            rHi = rMid;
    }
    return 0.5 * (rLo + rHi);
}


///////////////////////////////////////////////////////////////
// HitSlice — per-depth-slice data for shower visualization
///////////////////////////////////////////////////////////////

struct HitSlice {
    Float_t x, y;       // hit position [cm]
    Float_t depth_cm;   // depth into calorimeter [cm]
    Float_t energy;     // slice energy [GeV]
    Float_t r90;        // 90% containment radius [cm]
    Float_t rc, rt, p;  // Grindhammer params (RM units) for debug
    Int_t   hitIdx;
};


///////////////////////////////////////////////////////////////
// PlotEventShower()
// 3D visualization of the modeled shower for a specific event
// in the TShower diagnostic tree.
//
// eventNum is 1-based: PlotEventShower(1) plots the first
// unique event, PlotEventShower(10) the tenth.
//
// Each depth slice is drawn as a TPolyLine3D ring at the 90%
// energy containment radius. Ring color indicates deposited
// energy (blue=low, red=high). Line width scales with energy
// fraction.
//
// Axes: X = depth into calorimeter [cm] (beam enters from left)
//       Y = horizontal position on face [cm] (hitLx)
//       Z = vertical position on face [cm] (hitLy)
///////////////////////////////////////////////////////////////

void PlotEventShower(Int_t eventNum = 1) {

    if (!gTShower) {
        printf("PlotEventShower: gTShower is null. "
               "Was writeShowerDiag enabled when the file was produced?\n");
        return;
    }

    Long64_t nEntries = gTShower->GetEntries();
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

    gTShower->SetBranchAddress("eventIdx",    &sh_eventIdx);
    gTShower->SetBranchAddress("hitIdx",      &sh_hitIdx);
    gTShower->SetBranchAddress("hitX",        &sh_hitX);
    gTShower->SetBranchAddress("hitY",        &sh_hitY);
    gTShower->SetBranchAddress("hitE",        &sh_hitE);
    gTShower->SetBranchAddress("nDepthSteps", &sh_nDepthSteps);
    gTShower->SetBranchAddress("depthT",      sh_depthT);
    gTShower->SetBranchAddress("sliceE",      sh_sliceE);
    gTShower->SetBranchAddress("RC",          sh_RC);
    gTShower->SetBranchAddress("RT",          sh_RT);
    gTShower->SetBranchAddress("p",           sh_p);

    // --- Build ordered list of unique eventIdx values ---
    std::vector<Int_t> uniqueList;
    for (Long64_t i = 0; i < nEntries; i++) {
        gTShower->GetEntry(i);
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
        gTShower->GetEntry(i);
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

    hCbar->GetXaxis()->SetLabelSize(0);
    hCbar->GetXaxis()->SetTickLength(0);
    hCbar->GetYaxis()->SetLabelSize(0);
    hCbar->GetYaxis()->SetTickLength(0);
    hCbar->GetYaxis()->SetTitle("");
    hCbar->Draw("COL");

    padCbar->Update();

    TGaxis *rightAxis = new TGaxis(
        1.0, 0.0, 1.0, maxEnergy,
        0.0, maxEnergy,
        510, "+L");
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
// PlotShowerParams()
// Four-panel plot of the Grindhammer-Peters shower parameters
// for a specific event:
//   Panel 1: RC vs tau (core radius)
//   Panel 2: RT vs tau (tail radius)
//   Panel 3: p vs tau (core weight)
//   Panel 4: 1/E dE/dt vs t (longitudinal profile)
//
// Panels 1-3 reproduce the style of Figs. 9-11 in
// arXiv:hep-ex/0001020. Panel 4 reproduces Fig. 4.
// All hits overlaid, distinguished by color.
// Uses the same 1-based event indexing as PlotEventShower().
///////////////////////////////////////////////////////////////

void PlotShowerParams(Int_t eventNum = 1) {
 
    if (!gTShower) {
        printf("PlotShowerParams: gTShower is null. "
               "Was writeShowerDiag enabled when the file was produced?\n");
        return;
    }
 
    Long64_t nEntries = gTShower->GetEntries();
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
 
    gTShower->SetBranchAddress("eventIdx",    &sh_eventIdx);
    gTShower->SetBranchAddress("hitIdx",      &sh_hitIdx);
    gTShower->SetBranchAddress("hitX",        &sh_hitX);
    gTShower->SetBranchAddress("hitY",        &sh_hitY);
    gTShower->SetBranchAddress("hitE",        &sh_hitE);
    gTShower->SetBranchAddress("nDepthSteps", &sh_nDepthSteps);
    gTShower->SetBranchAddress("depthT",      sh_depthT);
    gTShower->SetBranchAddress("sliceE",      sh_sliceE);
    gTShower->SetBranchAddress("RC",          sh_RC);
    gTShower->SetBranchAddress("RT",          sh_RT);
    gTShower->SetBranchAddress("p",           sh_p);
 
    // --- Build unique event list ---
    std::vector<Int_t> uniqueList;
    for (Long64_t i = 0; i < nEntries; i++) {
        gTShower->GetEntry(i);
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
 
    std::vector<Int_t>   hitIdxList;
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
        gTShower->GetEntry(i);
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
 
        if ((Int_t)hitTList.size() <= currentHitSlot)
            hitTList.push_back(T);
 
        for (Int_t s = 0; s < sh_nDepthSteps; s++) {
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
    Double_t tauMax   = 0.0;
    Double_t rcMax    = 0.0;
    Double_t rtMax    = 0.0;
    Double_t depthMax = 0.0;
    for (Int_t ih = 0; ih < nHits; ih++) {
        for (size_t s = 0; s < tauVecs[ih].size(); s++) {
            if (tauVecs[ih][s]   > tauMax)   tauMax   = tauVecs[ih][s];
            if (rcVecs[ih][s]    > rcMax)    rcMax    = rcVecs[ih][s];
            if (rtVecs[ih][s]    > rtMax)    rtMax    = rtVecs[ih][s];
            if (depthVecs[ih][s] > depthMax) depthMax = depthVecs[ih][s];
        }
    }
 
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
 
    TLegend *leg = new TLegend(0.125, 0.525, 0.5, 0.875);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetTextSize(0.040);
    leg->SetTextFont(42);
 
    for (Int_t ih = 0; ih < nHits; ih++) {
        Int_t nPts = (Int_t)tauVecs[ih].size();
        if (nPts == 0) continue;
        TGraph *g = new TGraph(nPts,
            &tauVecs[ih][0], &rcVecs[ih][0]);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.8);
        g->SetMarkerColor(hitColors[ih % 6]);
        g->SetLineColor(hitColors[ih % 6]);
        g->Draw("PL same");
        leg->AddEntry(g,
            TString::Format("Hit %d  E=%.3f GeV",
                hitIdxList[ih], hitEList[ih]).Data(), "P");
    }
    leg->Draw();
 
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
        g->SetMarkerSize(0.8);
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
        g->SetMarkerSize(0.8);
        g->SetMarkerColor(hitColors[ih % 6]);
        g->SetLineColor(hitColors[ih % 6]);
        g->Draw("PL same");
    }
 
    // --- Panel 4: Longitudinal profile 1/E dE/dt vs t ---
    // Data points: sliceWt = sliceE/hitE per 1 X0 slice.
    // Analytic curve: gamma distribution (Eq. 2, arXiv:hep-ex/0001020).
    cParams->cd(4);
 
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
 
        TGraph *gData = new TGraph(nPts,
            &depthVecs[ih][0], &sliceWtVecs[ih][0]);
        gData->SetMarkerStyle(20);
        gData->SetMarkerSize(0.8);
        gData->SetMarkerColor(hitColors[ih % 6]);
        gData->SetLineColor(hitColors[ih % 6]);
        gData->Draw("P same");
 
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
// BuildEventIndex()
// Called once at load time by the main macro. Walks TShower to
// build gShowerEventList (ordered unique eventIdx values), then
// looks each one up in TDetResp to cache pmtELeft/pmtERight in
// gEventIndex.
//
// NOTE: Every eventIdx in TShower must have a corresponding
// entry in TDetResp -- this is guaranteed by the output writer
// in MolPolDetectorResponse.C, where TDetResp and TShower are
// filled in the same event loop. A mismatch should never occur
// under normal operation and indicates file corruption or a
// logic error upstream.
///////////////////////////////////////////////////////////////

void BuildEventIndex() {

    gShowerEventList.clear();
    gEventIndex.clear();

    if (!gTShower || !gTDetResp) return;

    // --- Build ordered unique eventIdx list from TShower ---
    Int_t sh_eventIdx;
    gTShower->SetBranchAddress("eventIdx", &sh_eventIdx);

    Long64_t nShower = gTShower->GetEntries();
    for (Long64_t i = 0; i < nShower; i++) {
        gTShower->GetEntry(i);
        Bool_t found = false;
        for (size_t u = 0; u < gShowerEventList.size(); u++) {
            if (gShowerEventList[u] == sh_eventIdx) { found = true; break; }
        }
        if (!found) gShowerEventList.push_back(sh_eventIdx);
    }

    // --- Look up arm energies and segment deposits for each eventIdx in TDetResp ---
    // Uses dr_* globals from SetupDetRespBranches(), reading throw 0 values.
    Long64_t nDetResp = gTDetResp->GetEntries();

    for (size_t u = 0; u < gShowerEventList.size(); u++) {
        Int_t evIdx = gShowerEventList[u];

        // eventIdx == TDetResp entry number (0-based)
        if (evIdx < 0 || evIdx >= (Int_t)nDetResp) {
            // This should never happen -- see note above.
            printf("WARNING: BuildEventIndex: eventIdx %d has no corresponding "
                   "TDetResp entry (nDetResp = %lld). File may be corrupt.\n",
                   evIdx, nDetResp);
            EventIndexEntry e;
            e.pmtELeft = 0.0; e.pmtERight = 0.0;
            for (Int_t s = 0; s < 8; s++) e.pmtE[s] = 0.0;
            gEventIndex[evIdx] = e;
            continue;
        }

        gTDetResp->GetEntry(evIdx);
        EventIndexEntry e;
        e.pmtELeft  = (Double_t)dr_pmtELeft[0];    // throw 0
        e.pmtERight = (Double_t)dr_pmtERight[0];   // throw 0
        for (Int_t s = 0; s < 8; s++) e.pmtE[s] = (Double_t)dr_pmtE[s];  // throw 0, segs 0-7
        gEventIndex[evIdx] = e;
    }

    printf("BuildEventIndex: indexed %d shower events.\n",
           (Int_t)gShowerEventList.size());
}


///////////////////////////////////////////////////////////////
// ListShowerEvents()
// Lists N shower events starting from position M (1-based),
// printing eventNum, eventIdx, and gate type determined by
// arm energy sums against energyCut.
//
// Gate types (same definitions as MolPolDetectorResponseAnalysis):
//   Coincidence : pmtELeft  >= cut AND pmtERight >= cut
//   Left        : pmtELeft  >= cut AND pmtERight <  cut
//   Right       : pmtERight >= cut AND pmtELeft  <  cut
//   None        : both arms below cut
//
// Arguments:
//   N          -- number of events to list
//   M          -- starting event number, 1-based (default 1)
//   energyCut  -- arm energy threshold in GeV (default 0.0)
///////////////////////////////////////////////////////////////

void ListShowerEvents(Int_t N, Int_t M = 1, Float_t energyCut = 0.0) {

    if (gShowerEventList.empty()) {
        printf("ListShowerEvents: event index is empty. "
               "Was the file loaded successfully?\n");
        return;
    }

    Int_t nTotal = (Int_t)gShowerEventList.size();

    if (M < 1 || M > nTotal) {
        printf("ListShowerEvents: starting position M=%d out of range [1, %d]\n",
               M, nTotal);
        return;
    }

    Int_t iEnd = M - 1 + N;
    if (iEnd > nTotal) {
        printf("ListShowerEvents: requested %d events from position %d, "
               "but only %d available. Listing to end.\n",
               N, M, nTotal - (M - 1));
        iEnd = nTotal;
    }

    printf("\n");
    printf("  %-10s  %-10s  %s\n", "eventNum", "eventIdx", "Type");
    printf("  %-10s  %-10s  %s\n", "--------", "--------", "----");

    for (Int_t i = M - 1; i < iEnd; i++) {
        Int_t evIdx = gShowerEventList[i];
        Int_t evNum = i + 1;

        std::map<Int_t, EventIndexEntry>::iterator it = gEventIndex.find(evIdx);
        if (it == gEventIndex.end()) {
            // Should never happen if BuildEventIndex ran cleanly.
            printf("  %-10d  %-10d  [index missing -- see BuildEventIndex warning]\n",
                   evNum, evIdx);
            continue;
        }

        Double_t eL = it->second.pmtELeft;
        Double_t eR = it->second.pmtERight;

        const char *type;
        if      (eL >= energyCut && eR >= energyCut) type = "Coincidence";
        else if (eL >= energyCut)                    type = "Left";
        else if (eR >= energyCut)                    type = "Right";
        else                                         type = "None";

        printf("  %-10d  %-10d  %s\n", evNum, evIdx, type);
    }

    printf("\n");
}



///////////////////////////////////////////////////////////////
// DrawPmtLabels()
// Draws a 4x2 grid of TPaveText boxes showing post-smearing
// energy deposit per PMT segment.
//
// Layout as seen from upstream:
//   Left column  (x >= 0): L1 (top) -> L4 (bottom)  [segs 0-3]
//   Right column (x <  0): R1 (top) -> R4 (bottom)  [segs 4-7]
//
// xLo,xHi,yLo,yHi define the grid region in normalized pad
// coordinates. Boxes share edges for a clean grid appearance.
// Reusable: can be called from PlotShowerFace, PlotEventShower,
// or any other function with a valid TVirtualPad.
///////////////////////////////////////////////////////////////

void DrawPmtLabels(Int_t eventIdx, TVirtualPad *pad,
                   Double_t xLo, Double_t xHi,
                   Double_t yLo, Double_t yHi) {

    if (!pad) return;

    std::map<Int_t, EventIndexEntry>::iterator it = gEventIndex.find(eventIdx);
    if (it == gEventIndex.end()) {
        printf("DrawPmtLabels: eventIdx %d not in index.\n", eventIdx);
        return;
    }

    const Double_t *pmtE = it->second.pmtE;

    // Segment indexing: col*4 + row
    //   col 0 = Left  (x >= 0): segs 0,1,2,3 = L1,L2,L3,L4
    //   col 1 = Right (x <  0): segs 4,5,6,7 = R1,R2,R3,R4
    const char *segName[8] = { "L1","L2","L3","L4","R1","R2","R3","R4" };

    Double_t boxW = (xHi - xLo) / 2.0;
    Double_t boxH = (yHi - yLo) / 4.0;

    pad->cd();

    for (Int_t col = 0; col < 2; col++) {
        for (Int_t row = 0; row < 4; row++) {
            Int_t seg = col * 4 + row;

            Double_t x1 = xLo + col * boxW;
            Double_t x2 = x1  + boxW;
            Double_t y2 = yHi - row * boxH;
            Double_t y1 = y2  - boxH;

            TPaveText *pt = new TPaveText(x1, y1, x2, y2, "NDC");
            pt->SetBorderSize(1);
            pt->SetFillColor(0);
            pt->SetTextSize(0.060);
            pt->SetTextFont(42);
            pt->AddText(segName[seg]);
            pt->AddText(TString::Format("%.4f GeV", pmtE[seg]).Data());
            pt->Draw();
        }
    }

    pad->Update();
}


///////////////////////////////////////////////////////////////
// PlotShowerFace()
// 2D contour map of shower energy density on the detector face
// (X x Y) for a specific event, integrated over all depth
// slices and all hits in the TShower tree.
//
// Histogram: TH2F, 1 mm bins, X: -9 to +9 cm, Y: -15 to +15 cm.
// Each bin accumulates E * f(r) * binW * binH [GeV] summed over
// all depth slices, where f(r) is the two-component
// Grindhammer-Peters radial PDF evaluated at the bin centre.
//
// Canvas layout:
//   Left  80% : TH2F density map with PMT grid overlaid
//   Right 20% :
//     Top    50% : 4x2 TPaveText PMT energy grid (DrawPmtLabels)
//     Bottom 50% : colorbar
//
// Uses the same 1-based event indexing as PlotEventShower().
///////////////////////////////////////////////////////////////

void PlotShowerFace(Int_t eventNum = 1) {

    if (!gTShower) {
        printf("PlotShowerFace: gTShower is null.\n");
        return;
    }

    Long64_t nEntries = gTShower->GetEntries();
    if (nEntries == 0) {
        printf("PlotShowerFace: TShower tree is empty.\n");
        return;
    }

    Int_t nUniqueEvents = (Int_t)gShowerEventList.size();
    if (eventNum < 1 || eventNum > nUniqueEvents) {
        printf("PlotShowerFace: eventNum %d out of range [1, %d]\n",
               eventNum, nUniqueEvents);
        return;
    }

    Int_t targetEventIdx = gShowerEventList[eventNum - 1];
    printf("PlotShowerFace: event #%d (eventIdx = %d)\n",
           eventNum, targetEventIdx);

    // --- Wire up TShower branches ---
    Int_t   sh_eventIdx, sh_hitIdx, sh_nDepthSteps;
    Float_t sh_hitX, sh_hitY, sh_hitE;
    Float_t sh_depthT[MAX_DEPTH_STEPS];
    Float_t sh_sliceE[MAX_DEPTH_STEPS];
    Float_t sh_RC[MAX_DEPTH_STEPS];
    Float_t sh_RT[MAX_DEPTH_STEPS];
    Float_t sh_p[MAX_DEPTH_STEPS];

    gTShower->SetBranchAddress("eventIdx",    &sh_eventIdx);
    gTShower->SetBranchAddress("hitIdx",      &sh_hitIdx);
    gTShower->SetBranchAddress("hitX",        &sh_hitX);
    gTShower->SetBranchAddress("hitY",        &sh_hitY);
    gTShower->SetBranchAddress("hitE",        &sh_hitE);
    gTShower->SetBranchAddress("nDepthSteps", &sh_nDepthSteps);
    gTShower->SetBranchAddress("depthT",      sh_depthT);
    gTShower->SetBranchAddress("sliceE",      sh_sliceE);
    gTShower->SetBranchAddress("RC",          sh_RC);
    gTShower->SetBranchAddress("RT",          sh_RT);
    gTShower->SetBranchAddress("p",           sh_p);

    // --- Build density histogram: 1 mm bins ---
    // X: -9 to +9 cm (180 bins), Y: -15 to +15 cm (300 bins)
    const Int_t    nBinsX = 180;
    const Int_t    nBinsY = 300;
    const Double_t xMin = -9.0,  xMax = 9.0;
    const Double_t yMin = -15.0, yMax = 15.0;
    const Double_t binW = (xMax - xMin) / nBinsX;  // 0.1 cm
    const Double_t binH = (yMax - yMin) / nBinsY;

    TString hName = TString::Format("hFace_ev%d", eventNum);
    TH2F *hFace = new TH2F(hName.Data(), "",
                            nBinsX, xMin, xMax,
                            nBinsY, yMin, yMax);
    hFace->SetStats(0);

    // --- Fill histogram ---
    // For each depth slice, evaluate the radial PDF at every bin
    // centre and accumulate E * f(r) * binW * binH into that bin.
    for (Long64_t i = 0; i < nEntries; i++) {
        gTShower->GetEntry(i);
        if (sh_eventIdx != targetEventIdx) continue;

        for (Int_t s = 0; s < sh_nDepthSteps; s++) {
            if (sh_sliceE[s] <= 0.0f) continue;
            if (!(sh_RC[s] >= 0.0f && sh_RC[s] < 10.0f)) continue;
            if (!(sh_RT[s] >= 0.0f && sh_RT[s] < 10.0f)) continue;
            if (!(sh_p[s]  >= 0.0f && sh_p[s]  <= 1.0f)) continue;

            Double_t p   = sh_p[s];
            Double_t Rc  = sh_RC[s] * PbConst::RM;  // RM -> cm
            Double_t Rt  = sh_RT[s] * PbConst::RM;
            Double_t Rc2 = Rc * Rc;
            Double_t Rt2 = Rt * Rt;
            Double_t E   = sh_sliceE[s];
            Double_t hx  = sh_hitX;
            Double_t hy  = sh_hitY;

            for (Int_t bx = 1; bx <= nBinsX; bx++) {
                Double_t x  = xMin + (bx - 0.5) * binW;
                Double_t dx = x - hx;
                Double_t dx2 = dx * dx;
                for (Int_t by = 1; by <= nBinsY; by++) {
                    Double_t y   = yMin + (by - 0.5) * binH;
                    Double_t dy  = y - hy;
                    Double_t r2  = dx2 + dy * dy;
                    Double_t r   = sqrt(r2);

                    Double_t fCore = (Rc > 0.0) ?
                        2.0 * r * Rc2 / ((r2 + Rc2) * (r2 + Rc2)) : 0.0;
                    Double_t fTail = (Rt > 0.0) ?
                        2.0 * r * Rt2 / ((r2 + Rt2) * (r2 + Rt2)) : 0.0;

                    Double_t val = E * (p * fCore + (1.0 - p) * fTail)
                                     * binW * binH;
                    if (val > 0.0)
                        hFace->AddBinContent(hFace->GetBin(bx, by), val);
                }
            }
        }
    }

    // --- Canvas ---
    TString cName = TString::Format("cFace_ev%d", eventNum);
    TCanvas *c = new TCanvas(cName.Data(),
        TString::Format("Event #%d Detector Face", eventNum).Data(),
        1100, 700);

    gStyle->SetPalette(55);
    gStyle->SetOptStat(0);

    // Left 80%: density map
    TPad *padMap = new TPad("padMap", "", 0.0, 0.0, 0.80, 1.0);
    padMap->SetRightMargin(0.02);
    padMap->SetLeftMargin(0.10);
    padMap->SetTopMargin(0.08);
    padMap->SetBottomMargin(0.10);
    padMap->Draw();
    padMap->cd();

    hFace->SetTitle(TString::Format(
        "Event #%d Shower Energy Density;X [cm];Y [cm]", eventNum).Data());
    hFace->GetXaxis()->SetTitleSize(0.05);
    hFace->GetXaxis()->SetTitleOffset(0.9);
    hFace->GetYaxis()->SetTitleSize(0.05);
    hFace->GetYaxis()->SetTitleOffset(0.9);
    hFace->Draw("COL");

    // Suppress default COLZ palette axis
    padMap->Update();
    TPaletteAxis *palette =
        (TPaletteAxis *)hFace->GetListOfFunctions()->FindObject("palette");
    if (palette) { palette->SetX1NDC(1.5); palette->SetX2NDC(1.6); }
    padMap->Modified();
    padMap->Update();

    // Overlay PMT segment grid (white dashed)
    for (Int_t iy = 0; iy < 5; iy++) {
        TLine *line = new TLine(xMin, DetGeom::rowYBound[iy],
                                xMax, DetGeom::rowYBound[iy]);
        line->SetLineColor(kWhite);
        line->SetLineWidth(1);
        line->SetLineStyle(2);
        line->Draw("same");
    }
    TLine *lrLine = new TLine(0.0, yMin, 0.0, yMax);
    lrLine->SetLineColor(kWhite);
    lrLine->SetLineWidth(1);
    lrLine->SetLineStyle(2);
    lrLine->Draw("same");
    padMap->Update();

    // Right side
    c->cd();

    // Top 50% of right 20%: PMT label grid
    TPad *padLabels = new TPad("padLabels", "", 0.80, 0.50, 1.0, 1.0);
    padLabels->SetMargin(0.02, 0.02, 0.02, 0.05);
    padLabels->Draw();
    DrawPmtLabels(targetEventIdx, padLabels, 0.02, 0.925, 0.05, 0.8375);

    // Colorbar drawn directly on the canvas in NDC coordinates:
    // color strip: x(0.85, 0.90), y(0.10, 0.45)
    // axis labels drawn to the right of the strip.
    c->cd();
    Double_t maxVal = hFace->GetMaximum();

    // Draw the color strip as a filled TH2F inside a TPad sized
    // exactly to the strip region in canvas NDC.
    TPad *padCbar = new TPad("padCbar", "", 0.83, 0.08, 0.98, 0.50);
    //padCbar->SetMargin(0.0, 0.0, 0.0, 0.0);
    padCbar->SetMargin(0.0, 0.55, 0.05, 0.0);
    padCbar->SetFillStyle(0);
    padCbar->Draw();
    padCbar->cd();

    TH2F *hCbar = new TH2F("hCbarFace", "",
                            1, 0.0, 1.0,
                            100, 0.0, maxVal);
    for (Int_t ib = 1; ib <= 100; ib++)
        hCbar->SetBinContent(1, ib, hCbar->GetYaxis()->GetBinCenter(ib));
    hCbar->GetXaxis()->SetLabelSize(0);
    hCbar->GetXaxis()->SetTickLength(0);
    hCbar->GetYaxis()->SetLabelSize(0);
    hCbar->GetYaxis()->SetTickLength(0);
    hCbar->Draw("COL");
    padCbar->Update();

    // TGaxis drawn on canvas in NDC, sitting just right of the strip
    //c->cd();
    //TGaxis *cbAxis = new TGaxis(0.90, 0.10, 0.90, 0.48,0.0, maxVal, 510, "+L");
    padCbar->cd();
    TGaxis *cbAxis = new TGaxis(1.0, 0.0, 1.0, maxVal ,0.0, maxVal, 510, "+L");
    cbAxis->SetLineWidth(2);
    cbAxis->SetTitle("Energy [GeV/bin]");
    cbAxis->SetTitleOffset(1.8);
    cbAxis->SetTitleSize(0.125);
    cbAxis->SetTitleFont(42);
    cbAxis->SetLabelSize(0.125);
    cbAxis->SetLabelFont(42);
    cbAxis->Draw();
    c->Update();

    gStyle->SetOptStat(1111);
    c->Update();
}


///////////////////////////////////////////////////////////////
// Main macro function
///////////////////////////////////////////////////////////////


// Help()
// Reprints the available function summary. Call at any time
// from the ROOT prompt.
///////////////////////////////////////////////////////////////

void Help() {
    printf("  Available functions:\n");
    printf("\n");
    printf("    ListShowerEvents(N, M=1, energyCut=0.0)\n");
    printf("      List N shower events starting from event M (1-based).\n");
    printf("      Prints eventNum, eventIdx, and gate type:\n");
    printf("        Coincidence : both arms >= energyCut\n");
    printf("        Left        : left arm only >= energyCut\n");
    printf("        Right       : right arm only >= energyCut\n");
    printf("        None        : both arms below energyCut\n");
    printf("\n");
    printf("    PlotEventShower(eventNum=1)\n");
    printf("      3D visualization of shower depth slices for one event.\n");
    printf("      Rings drawn at 90%% energy containment radius per slice.\n");
    printf("      Color encodes deposited energy (blue=low, red=high).\n");
    printf("\n");
    printf("    PlotShowerParams(eventNum=1)\n");
    printf("      Four-panel plot of Grindhammer-Peters shower parameters:\n");
    printf("        Panel 1: R_C vs tau  (core radius)\n");
    printf("        Panel 2: R_T vs tau  (tail radius)\n");
    printf("        Panel 3: p   vs tau  (core weight)\n");
    printf("        Panel 4: 1/E dE/dt vs t (longitudinal profile,\n");
    printf("                 data points + analytic gamma curve)\n");
    printf("\n");
    printf("    PlotShowerFace(eventNum=1)\n");
    printf("      2D energy density map on the detector face (X x Y).\n");
    printf("      Colorbar on right, PMT segment energy grid top-right.\n");
    printf("\n");
    printf("  eventNum is 1-based in all functions.\n");
    printf("\n");
    printf("  Examples:\n");
    printf("    ListShowerEvents(20)\n");
    printf("    ListShowerEvents(20, 5, 0.5)\n");
    printf("    PlotEventShower(1)\n");
    printf("    PlotShowerParams(3)\n");
    printf("    PlotShowerFace(1)\n");
    printf("=============================================================\n");
    printf("\n");
}

void MolPolDetectorResponseShowerInspection(const char *filename) {

    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        printf("ERROR: Could not open file: %s\n", filename);
        return;
    }

    gTDetResp = (TTree *)f->Get("TDetResp");
    if (!gTDetResp) {
        printf("ERROR: TDetResp tree not found in %s\n", filename);
        printf("       This does not appear to be a valid detresponse output file.\n");
        f->Close();
        return;
    }

    // Wire up TDetResp branches using the standard reader
    SetupDetRespBranches(gTDetResp);

    gTShower = (TTree *)f->Get("TShower");
    if (!gTShower) {
        printf("ERROR: TShower tree not found in %s\n", filename);
        printf("       Re-run MolPolDetectorResponse.C with writeShowerDiag = true.\n");
        f->Close();
        return;
    }

    // Build event index at load time for fast lookup.
    BuildEventIndex();

    Int_t nShowerEvents = (Int_t)gShowerEventList.size();
    Long64_t nDetResp   = gTDetResp->GetEntries();

    printf("\n");
    printf("=============================================================\n");
    printf("  MolPolDetectorResponseShowerInspection\n");
    printf("=============================================================\n");
    printf("  File        : %s\n", filename);
    printf("  TDetResp    : %lld entries\n", nDetResp);
    printf("  TShower     : %d unique events indexed\n", nShowerEvents);
    printf("\n");
    Help();
}

#endif // MOLPOL_DETECTOR_RESPONSE_SHOWER_INSPECTION_C