///////////////////////////////////////////////////////////////
// MolPolDetectorResponse.C
// ROOT macro for modeling the detector response of the MolPol
// lead-fiber calorimeter from flux-plane (detector 9) hit data.
//
// Models electromagnetic shower development in Pb using the
// Grindhammer-Peters parameterization and distributes deposited
// energy across the 2x4 PMT segment grid. For each event,
// produces DR_NTHROWS (20) MC samples of the PMT response
// with independent fiber penetration depths and energy smearing.
//
// Shower parameterization reference:
//   G. Grindhammer and S. Peters,
//   "The Parameterized Simulation of Electromagnetic Showers
//    in Homogeneous and Sampling Calorimeters",
//   arXiv:hep-ex/0001020 (2000).
//
// Energy resolution reference:
//   R. Livan, M. Vercesi, R. Wigmans,
//   "Scintillating-Fibre Calorimetry" (1995).
//
// Usage:
//   root 'MolPolDetectorResponse.C("input.root")'
//   root 'MolPolDetectorResponse.C("input.root", true)'           // with shower diagnostics
//   root 'MolPolDetectorResponse.C("input.root", false, 12345)'   // custom RNG seed
//   root 'MolPolDetectorResponse.C("filelist.txt")'
//
///////////////////////////////////////////////////////////////

#ifndef MOLPOL_DETECTOR_RESPONSE_C
#define MOLPOL_DETECTOR_RESPONSE_C

#include "MolPolData.h"
#include "MolPolAnalysis.h"
#include "MolPolDetectorResponse.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"

#include <iostream>
#include <iomanip>

// Maximum number of events for which shower diagnostics are stored
const Int_t MAX_DIAG_EVENTS = 10000;

///////////////////////////////////////////////////////////////
// Main macro function
///////////////////////////////////////////////////////////////

void MolPolDetectorResponse(const char *fileList,
                             Bool_t writeShowerDiag = false,
                             UInt_t seed = 1800707365) {

    TChain *T = SetupMolpolChain(fileList);
    if (!T) return;

    SetupMolpolBranches(T);

    // Print material constants
    printf("=== MolPolDetectorResponse ===\n");
    printf("Lead (Pb) shower parameters:\n");
    printf("  Z   = %.1f\n", PbConst::Z);
    printf("  A   = %.1f\n", PbConst::A);
    printf("  rho = %.2f g/cm^3\n", PbConst::rho);
    printf("  X0  = %.4f cm (%.4f g/cm^2)\n", PbConst::X0, PbConst::X0g);
    printf("  Ec  = %.4f MeV (Eq. 5)\n", PbConst::Ec);
    printf("  RM  = %.4f cm\n", PbConst::RM);
    printf("  Es  = %.1f MeV\n", PbConst::Es);
    printf("Calorimeter:\n");
    printf("  Depth    = %.1f cm (%.1f X0)\n",
           DetGeom::caloDepth, DetGeom::caloDepthX0);
    printf("  Resolution: sigma/E = %.1f%%/sqrt(E) + %.1f%%\n",
           CaloResolution::a * 100.0, CaloResolution::b * 100.0);
    printf("MC throws per event: %d\n", DR_NTHROWS);
    printf("Fiber model:\n");
    printf("  Fiber X0:        %.1f cm\n", DetGeom::fiberX0);
    printf("  Fiber diameter:  %.1f mm\n", DetGeom::fiberDiam * 10.0);
    printf("  Fiber vol frac:  %.2f\n", DetGeom::fiberVolFrac);
    printf("  Detector angle:  %.1f deg\n", DetGeom::detAngleDeg);
    printf("Shower diagnostics: %s\n", writeShowerDiag ? "ON" : "OFF");
    if (writeShowerDiag)
        printf("  Max diagnostic events: %d\n", MAX_DIAG_EVENTS);

    // Random number generator (seeded for reproducibility)
    TRandom3 *rng = new TRandom3(seed);
    printf("RNG seed: %u\n", seed);

    // Build output filename
    TString outName(fileList);
    if (outName.EndsWith(".root"))
        outName.ReplaceAll(".root", "_detresponse.root");
    else
        outName = "detresponse_output.root";

    printf("Output file: %s\n", outName.Data());

    TFile *fOut = new TFile(outName.Data(), "RECREATE");

    // --- Output tree: per-event PMT response ---
    TTree *tOut = new TTree("TDetResp", "Detector response per event");

    // Event-level physics branches (copied through)
    Double_t out_evXs, out_evAsym, out_evUnpolWght;
    Double_t out_evPolPlusWghtZ;
    Double_t out_evPolMinusWghtZ;

    tOut->Branch("evXs",            &out_evXs,            "evXs/D");
    tOut->Branch("evAsym",          &out_evAsym,          "evAsym/D");
    tOut->Branch("evUnpolWght",     &out_evUnpolWght,     "evUnpolWght/D");
    tOut->Branch("evPolPlusWghtZ",  &out_evPolPlusWghtZ,  "evPolPlusWghtZ/D");
    tOut->Branch("evPolMinusWghtZ", &out_evPolMinusWghtZ, "evPolMinusWghtZ/D");

    // PMT energy deposits: DR_NTHROWS throws x 8 segments (Float_t)
    const Int_t nPmtVals = DR_NTHROWS * 8;
    Float_t  out_pmtE[nPmtVals];          // [throw*8 + seg]
    Float_t  out_pmtETotal[DR_NTHROWS];
    Float_t  out_pmtELeft[DR_NTHROWS];
    Float_t  out_pmtERight[DR_NTHROWS];
    Float_t  out_hitZ[DR_NTHROWS];        // shower start depth per throw [cm]
    Int_t    out_nHitsDet9;

    tOut->Branch("pmtE",       out_pmtE,
        TString::Format("pmtE[%d]/F", nPmtVals).Data());
    tOut->Branch("pmtETotal",  out_pmtETotal,
        TString::Format("pmtETotal[%d]/F", DR_NTHROWS).Data());
    tOut->Branch("pmtELeft",   out_pmtELeft,
        TString::Format("pmtELeft[%d]/F", DR_NTHROWS).Data());
    tOut->Branch("pmtERight",  out_pmtERight,
        TString::Format("pmtERight[%d]/F", DR_NTHROWS).Data());
    tOut->Branch("hitZ",       out_hitZ,
        TString::Format("hitZ[%d]/F", DR_NTHROWS).Data());
    tOut->Branch("nHitsDet9",  &out_nHitsDet9,   "nHitsDet9/I");

    // Per-hit momentum vectors (variable-length, indexed by nHitsDet9)
    Double_t out_hitPx[DR_MAXNHIT];
    Double_t out_hitPy[DR_MAXNHIT];
    Double_t out_hitPz[DR_MAXNHIT];
    Double_t out_hitE[DR_MAXNHIT];

    tOut->Branch("hitPx", out_hitPx, "hitPx[nHitsDet9]/D");
    tOut->Branch("hitPy", out_hitPy, "hitPy[nHitsDet9]/D");
    tOut->Branch("hitPz", out_hitPz, "hitPz[nHitsDet9]/D");
    tOut->Branch("hitE",  out_hitE,  "hitE[nHitsDet9]/D");

    // --- Optional shower diagnostic tree ---
    TTree *tDiag = 0;
    Long64_t nDiagEventsFilled = 0;
    if (writeShowerDiag) {
        tDiag = new TTree("TShower", "Shower profile diagnostics (per hit)");
        tDiag->Branch("eventIdx",    &diag_eventIdx,    "eventIdx/I");
        tDiag->Branch("hitIdx",      &diag_hitIdx,      "hitIdx/I");
        tDiag->Branch("hitX",        &diag_hitX,        "hitX/F");
        tDiag->Branch("hitY",        &diag_hitY,        "hitY/F");
        tDiag->Branch("hitE",        &diag_hitE,        "hitE/F");
        tDiag->Branch("nDepthSteps", &diag_nDepthSteps, "nDepthSteps/I");

        TString leafDepthT  = TString::Format("depthT[%d]/F",  MAX_DEPTH_STEPS);
        TString leafSliceE  = TString::Format("sliceE[%d]/F",  MAX_DEPTH_STEPS);
        TString leafRC      = TString::Format("RC[%d]/F",      MAX_DEPTH_STEPS);
        TString leafRT      = TString::Format("RT[%d]/F",      MAX_DEPTH_STEPS);
        TString leafP       = TString::Format("p[%d]/F",       MAX_DEPTH_STEPS);
        TString leafPmtFrac = TString::Format("pmtFrac[%d]/F", MAX_DEPTH_STEPS * 8);

        tDiag->Branch("depthT",  diag_depthT,  leafDepthT.Data());
        tDiag->Branch("sliceE",  diag_sliceE,  leafSliceE.Data());
        tDiag->Branch("RC",      diag_RC,      leafRC.Data());
        tDiag->Branch("RT",      diag_RT,      leafRT.Data());
        tDiag->Branch("p",       diag_p,       leafP.Data());
        tDiag->Branch("pmtFrac", diag_pmtFrac, leafPmtFrac.Data());

        TString leafHitZ = TString::Format("hitZ[%d]/F", DR_NTHROWS);
        tDiag->Branch("hitZ", diag_hitZ, leafHitZ.Data());
    }

    // --- First pass: determine beam energy ---
    Long64_t nEntries = T->GetEntries();
    printf("Entries: %lld\n", nEntries);

    Double_t maxSumEvP = 0.0;
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        T->GetEntry(iEntry);
        if (hitFluxTrackCoinc()) {
            Double_t sumEvP = evP[0] + evP[1];
            if (sumEvP > maxSumEvP) maxSumEvP = sumEvP;
        }
    }
    printf("Beam energy (maxSumEvP): %.3f GeV\n", maxSumEvP);

    Long64_t nHitsProcessed = 0;
    Long64_t nHitsSkippedNonEM = 0;
    Long64_t nEnergyViolations = 0;
    Long64_t nOverBeamEnergy = 0;
    Long64_t nCompletePassThrough = 0;
    Double_t totalEnergyIn = 0.0;       // sum of all hit energies entering det 9
    Double_t totalEnergyDeposited = 0.0; // sum of all pmtETotal (throw 0, pre-smear)

    // Diagnostic histograms for fiber penetration model
    TH1F *hZstart = new TH1F("hZstart",
        "Shower Start Depth (z_{start});z_{start} [cm];Throws",
        500, 0.0, 50.0);
    TH1F *hZstartLead = new TH1F("hZstartLead",
        "z_{start} Lead;z_{start} [cm];Throws",
        500, 0.0, 50.0);
    TH1F *hZstartFiber = new TH1F("hZstartFiber",
        "z_{start} Fiber;z_{start} [cm];Throws",
        500, 0.0, 50.0);
    hZstart->SetLineColor(kBlack);
    hZstart->SetLineWidth(2);
    hZstartLead->SetLineColor(kRed);
    hZstartLead->SetLineWidth(2);
    hZstartFiber->SetLineColor(kBlue);
    hZstartFiber->SetLineWidth(2);
    TH1F *hZavailable = new TH1F("hZavailable",
        "Available Fiber Depth (z_{available});z_{available} [cm];Throws",
        500, 0.0, 50.0);
    TH2F *hZstartVsZavail = new TH2F("hZstartVsZavail",
        "z_{start} vs z_{available} (fiber hits);z_{available} [cm];z_{start} [cm]",
        200, 0.0, 50.0, 200, 0.0, 50.0);
    TH1F *hDtrans = new TH1F("hDtrans",
        "Fiber Transverse Distance (d_{trans});d_{trans} [cm];Throws",
        200, 0.0, DetGeom::fiberDiam);

    // Theta_incident histograms by particle type
    TH1F *hThetaIncAll = new TH1F("hThetaIncAll",
        "#theta_{incident};#theta_{incident} [deg];Hits", 500, 0.0, 25.0);
    TH1F *hThetaIncMoller = new TH1F("hThetaIncMoller",
        "#theta_{inc} Gen Mollers;#theta_{incident} [deg];Hits", 500, 0.0, 25.0);
    TH1F *hThetaIncOtherE = new TH1F("hThetaIncOtherE",
        "#theta_{inc} Other e-;#theta_{incident} [deg];Hits", 500, 0.0, 25.0);
    TH1F *hThetaIncGamma = new TH1F("hThetaIncGamma",
        "#theta_{inc} #gamma;#theta_{incident} [deg];Hits", 500, 0.0, 25.0);
    TH1F *hThetaIncProton = new TH1F("hThetaIncProton",
        "#theta_{inc} proton;#theta_{incident} [deg];Hits", 500, 0.0, 25.0);
    hThetaIncAll->SetLineColor(kBlack);       hThetaIncAll->SetLineWidth(2);
    hThetaIncMoller->SetLineColor(kRed);      hThetaIncMoller->SetLineWidth(2);
    hThetaIncOtherE->SetLineColor(kBlue);     hThetaIncOtherE->SetLineWidth(2);
    hThetaIncGamma->SetLineColor(kGreen+2);   hThetaIncGamma->SetLineWidth(2);
    hThetaIncProton->SetLineColor(kMagenta);  hThetaIncProton->SetLineWidth(2);

    // Energy diagnostic histograms (2x2 canvas)
    TH2F *hEdetVsEactual = new TH2F("hEdetVsEactual",
        "E_{detected} vs E_{actual};E_{actual} [GeV];E_{detected} [GeV]",
        200, 0.0, 12.0, 200, 0.0, 12.0);
    TH1F *hLeakageLeft = new TH1F("hLeakageLeft",
        "Shower Leakage per Arm;Leakage [GeV];Events",
        200, -0.5, 2.0);
    hLeakageLeft->SetLineColor(kRed);
    hLeakageLeft->SetLineWidth(2);
    TH1F *hLeakageRight = new TH1F("hLeakageRight",
        "Shower Leakage per Arm;Leakage [GeV];Events",
        200, -0.5, 2.0);
    hLeakageRight->SetLineColor(kBlue);
    hLeakageRight->SetLineWidth(2);
    TH2F *hLeakageVsTheta = new TH2F("hLeakageVsTheta",
        "Shower Leakage vs #theta_{incident};#theta_{incident} [deg];Leakage [GeV]",
        200, 0.0, 25.0, 200, -0.5, 3.0);
    TH1F *hEratio = new TH1F("hEratio",
        "Energy Ratio (E_{det}/E_{actual});E_{det}/E_{actual};Events",
        200, 0.0, 1.2);

    // Scatter plot data for panel 6: shower vs pass-through
    // (fiber hits only, shower points capped to avoid memory issues)
    const Int_t MAX_SCATTER_PTS = 200000;
    std::vector<Double_t> scatShowerX, scatShowerY;   // z_avail, z_start_raw
    std::vector<Double_t> scatPassThruX, scatPassThruY;

    // --- Second pass: shower processing ---
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        T->GetEntry(iEntry);

        if (iEntry % 10000 == 0 && iEntry > 0)
            printf("  Entry %lld / %lld (%.1f%%)\n",
                   iEntry, nEntries, 100.0 * iEntry / nEntries);

        // Copy event-level branches
        out_evXs            = evXs;
        out_evAsym          = evAsym;
        out_evUnpolWght     = evUnpolWght;
        out_evPolPlusWghtZ  = evPolPlusWghtZ;
        out_evPolMinusWghtZ = evPolMinusWghtZ;

        // Zero all output arrays
        for (Int_t v = 0; v < nPmtVals; v++) out_pmtE[v] = 0.0f;
        for (Int_t k = 0; k < DR_NTHROWS; k++) {
            out_pmtETotal[k] = 0.0f;
            out_pmtELeft[k]  = 0.0f;
            out_pmtERight[k] = 0.0f;
            out_hitZ[k]      = 0.0f;
        }
        out_nHitsDet9 = 0;

        // Diagnostics
        Bool_t diagThisEvent = (tDiag != 0 && nDiagEventsFilled < MAX_DIAG_EVENTS);

        // --- Stage 1: Collect EM hits, compute shower profiles ---
        // Store profiles for all EM hits (up to DR_MAXNHIT)
        ShowerProfile profiles[DR_MAXNHIT];
        Double_t sumHitE = 0.0;
        Double_t sumHitELeft = 0.0;
        Double_t sumHitERight = 0.0;
        Int_t nProfiles = 0;

        for (Int_t i = 0; i < hitN; i++) {
            if (hitDet[i] != 9) continue;

            // Compute theta_incident for all particles on det 9
            {
                Double_t pT_i = sqrt(hitPx[i]*hitPx[i] + hitPy[i]*hitPy[i]);
                Double_t theta_p = atan2(pT_i, hitPz[i]);
                Double_t theta_inc = fabs(theta_p - DetGeom::detAngleRad);
                Double_t theta_deg = theta_inc * 180.0 / 3.14159265358979;

                hThetaIncAll->Fill(theta_deg);

                Int_t pid  = hitPid[i];
                Int_t trid = hitTrid[i];
                if (trid == 1 || trid == 2)
                    hThetaIncMoller->Fill(theta_deg);
                else if (pid == 11 || pid == -11)
                    hThetaIncOtherE->Fill(theta_deg);
                else if (pid == 22)
                    hThetaIncGamma->Fill(theta_deg);
                else if (pid == 2212)
                    hThetaIncProton->Fill(theta_deg);
            }

            // EM-only filter
            Int_t pid = hitPid[i];
            if (pid != 11 && pid != -11 && pid != 22) {
                nHitsSkippedNonEM++;
                continue;
            }

            if (out_nHitsDet9 >= DR_MAXNHIT) continue;  // safety cap

            Double_t hx_cm = hitLx[i] * 100.0;
            Double_t hy_cm = hitLy[i] * 100.0;

            sumHitE += hitE[i];
            if (hx_cm >= 0.0)
                sumHitELeft += hitE[i];
            else
                sumHitERight += hitE[i];

            // Store hit momentum and energy
            out_hitPx[out_nHitsDet9] = hitPx[i];
            out_hitPy[out_nHitsDet9] = hitPy[i];
            out_hitPz[out_nHitsDet9] = hitPz[i];
            out_hitE[out_nHitsDet9]  = hitE[i];

            // Compute shower profile
            ComputeShowerProfile(hitE[i], hx_cm, hy_cm, profiles[nProfiles]);

            nProfiles++;
            out_nHitsDet9++;
            nHitsProcessed++;
        }

        // Track diagnostic event count
        if (diagThisEvent && out_nHitsDet9 > 0)
            nDiagEventsFilled++;

        // --- Stage 2: MC throws for fiber penetration + accumulation ---
        Bool_t throwHasPassThrough[DR_NTHROWS];
        Bool_t eventHasAnyPassThrough = false;

        // Per-hit-per-throw start depths for TShower diagnostic
        Float_t hitZperHit[DR_MAXNHIT][DR_NTHROWS];
        for (Int_t h = 0; h < nProfiles; h++)
            for (Int_t kk = 0; kk < DR_NTHROWS; kk++)
                hitZperHit[h][kk] = 0.0f;

        for (Int_t k = 0; k < DR_NTHROWS; k++) {
            // Zero pmtE for this throw
            Double_t throwPmtE[8] = {0,0,0,0,0,0,0,0};
            throwHasPassThrough[k] = false;

            // Accumulate all hits with their individual start depths
            Double_t maxStartDepth = 0.0;

            for (Int_t h = 0; h < nProfiles; h++) {
                Double_t z_avail = -1.0;
                Double_t z_start_raw = -1.0;
                Double_t d_trans = -1.0;
                Double_t startDepth = ComputeFiberPenetration(
                    out_hitPx[h], out_hitPy[h], out_hitPz[h], rng,
                    &z_avail, &z_start_raw, &d_trans);

                // Store per-hit start depth for TShower diagnostic
                hitZperHit[h][k] = (Float_t)startDepth;

                // Fill diagnostic histograms (all throws, all hits)
                if (startDepth >= 0.0) {
                    hZstart->Fill(startDepth);
                    if (z_avail >= 0.0)
                        hZstartFiber->Fill(startDepth);
                    else
                        hZstartLead->Fill(startDepth);
                }
                if (z_avail >= 0.0) {
                    hZavailable->Fill(z_avail);
                    if (z_start_raw >= 0.0)
                        hZstartVsZavail->Fill(z_avail, z_start_raw);
                }
                if (d_trans >= 0.0) {
                    hDtrans->Fill(d_trans);
                }

                if (startDepth < 0.0) {
                    throwHasPassThrough[k] = true;
                    eventHasAnyPassThrough = true;
                    // Collect pass-through scatter point (fiber only)
                    if (z_avail >= 0.0) {
                        scatPassThruX.push_back(z_avail);
                        scatPassThruY.push_back(z_start_raw);
                    }
                    continue;
                }

                // Collect shower scatter point (fiber only, capped)
                if (z_avail >= 0.0 &&
                    (Int_t)scatShowerX.size() < MAX_SCATTER_PTS) {
                    scatShowerX.push_back(z_avail);
                    scatShowerY.push_back(z_start_raw);
                }

                AccumulateShowerInBlock(profiles[h], startDepth, throwPmtE);
                if (startDepth > maxStartDepth) maxStartDepth = startDepth;
            }

            // Store start depth for this throw
            out_hitZ[k] = (Float_t)maxStartDepth;

            // Copy to output (before smearing, for conservation check on throw 0)
            for (Int_t s = 0; s < 8; s++)
                out_pmtE[k * 8 + s] = (Float_t)throwPmtE[s];
        }

        // --- Pass-through tracking for this event ---
        if (eventHasAnyPassThrough) {
            nCompletePassThrough++;
        }

        // --- Fill TShower diagnostics (deferred from Stage 1) ---
        // Now that we have the per-throw start depths, fill TShower
        // for each hit with the profile + hitZ data.
        if (diagThisEvent) {
            for (Int_t h = 0; h < nProfiles; h++) {
                FillShowerDiag(profiles[h], tDiag,
                               (Int_t)iEntry, h, hitZperHit[h]);
            }
        }

        // --- Energy conservation check on throw 0 (before smearing) ---
        if (out_nHitsDet9 > 0) {
            Double_t throw0Total = 0.0;
            for (Int_t s = 0; s < 8; s++) throw0Total += out_pmtE[s];

            totalEnergyIn += sumHitE;
            totalEnergyDeposited += throw0Total;

            if (throw0Total > sumHitE * 1.01) {
                nEnergyViolations++;
                if (nEnergyViolations <= 20) {
                    printf("  WARNING: Energy violation at eventIdx %lld: "
                           "sumHitE = %.4f GeV, pmtETotal = %.4f GeV "
                           "(ratio = %.3f, nHits = %d)\n",
                           iEntry, sumHitE, throw0Total,
                           throw0Total / sumHitE, out_nHitsDet9);
                }
                if (nEnergyViolations == 20)
                    printf("  (suppressing further energy violation warnings)\n");
            }

            if (maxSumEvP > 0.0 && throw0Total > maxSumEvP * 1.01) {
                nOverBeamEnergy++;
                if (nOverBeamEnergy <= 20) {
                    printf("  INFO: pmtETotal > beam energy at eventIdx %lld: "
                           "beam = %.3f GeV, sumHitE = %.4f GeV, "
                           "pmtETotal = %.4f GeV, nHits = %d\n",
                           iEntry, maxSumEvP, sumHitE, throw0Total,
                           out_nHitsDet9);
                }
                if (nOverBeamEnergy == 20)
                    printf("  (suppressing further beam energy warnings)\n");
            }
        }

        // --- Energy diagnostic histograms (throw 0, pre-smear) ---
        if (out_nHitsDet9 > 0) {
            Double_t detTotal = 0.0;
            Double_t detLeft = 0.0;
            Double_t detRight = 0.0;
            for (Int_t s = 0; s < 4; s++) detLeft  += out_pmtE[s];
            for (Int_t s = 4; s < 8; s++) detRight += out_pmtE[s];
            detTotal = detLeft + detRight;

            // Panel 1: detected vs actual
            hEdetVsEactual->Fill(sumHitE, detTotal);

            // Panel 2: leakage per arm
            if (sumHitELeft > 0.0)
                hLeakageLeft->Fill(sumHitELeft - detLeft);
            if (sumHitERight > 0.0)
                hLeakageRight->Fill(sumHitERight - detRight);

            // Panel 3: leakage vs theta_incident (per-hit)
            Double_t totalLeakage = sumHitE - detTotal;
            for (Int_t h = 0; h < out_nHitsDet9; h++) {
                Double_t pT_h = sqrt(out_hitPx[h]*out_hitPx[h] +
                                     out_hitPy[h]*out_hitPy[h]);
                Double_t theta_p = atan2(pT_h, out_hitPz[h]);
                Double_t theta_i = fabs(theta_p - DetGeom::detAngleRad);
                Double_t theta_i_deg = theta_i * 180.0 / 3.14159265358979;
                hLeakageVsTheta->Fill(theta_i_deg, totalLeakage);
            }

            // Panel 4: energy ratio
            if (sumHitE > 0.0)
                hEratio->Fill(detTotal / sumHitE);
        }

        // --- Per-PMT energy smearing, independently per throw ---
        for (Int_t k = 0; k < DR_NTHROWS; k++) {
            for (Int_t s = 0; s < 8; s++) {
                Float_t &e = out_pmtE[k * 8 + s];
                if (e > 0.0f) {
                    Double_t sigma = CaloResolution::Sigma((Double_t)e);
                    e = (Float_t)rng->Gaus((Double_t)e, sigma);
                    if (e < 0.0f) e = 0.0f;
                }
            }

            // Compute sums from smeared values
            out_pmtELeft[k]  = 0.0f;
            out_pmtERight[k] = 0.0f;
            for (Int_t s = 0; s < 4; s++) out_pmtELeft[k]  += out_pmtE[k * 8 + s];
            for (Int_t s = 4; s < 8; s++) out_pmtERight[k] += out_pmtE[k * 8 + s];
            out_pmtETotal[k] = out_pmtELeft[k] + out_pmtERight[k];
        }

        tOut->Fill();
    }

    // --- Summary ---
    printf("\n=== Processing Summary ===\n");
    printf("Beam energy (maxSumEvP): %.3f GeV\n", maxSumEvP);
    printf("Events processed:        %lld\n", nEntries);
    printf("Detector hits showered:  %lld\n", nHitsProcessed);
    printf("Non-EM hits skipped:     %lld\n", nHitsSkippedNonEM);
    printf("Energy violations:       %lld (pmtETotal > sumHitE)\n", nEnergyViolations);
    printf("Over beam energy:        %lld (pmtETotal > beam)\n", nOverBeamEnergy);
    printf("Complete pass-throughs:  %lld\n", nCompletePassThrough);
    printf("MC throws per event:     %d\n", DR_NTHROWS);
    if (totalEnergyIn > 0.0) {
        Double_t leakageFrac = 1.0 - totalEnergyDeposited / totalEnergyIn;
        printf("Energy leakage fraction: %.4f (throw 0, pre-smear)\n", leakageFrac);
    }
    if (tDiag)
        printf("Showering events saved:  %lld / %d max\n",
               nDiagEventsFilled, MAX_DIAG_EVENTS);

    // Write and close
    fOut->cd();
    tOut->Write();
    if (tDiag) tDiag->Write();

    // Write diagnostic histograms
    hZstart->Write();
    hZstartLead->Write();
    hZstartFiber->Write();
    hZavailable->Write();
    hZstartVsZavail->Write();
    hDtrans->Write();
    hThetaIncAll->Write();
    hThetaIncMoller->Write();
    hThetaIncOtherE->Write();
    hThetaIncGamma->Write();
    hThetaIncProton->Write();

    // Write energy diagnostic histograms
    hEdetVsEactual->Write();
    hLeakageLeft->Write();
    hLeakageRight->Write();
    hLeakageVsTheta->Write();
    hEratio->Write();

    // Detach histograms from file before closing so they survive for drawing
    hZstart->SetDirectory(0);
    hZstartLead->SetDirectory(0);
    hZstartFiber->SetDirectory(0);
    hZavailable->SetDirectory(0);
    hZstartVsZavail->SetDirectory(0);
    hDtrans->SetDirectory(0);
    hThetaIncAll->SetDirectory(0);
    hThetaIncMoller->SetDirectory(0);
    hThetaIncOtherE->SetDirectory(0);
    hThetaIncGamma->SetDirectory(0);
    hThetaIncProton->SetDirectory(0);
    hEdetVsEactual->SetDirectory(0);
    hLeakageLeft->SetDirectory(0);
    hLeakageRight->SetDirectory(0);
    hLeakageVsTheta->SetDirectory(0);
    hEratio->SetDirectory(0);

    fOut->Close();

    // Draw diagnostic histograms (3x2 layout)
    gStyle->SetPalette(55);

    TCanvas *cFiberDiag = new TCanvas("cFiberDiag",
        "Fiber Penetration Diagnostics", 1500, 900);
    cFiberDiag->Divide(3, 2);

    // Panel 1: z_start (total/lead/fiber)
    cFiberDiag->cd(1);
    gPad->SetLogy();
    hZstart->SetMinimum(0.9);
    hZstart->Draw();
    hZstartLead->Draw("same");
    hZstartFiber->Draw("same");
    gPad->Update();
    TPaveStats *st1 = (TPaveStats*)hZstart->FindObject("stats");
    if (st1) { st1->SetX1NDC(0.675); st1->SetY1NDC(0.15);
               st1->SetX2NDC(0.875); st1->SetY2NDC(0.325); }
    gPad->Modified();
    TLegend *legZ = new TLegend(0.66, 0.66, 0.9, 0.9);
    legZ->AddEntry(hZstart, "Total", "l");
    legZ->AddEntry(hZstartLead, "Lead", "l");
    legZ->AddEntry(hZstartFiber, "Fiber", "l");
    legZ->SetBorderSize(0);
    legZ->SetFillStyle(0);
    legZ->Draw();

    // Panel 2: z_available
    cFiberDiag->cd(2);
    gPad->SetLogy();
    hZavailable->SetMinimum(0.9);
    hZavailable->Draw();
    gPad->Update();
    TPaveStats *st2 = (TPaveStats*)hZavailable->FindObject("stats");
    if (st2) { st2->SetX1NDC(0.675); st2->SetY1NDC(0.15);
               st2->SetX2NDC(0.875); st2->SetY2NDC(0.325); }
    gPad->Modified();

    // Panel 3: z_start_raw vs z_available (2D COLZ with L-corner boundary)
    cFiberDiag->cd(3);
    gPad->SetRightMargin(0.15);
    hZstartVsZavail->SetTitle("z_{start,raw}^{MC} vs z_{avail}");
    hZstartVsZavail->SetStats(kFALSE);
    hZstartVsZavail->Draw("COLZ");
    Double_t diagMax = std::max(hZstartVsZavail->GetXaxis()->GetXmax(),
                                hZstartVsZavail->GetYaxis()->GetXmax());
    // Pass-through boundary: L-corner at (30,30)
    TLine *lPTH = new TLine(DetGeom::caloDepth, DetGeom::caloDepth,
                             diagMax, DetGeom::caloDepth);
    lPTH->SetLineColor(kMagenta);
    lPTH->SetLineWidth(2);
    lPTH->SetLineStyle(2);
    lPTH->Draw("same");
    TLine *lPTV = new TLine(DetGeom::caloDepth, DetGeom::caloDepth,
                             DetGeom::caloDepth, diagMax);
    lPTV->SetLineColor(kMagenta);
    lPTV->SetLineWidth(2);
    lPTV->SetLineStyle(2);
    lPTV->Draw("same");

    // Panel 4: d_trans
    cFiberDiag->cd(4);
    hDtrans->Draw();
    gPad->Update();
    TPaveStats *st4 = (TPaveStats*)hDtrans->FindObject("stats");
    if (st4) { st4->SetX1NDC(0.15); st4->SetY1NDC(0.7);
               st4->SetX2NDC(0.35); st4->SetY2NDC(0.875); }
    gPad->Modified();

    // Panel 5: theta_incident by particle type
    cFiberDiag->cd(5);
    gPad->SetLogy();
    hThetaIncAll->SetMinimum(0.9);
    hThetaIncAll->SetTitle("#theta_{incident} by Particle Type");
    hThetaIncAll->SetStats(kFALSE);
    hThetaIncAll->Draw();
    hThetaIncMoller->Draw("same");
    hThetaIncOtherE->Draw("same");
    hThetaIncGamma->Draw("same");
    hThetaIncProton->Draw("same");
    TLegend *legT = new TLegend(0.55, 0.55, 0.88, 0.88);
    legT->AddEntry(hThetaIncAll, "All", "l");
    legT->AddEntry(hThetaIncMoller, "Gen Mollers", "l");
    legT->AddEntry(hThetaIncOtherE, "Other e-", "l");
    legT->AddEntry(hThetaIncGamma, "#gamma", "l");
    legT->AddEntry(hThetaIncProton, "proton", "l");
    legT->SetBorderSize(0);
    legT->SetFillStyle(0);
    legT->Draw();

    // Panel 6: Scatter plot — shower (black) vs pass-through (red)
    cFiberDiag->cd(6);
    TH1F *hScatFrame = gPad->DrawFrame(0.0, 0.0, 50.0, 50.0);
    hScatFrame->SetTitle(
        "Shower vs Pass-through (fiber);z_{available} [cm];z_{start,raw}^{MC} [cm]");

    TGraph *gShower = 0;
    TGraph *gPassThru = 0;

    // Shower points (black)
    if (!scatShowerX.empty()) {
        gShower = new TGraph((Int_t)scatShowerX.size(),
                              &scatShowerX[0], &scatShowerY[0]);
        gShower->SetMarkerStyle(6);
        gShower->SetMarkerColor(kBlack);
        gShower->Draw("P same");
    }

    // Pass-through points (red)
    if (!scatPassThruX.empty()) {
        gPassThru = new TGraph((Int_t)scatPassThruX.size(),
                                &scatPassThruX[0], &scatPassThruY[0]);
        gPassThru->SetMarkerStyle(6);
        gPassThru->SetMarkerColor(kRed);
        gPassThru->Draw("P same");
    }

    // Legend
    TLegend *legS = new TLegend(0.15, 0.75, 0.45, 0.88);
    if (gShower)   legS->AddEntry(gShower, "Showered", "p");
    if (gPassThru) legS->AddEntry(gPassThru, "Pass-through", "p");
    legS->SetBorderSize(1);
    legS->SetFillStyle(1001);
    legS->SetFillColor(kWhite);
    legS->Draw();

    cFiberDiag->Update();

    // --- Energy diagnostic canvas (2x2) ---
    TCanvas *cEnergyDiag = new TCanvas("cEnergyDiag",
        "Energy Diagnostics", 1200, 1000);
    cEnergyDiag->Divide(2, 2);

    // Panel 1: E_detected vs E_actual
    cEnergyDiag->cd(1);
    gPad->SetRightMargin(0.15);
    hEdetVsEactual->SetStats(kFALSE);
    hEdetVsEactual->Draw("COLZ");
    // Diagonal: E_det = E_actual
    Double_t eMax = std::max(hEdetVsEactual->GetXaxis()->GetXmax(),
                              hEdetVsEactual->GetYaxis()->GetXmax());
    TLine *lEdiag = new TLine(0.0, 0.0, eMax, eMax);
    lEdiag->SetLineColor(kBlack);
    lEdiag->SetLineWidth(2);
    lEdiag->SetLineStyle(2);
    lEdiag->Draw("same");

    // Panel 2: Leakage per arm
    cEnergyDiag->cd(2);
    Double_t yMaxLeak = std::max(hLeakageLeft->GetMaximum(),
                                  hLeakageRight->GetMaximum()) * 1.15;
    hLeakageLeft->SetMaximum(yMaxLeak);
    hLeakageLeft->Draw("HISTE");
    hLeakageRight->Draw("HISTE same");
    TLegend *legLeak = new TLegend(0.55, 0.7, 0.88, 0.88);
    legLeak->AddEntry(hLeakageLeft, "Left arm", "l");
    legLeak->AddEntry(hLeakageRight, "Right arm", "l");
    legLeak->SetBorderSize(0);
    legLeak->SetFillStyle(0);
    legLeak->Draw();

    // Panel 3: Leakage vs theta_incident
    cEnergyDiag->cd(3);
    gPad->SetRightMargin(0.15);
    hLeakageVsTheta->SetStats(kFALSE);
    // Auto-range x-axis to last populated bin
    Int_t lastBinX = 1;
    for (Int_t bx = hLeakageVsTheta->GetNbinsX(); bx >= 1; bx--) {
        Bool_t hasContent = false;
        for (Int_t by = 1; by <= hLeakageVsTheta->GetNbinsY(); by++) {
            if (hLeakageVsTheta->GetBinContent(bx, by) > 0.0) {
                hasContent = true;
                break;
            }
        }
        if (hasContent) { lastBinX = bx; break; }
    }
    Double_t xMaxTheta = hLeakageVsTheta->GetXaxis()->GetBinUpEdge(lastBinX + 2);
    hLeakageVsTheta->GetXaxis()->SetRangeUser(0.0, xMaxTheta);
    hLeakageVsTheta->Draw("COLZ");

    // Panel 4: Energy ratio
    cEnergyDiag->cd(4);
    hEratio->Draw("HISTE");

    cEnergyDiag->Update();

    printf("Output written to: %s\n", outName.Data());

    delete T;
    delete rng;
}

#endif // MOLPOL_DETECTOR_RESPONSE_C