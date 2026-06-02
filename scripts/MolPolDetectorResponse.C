///////////////////////////////////////////////////////////////
// MolPolDetectorResponse.C
// Eric King 06/02/2026 (updated)
//
// ROOT macro for modeling the detector response of the MolPol
// lead-fiber calorimeter from flux-plane (detector 9) hit data.
//
// Generates the Detector Response and Shower Diagnostics trees
// in a _detresponse.root_ from the MolPol simulation data ROOT file.
//
// Models electromagnetic shower development in Pb using the
// Grindhammer-Peters parameterization and distributes deposited
// energy across the 2x4 PMT segment grid. Produces an output
// ROOT file with per-event PMT energy deposits suitable for
// downstream asymmetry and rate analysis.
//
// Shower parameterization reference:
//   G. Grindhammer and S. Peters,
//   "The Parameterized Simulation of Electromagnetic Showers
//    in Homogeneous and Sampling Calorimeters",
//   arXiv:hep-ex/0001020 (2000).
//   https://arxiv.org/abs/hep-ex/0001020
//
// TODO (potential future improvements):
//   - Momentum-vector-tilted shower axis (non-normal incidence)
//     This is likely the source of the tails seen in polarimeter ADC spectra.
//   - Sampling calorimeter corrections (Appendix A.2 of
//     Grindhammer-Peters) for effective lead-fiber Moliere
//     radius, replacing the current pure-Pb approximation
//
// Usage:
//   root 'MolPolDetectorResponse.C("input.root")'
//   root 'MolPolDetectorResponse.C("input.root", true)'   // with shower diagnostics
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

#include <iostream>
#include <iomanip>

// Maximum number of events for which shower diagnostics are stored
const Int_t MAX_DIAG_EVENTS = 10000;

///////////////////////////////////////////////////////////////
// Main macro function
///////////////////////////////////////////////////////////////

void MolPolDetectorResponse(const char *fileList,
                             Bool_t writeShowerDiag = kFALSE) {

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
    printf("Shower diagnostics: %s\n", writeShowerDiag ? "ON" : "OFF");
    if (writeShowerDiag)
        printf("  Max diagnostic events: %d\n", MAX_DIAG_EVENTS);

    // Random number generator for energy smearing
    TRandom3 *rng = new TRandom3();

    // Setup input chain
    TChain *T = SetupMolpolChain(fileList);
    if (!T) return;
    SetupMolpolBranches(T);

    // Build output filename from first file in chain
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
    Double_t out_evPolPlusWghtX, out_evPolPlusWghtY, out_evPolPlusWghtZ;
    Double_t out_evPolMinusWghtX, out_evPolMinusWghtY, out_evPolMinusWghtZ;

    tOut->Branch("evXs",            &out_evXs,            "evXs/D");
    tOut->Branch("evAsym",          &out_evAsym,          "evAsym/D");
    tOut->Branch("evUnpolWght",     &out_evUnpolWght,     "evUnpolWght/D");
    tOut->Branch("evPolPlusWghtX",  &out_evPolPlusWghtX,  "evPolPlusWghtX/D");
    tOut->Branch("evPolPlusWghtY",  &out_evPolPlusWghtY,  "evPolPlusWghtY/D");
    tOut->Branch("evPolPlusWghtZ",  &out_evPolPlusWghtZ,  "evPolPlusWghtZ/D");
    tOut->Branch("evPolMinusWghtX", &out_evPolMinusWghtX, "evPolMinusWghtX/D");
    tOut->Branch("evPolMinusWghtY", &out_evPolMinusWghtY, "evPolMinusWghtY/D");
    tOut->Branch("evPolMinusWghtZ", &out_evPolMinusWghtZ, "evPolMinusWghtZ/D");

    // PMT energy deposits
    Double_t out_pmtE[DetGeom::NPMT];
    Double_t out_pmtETotal;
    Double_t out_pmtELeft;
    Double_t out_pmtERight;
    Int_t    out_nHitsDet9;

    tOut->Branch("pmtE",       out_pmtE,        "pmtE[8]/D");
    tOut->Branch("pmtETotal",  &out_pmtETotal,   "pmtETotal/D");
    tOut->Branch("pmtELeft",   &out_pmtELeft,    "pmtELeft/D");
    tOut->Branch("pmtERight",  &out_pmtERight,   "pmtERight/D");
    tOut->Branch("nHitsDet9",  &out_nHitsDet9,   "nHitsDet9/I");

    // --- Optional shower diagnostic tree (Float_t to save space) ---
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
    }

    // --- Event loop ---
    Long64_t nEntries = T->GetEntries();
    printf("Processing %lld entries\n", nEntries);

    // --- First pass: determine beam energy ---
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

    // --- Second pass: shower processing ---
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        T->GetEntry(iEntry);

        // Progress reporting
        if (iEntry % 10000 == 0 && iEntry > 0)
            printf("  Entry %lld / %lld (%.1f%%)\n",
                   iEntry, nEntries, 100.0 * iEntry / nEntries);

        // Copy event-level branches
        out_evXs            = evXs;
        out_evAsym          = evAsym;
        out_evUnpolWght     = evUnpolWght;
        out_evPolPlusWghtX  = evPolPlusWghtX;
        out_evPolPlusWghtY  = evPolPlusWghtY;
        out_evPolPlusWghtZ  = evPolPlusWghtZ;
        out_evPolMinusWghtX = evPolMinusWghtX;
        out_evPolMinusWghtY = evPolMinusWghtY;
        out_evPolMinusWghtZ = evPolMinusWghtZ;

        // Zero PMT accumulator
        for (Int_t s = 0; s < DetGeom::NPMT; s++) out_pmtE[s] = 0.0;
        out_pmtETotal = 0.0;
        out_pmtELeft  = 0.0;
        out_pmtERight = 0.0;
        out_nHitsDet9 = 0;

        // Track sum of input hit energies for conservation check
        Double_t sumHitE = 0.0;
        // Store hit info for debug output (up to 10 per event)
        const Int_t MAX_DEBUG_HITS = 10;
        Double_t dbgHitX[MAX_DEBUG_HITS], dbgHitY[MAX_DEBUG_HITS], dbgHitE[MAX_DEBUG_HITS];
        Int_t    dbgHitPid[MAX_DEBUG_HITS], dbgHitTrid[MAX_DEBUG_HITS];
        Int_t nDbgHits = 0;

        // Determine whether this event gets diagnostics
        Bool_t diagThisEvent = (tDiag != 0 && nDiagEventsFilled < MAX_DIAG_EVENTS);
        TTree *diagTreeForEvent = diagThisEvent ? tDiag : 0;

        // Loop over hits on detector 9
        Int_t hitIdxDet9 = 0;
        for (Int_t i = 0; i < hitN; i++) {
            if (hitDet[i] != 9) continue;

            // Accept only electromagnetic particles: e- (11), e+ (-11), gamma (22).
            // All other particles (protons, pions, kaons, etc.) are excluded because:
            //   1. The Grindhammer-Peters parameterization models electromagnetic
            //      showers only; hadronic showers have different profiles.
            //   2. Using hitE for massive particles includes rest mass energy,
            //      which inflates the deposited energy beyond what the beam
            //      contributed as kinetic energy. For example, a proton with
            //      0.5 GeV kinetic energy has hitE = 1.44 GeV due to its
            //      0.938 GeV rest mass.
            Int_t pid = hitPid[i];
            if (pid != 11 && pid != -11 && pid != 22) {
                nHitsSkippedNonEM++;
                continue;
            }

            out_nHitsDet9++;

            // Convert hit position from meters to cm
            Double_t hx_cm = hitLx[i] * 100.0;
            Double_t hy_cm = hitLy[i] * 100.0;

            sumHitE += hitE[i];
            if (nDbgHits < MAX_DEBUG_HITS) {
                dbgHitX[nDbgHits]    = hx_cm;
                dbgHitY[nDbgHits]    = hy_cm;
                dbgHitE[nDbgHits]    = hitE[i];
                dbgHitPid[nDbgHits]  = hitPid[i];
                dbgHitTrid[nDbgHits] = hitTrid[i];
                nDbgHits++;
            }

            // Process shower for this hit
            ProcessHitShower(hitE[i], hx_cm, hy_cm, out_pmtE,
                             diagTreeForEvent, (Int_t)iEntry, hitIdxDet9);

            hitIdxDet9++;
            nHitsProcessed++;
        }

        // Track diagnostic event count
        if (diagThisEvent && hitIdxDet9 > 0)
            nDiagEventsFilled++;

        // Compute left/right/total PMT energy sums
        // Left = L1(0) + L2(1) + L3(2) + L4(3)
        // Right = R1(4) + R2(5) + R3(6) + R4(7)
        for (Int_t s = 0; s < 4; s++) out_pmtELeft  += out_pmtE[s];
        for (Int_t s = 4; s < 8; s++) out_pmtERight += out_pmtE[s];
        out_pmtETotal = out_pmtELeft + out_pmtERight;

        // Energy conservation check: pmtETotal vs sumHitE
        if (out_nHitsDet9 > 0 && out_pmtETotal > sumHitE * 1.01) {
            nEnergyViolations++;
            if (nEnergyViolations <= 20) {
                printf("  WARNING: Energy violation at eventIdx %lld: "
                       "sumHitE = %.4f GeV, pmtETotal = %.4f GeV "
                       "(ratio = %.3f, nHits = %d)\n",
                       iEntry, sumHitE, out_pmtETotal,
                       out_pmtETotal / sumHitE, out_nHitsDet9);
                printf("    pmtE: L1=%.4f L2=%.4f L3=%.4f L4=%.4f "
                       "R1=%.4f R2=%.4f R3=%.4f R4=%.4f\n",
                       out_pmtE[0], out_pmtE[1], out_pmtE[2], out_pmtE[3],
                       out_pmtE[4], out_pmtE[5], out_pmtE[6], out_pmtE[7]);
                for (Int_t ih = 0; ih < nDbgHits; ih++) {
                    Int_t col = (dbgHitX[ih] >= 0.0) ? 0 : 1;
                    Int_t row = 3;
                    for (Int_t r = 0; r < 4; r++) {
                        if (dbgHitY[ih] >= DetGeom::rowYBound[r + 1]) {
                            row = r; break;
                        }
                    }
                    Int_t seg = col * 4 + row;
                    printf("    hit[%d]: x=%.2f cm, y=%.2f cm, E=%.4f GeV, "
                           "pid=%d, trid=%d -> %s (seg %d)\n",
                           ih, dbgHitX[ih], dbgHitY[ih], dbgHitE[ih],
                           dbgHitPid[ih], dbgHitTrid[ih],
                           DetGeom::PmtSegName[seg], seg);
                }
            }
            if (nEnergyViolations == 20)
                printf("  (suppressing further energy violation warnings)\n");
        }

        // Beam energy check: flag events where sumHitE > beam energy
        if (out_nHitsDet9 > 0 && maxSumEvP > 0.0 &&
            out_pmtETotal > maxSumEvP * 1.01) {
            nOverBeamEnergy++;
            if (nOverBeamEnergy <= 20) {
                printf("  INFO: pmtETotal > beam energy at eventIdx %lld: "
                       "beam = %.3f GeV, sumHitE = %.4f GeV, "
                       "pmtETotal = %.4f GeV, nHits = %d\n",
                       iEntry, maxSumEvP, sumHitE, out_pmtETotal,
                       out_nHitsDet9);
                // Print each hit with particle info
                for (Int_t ih = 0; ih < nDbgHits; ih++) {
                    const char *pname = "other";
                    switch (dbgHitPid[ih]) {
                        case  11:  pname = "e-";    break;
                        case -11:  pname = "e+";    break;
                        case  22:  pname = "gamma"; break;
                        case 211:  pname = "pi+";   break;
                        case -211: pname = "pi-";   break;
                        case 2212: pname = "p";     break;
                    }
                    printf("    hit[%d]: %s (pid=%d, trid=%d), "
                           "E=%.4f GeV, x=%.2f cm, y=%.2f cm\n",
                           ih, pname, dbgHitPid[ih], dbgHitTrid[ih],
                           dbgHitE[ih], dbgHitX[ih], dbgHitY[ih]);
                }
                if (nDbgHits < out_nHitsDet9)
                    printf("    ... (%d more hits not shown)\n",
                           out_nHitsDet9 - nDbgHits);
            }
            if (nOverBeamEnergy == 20)
                printf("  (suppressing further beam energy warnings)\n");
        }

        // --- Per-PMT energy smearing ---
        // Applied after conservation checks to model independent
        // sampling fluctuations and photostatistics in each PMT.
        // sigma/E = 12.9%/sqrt(E) + 1.2% (E in GeV)
        // Negative results clamped to zero.
        for (Int_t s = 0; s < DetGeom::NPMT; s++) {
            if (out_pmtE[s] > 0.0) {
                Double_t sigma = CaloResolution::Sigma(out_pmtE[s]);
                out_pmtE[s] = rng->Gaus(out_pmtE[s], sigma);
                if (out_pmtE[s] < 0.0) out_pmtE[s] = 0.0;
            }
        }

        // Recompute sums from smeared values
        out_pmtELeft  = 0.0;
        out_pmtERight = 0.0;
        for (Int_t s = 0; s < 4; s++) out_pmtELeft  += out_pmtE[s];
        for (Int_t s = 4; s < 8; s++) out_pmtERight += out_pmtE[s];
        out_pmtETotal = out_pmtELeft + out_pmtERight;

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
    if (tDiag)
        printf("Showering events saved:  %lld / %d max\n",
               nDiagEventsFilled, MAX_DIAG_EVENTS);

    // Write and close
    fOut->cd();
    tOut->Write();
    if (tDiag) tDiag->Write();
    fOut->Close();

    printf("Output written to: %s\n", outName.Data());

    delete T;
    delete rng;
}

#endif // MOLPOL_DETECTOR_RESPONSE_C