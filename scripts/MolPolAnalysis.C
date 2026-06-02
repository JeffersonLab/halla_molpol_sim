///////////////////////////////////////////////////////////////
// MolPolAnalysis.C -- Eric King 06/02/2026 (updated)
// ROOT macro for analyzing MolPol simulation data: coincidence,
// left singles, and right singles weighted asymmetry / analyzing power.
//
// Usage (from this directory):
//   root 'MolPolAnalysis.C("path/to/file.root")'
//   root 'MolPolAnalysis.C("filelist.txt", false)'
//   root 'MolPolAnalysis.C("file.root", false, 2.0)'
///////////////////////////////////////////////////////////////

#include "MolPolAnalysis.h"
#include "MollerCrossSection.C"

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TPaveStats.h"

#include <iostream>
#include <iomanip>
#include <cmath>

// Channel enumeration for the asymmetry analysis, maps to integer values for the accumulator array.
enum MolPolAsymChannel {
    kMolPolAsymCoinc = 0,
    kMolPolAsymLeft,
    kMolPolAsymRight,
    kMolPolAsymN
};

// Function to return common names for the asymmetry analysis channels
const char *MolPolAsymChannelName(Int_t ch) {
    switch (ch) {
        case kMolPolAsymCoinc: return "Coincidence";
        case kMolPolAsymLeft:  return "Left singles";
        case kMolPolAsymRight: return "Right singles";
        default:               return "?";
    }
}

// Structure to accumulate the asymmetry analysis results
struct MolPolAsymAccum {
    Int_t    nPass;           // unweighted count of accepted MC events
    Double_t weightAccepted;  // sum of evUnpolWght for accepted events
    Double_t sumMollerEvXs;   // sum of evXs for accepted events
    Double_t summedAsym;
    Double_t summedAsym2;
    Double_t sumPolWtPos;
    Double_t sumPolWtNeg;
};

// Function to accumulate the asymmetry analysis results for a given channel
void MolPolAccumulateChannel(MolPolAsymAccum &a, Double_t eventAzz) {
    a.nPass++;
    a.weightAccepted += evUnpolWght;
    a.sumMollerEvXs += evXs;
    a.summedAsym  += eventAzz * (evPolPlusWghtZ + evPolMinusWghtZ);
    a.summedAsym2 += eventAzz * eventAzz * (evPolPlusWghtZ + evPolMinusWghtZ);
    a.sumPolWtPos += evPolPlusWghtZ;
    a.sumPolWtNeg += evPolMinusWghtZ;
}

// Place the ROOT stats box in NDC (lower-left to upper-right corners).
void MolPolSetStatsBoxNDC(TH1 *h, Double_t x1, Double_t y1, Double_t x2, Double_t y2) {
    if (!h || !gPad) return;
    gPad->Update();
    TPaveStats *stats = (TPaveStats *)h->GetListOfFunctions()->FindObject("stats");
    if (stats) {
        stats->SetX1NDC(x1);
        stats->SetY1NDC(y1);
        stats->SetX2NDC(x2);
        stats->SetY2NDC(y2);
        gPad->Modified();
        gPad->Update();
    }
}

// Main macro function to analyze the asymmetry data
// useTrackGates: true (default) = hitFluxTrack* ; false = hitFluxP* with momentumCut
void MolPolAnalysis(const char *fileList, 
                    Bool_t useTrackGates = kTRUE, 
                    Double_t momentumCut = 0.0) {

    TChain *T = SetupMolpolChain(fileList);
    if (!T) return;

    SetupMolpolBranches(T);

    gStyle->SetPalette(55);

    // PMT row configuration (all active by default), deactivate PMT row with function call below.
     SetPmtRowInactive(4);SetPmtRowInactive(3);

    const Double_t targetPolarization = 0.08014;
    const MolPolGateMode gateMode =
        useTrackGates ? kMolPolGateTrack : kMolPolGateMomentum;

    printf("=== MolPolAnalysis ===\n");
    printf("File: %s\n", fileList);
    printf("Gate mode: %s\n", MolPolGateModeName(gateMode));
    if (gateMode == kMolPolGateMomentum)
        printf("Momentum cut: %g GeV\n", momentumCut);
    PrintActiveRows();

    // Declare the accumulator array for the asymmetry analysis
    MolPolAsymAccum accum[kMolPolAsymN];
    // Initialize the accumulator array values to zero
    for (Int_t ch = 0; ch < kMolPolAsymN; ch++) {
        accum[ch].nPass = 0;
        accum[ch].weightAccepted = 0.0;
        accum[ch].sumMollerEvXs = 0.0;
        accum[ch].summedAsym = 0.0;
        accum[ch].summedAsym2 = 0.0;
        accum[ch].sumPolWtPos = 0.0;
        accum[ch].sumPolWtNeg = 0.0;
    }

    // Total weights over all MC entries in the chain
    Double_t weightTotalThrown = 0.0;
    Double_t totalMollerEvXs = 0.0;

    // Declare the variables to track the min and max values of evThcom and evPhcom to determine the
    Double_t maxThcom = 0.0;
    Double_t minThcom = 180.0;
    Double_t maxPhcom = 0.0;
    Double_t minPhcom = 360.0;
    Double_t maxSumEvP= 0.0;

    // Histograms for the asymmetry analysis
    TH1F *hAzzCoinc = new TH1F("hAzzCoinc", "Analyzing Power for Coincidence;Analyzing Power;Weighted Events", 360, 0.0, 0.777777778);
    TH1F *hAzzLeft = new TH1F("hAzzLeft", "Analyzing Power for Left;Analyzing Power;Weighted Events", 360, 0.0, 0.777777778);
    TH1F *hAzzRight = new TH1F("hAzzRight", "Analyzing Power for Right;Analyzing Power;Weighted Events", 360, 0.0, 0.777777778);

    TH2F *hThetaPhiCoinc = new TH2F("hThetaPhiCoinc", "MC #theta_{COM} vs. #phi for Coincidence;#theta_{COM};#phi", 360, 0.0, 180.0, 100, -180., 180.);
    TH2F *hThetaPhiLeft = new TH2F("hThetaPhiLeft", "MC #theta_{COM} vs. #phi for Left;#theta_{COM};#phi", 360, 0.0, 180.0, 100, -180., 180.);
    TH2F *hThetaPhiRight = new TH2F("hThetaPhiRight", "MC #theta_{COM} vs. #phi for Right;#theta_{COM};#phi", 360, 0.0, 180.0, 100, -180., 180.);

    Double_t hitLxMin = -9.0;
    Double_t hitLxMax = 9.0;
    Double_t hitLyMin = -15.0;
    Double_t hitLyMax = 15.0;
    Int_t    hitLxBin = 90;
    Int_t    hitLyBin = 150;
    // Hit detector 9 Lx vs Ly maps for singles and coincidence
    TH2F *hHitLxLyETotalSinglesMap = new TH2F("hHitLxLyETotalSinglesMap", "All Events (#SigmaE Fraction [%]);Lx [cm];Ly [cm]", hitLxBin, hitLxMin, hitLxMax, hitLyBin, hitLyMin, hitLyMax);
    TH2F *hHitLxLyETotalCoincMap = new TH2F("hHitLxLyETotalCoincMap", "Coincidence Events (#SigmaE Fraction [%]);Lx [cm];Ly [cm]", hitLxBin, hitLxMin, hitLxMax, hitLyBin, hitLyMin, hitLyMax);
    TH2F *hHitLxLyEBarSinglesMap = new TH2F("hHitLxLyEBarSinglesMap", "All Events (#bar{E} Per Bin [GeV]);Lx [cm];Ly [cm]", hitLxBin, hitLxMin, hitLxMax, hitLyBin, hitLyMin, hitLyMax);
    TH2F *hHitLxLyEBarCoincMap = new TH2F("hHitLxLyEBarCoincMap", "Coincidence Events (#bar{E} Per Bin [GeV]);Lx [cm];Ly [cm]", hitLxBin, hitLxMin, hitLxMax, hitLyBin, hitLyMin, hitLyMax);

    // Detector 9 momentum flux (all entries; independent of analysis gate mode)
    Float_t fluxPMin = 0.0;
    Float_t fluxPMax = 11.0;
    Int_t   fluxPBin = 220;
    TH1F *hFluxPLeft = new TH1F("hFluxPLeft",
        "Momentum flux on detector 9 (left);#Sigma |p| (GeV);Weighted events", fluxPBin, fluxPMin, fluxPMax);
    TH1F *hFluxPRight = new TH1F("hFluxPRight",
        "Momentum flux on detector 9 (right);#Sigma |p| (GeV);Weighted events", fluxPBin, fluxPMin, fluxPMax);
    TH1F *hFluxPSum = new TH1F("hFluxPSum",
        "Momentum flux on detector 9 (left + right);#Sigma |p| (GeV);Weighted events", fluxPBin, fluxPMin, fluxPMax);
    TH2F *hFluxPLR = new TH2F("hFluxPLR",
        "Momentum flux on detector 9 (left vs right);#Sigma |p| left (GeV);#Sigma |p| right (GeV)",
        fluxPBin, 0.0, fluxPMax, fluxPBin, 0.0, fluxPMax);

    // Get the number of entries in the tree
    Long64_t nEntries = T->GetEntries();
    printf("Entries: %lld\n", nEntries);

    // Standard event loop over number of tree/chain entries
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        T->GetEntry(iEntry);

        // Collect the min and max values of evThcom and evPhcom to determine the
        // monte carlo'ed theta and phi range in order to calculate simulated cross section.
        if(evThcom > maxThcom) maxThcom = evThcom;
        if(evThcom < minThcom) minThcom = evThcom;
        if(evPhcom > maxPhcom) maxPhcom = evPhcom;
        if(evPhcom < minPhcom) minPhcom = evPhcom;
        // Here we need to calculated the sum of evP for events which result in a track-qualified coincidence.
        Double_t sumEvP = 0.0;
        if( hitFluxTrackCoinc() ){
            sumEvP += evP[0] + evP[1];
            if(sumEvP > maxSumEvP) maxSumEvP = sumEvP;
        }

        // Momentum flux on det 9 (non-neutron); filled for every entry, not gated
        //{
            Double_t fluxLeft  = hitFluxPLeftSum();
            Double_t fluxRight = hitFluxPRightSum();
            //fluxLeft = log10(fluxLeft);
            //fluxRight = log10(fluxRight);
            //if(fluxLeft + fluxRight < 0.1)
            //std::cout << "iEntry: " << iEntry << " fluxLeft: " << fluxLeft << " fluxRight: " << fluxRight << " sum: " << fluxLeft + fluxRight << std::endl;
            Double_t fluxCut = momentumCut;
            if(fluxLeft > fluxCut) hFluxPLeft->Fill(fluxLeft, 1);
            if(fluxRight > fluxCut) hFluxPRight->Fill(fluxRight, 1);
            if(fluxLeft + fluxRight > fluxCut) hFluxPSum->Fill(fluxLeft + fluxRight, 1);
            if(fluxLeft > fluxCut || fluxRight > fluxCut) hFluxPLR->Fill(fluxLeft, fluxRight, 1);
        //}

        // Total weights for all thrown MC events
        weightTotalThrown += evUnpolWght;
        totalMollerEvXs += evXs;
        // Calculate the analyzing power for the current simulated event
        Double_t eventAzz = MolPolEventAsymmetry();
        // Accumulate the analyzing power for the current simulated event for individual arms and coincidence
        if (MolPolPassCoinc(gateMode, momentumCut)) MolPolAccumulateChannel(accum[kMolPolAsymCoinc], eventAzz);
        if (MolPolPassLeft(gateMode, momentumCut))  MolPolAccumulateChannel(accum[kMolPolAsymLeft], eventAzz);
        if (MolPolPassRight(gateMode, momentumCut)) MolPolAccumulateChannel(accum[kMolPolAsymRight], eventAzz);

        // Fill the histograms with the analyzing power for the current simulated event
        if (MolPolPassCoinc(gateMode, momentumCut) && eventAzz > 0.0) hAzzCoinc->Fill(eventAzz, (evPolPlusWghtZ + evPolMinusWghtZ));
        if (MolPolPassLeft(gateMode, momentumCut) && eventAzz > 0.0) hAzzLeft->Fill(eventAzz, (evPolPlusWghtZ + evPolMinusWghtZ));
        if (MolPolPassRight(gateMode, momentumCut) && eventAzz > 0.0) hAzzRight->Fill(eventAzz, (evPolPlusWghtZ + evPolMinusWghtZ));

        // Fill the energy deposit heat maps for coincidence and singles events
        if (MolPolPassCoinc(gateMode, momentumCut)) {
            for (Int_t i = 0; i < hitN; i++) {
                if(IsActivePmtRow(i)) hHitLxLyETotalCoincMap->Fill(hitLx[i]*100.0, hitLy[i]*100.0, hitE[i]);
                if(IsActivePmtRow(i)) hHitLxLyEBarCoincMap->Fill(hitLx[i]*100.0, hitLy[i]*100.0);
            }
        }
        if (MolPolPassLeft(gateMode, momentumCut) || MolPolPassRight(gateMode, momentumCut)) {
            for (Int_t i = 0; i < hitN; i++) {
                if(IsActivePmtRow(i)) hHitLxLyETotalSinglesMap->Fill(hitLx[i]*100.0, hitLy[i]*100.0, hitE[i]);
                if(IsActivePmtRow(i)) hHitLxLyEBarSinglesMap->Fill(hitLx[i]*100.0, hitLy[i]*100.0);
            }
        }

        // Fill the histograms with the MC #theta_{COM} and #phi for the current simulated event
        if (MolPolPassCoinc(gateMode, momentumCut)) hThetaPhiCoinc->Fill(evThcom, evPhcom);
        if (MolPolPassLeft(gateMode, momentumCut))  hThetaPhiLeft->Fill(evThcom, evPhcom);
        if (MolPolPassRight(gateMode, momentumCut)) hThetaPhiRight->Fill(evThcom, evPhcom);

    }
    
    // Calculate the #bar{E} Per Bin for the hitLxLySinglesCountMap and hitLxLyCoincCountMap
    hHitLxLyEBarSinglesMap->Divide(hHitLxLyETotalSinglesMap,hHitLxLyEBarSinglesMap,1.0,1.0);
    hHitLxLyEBarCoincMap->Divide(hHitLxLyETotalCoincMap,hHitLxLyEBarCoincMap,1.0,1.0);

    // Scale the histograms by the sum of the polarized weights for the current simulated event
    hAzzCoinc->Scale(1.0/(evPolPlusWghtZ + evPolMinusWghtZ));
    hAzzLeft->Scale(1.0/(evPolPlusWghtZ + evPolMinusWghtZ));
    hAzzRight->Scale(1.0/(evPolPlusWghtZ + evPolMinusWghtZ));
    
    // Histogram axis ranges (1D: x only; 2D: x and y)
    MolPolSetXAxisFromContent1D(hAzzCoinc);
    MolPolSetXAxisFromContent1D(hAzzLeft);
    MolPolSetXAxisFromContent1D(hAzzRight);
    MolPolSetAxesFromContent2D(hThetaPhiCoinc);
    MolPolSetAxesFromContent2D(hThetaPhiLeft);
    MolPolSetAxesFromContent2D(hThetaPhiRight);
    MolPolSetXAxisFromContent1D(hFluxPLeft);
    MolPolSetXAxisFromContent1D(hFluxPRight);
    MolPolSetXAxisFromContent1D(hFluxPSum);
    MolPolSetAxesFromContent2D(hFluxPLR);

    TCanvas *C1  = new TCanvas("C1", "C1", 1400, 400);
    C1->Divide(3, 1);
    const Double_t statX1 = 0.4, statY1 = 0.75, statX2 = 0.6, statY2 = 0.9;
    C1->cd(1);
    hAzzCoinc->SetMinimum(0.0);
    hAzzCoinc->Draw("HISTE");
    MolPolSetStatsBoxNDC(hAzzCoinc, statX1, statY1, statX2, statY2);
    C1->cd(2);
    hAzzLeft->SetMinimum(0.0);
    hAzzLeft->Draw("HISTE");
    MolPolSetStatsBoxNDC(hAzzLeft, statX1, statY1, statX2, statY2);
    C1->cd(3);
    hAzzRight->SetMinimum(0.0);
    hAzzRight->Draw("HISTE");
    MolPolSetStatsBoxNDC(hAzzRight, statX1, statY1, statX2, statY2);

    TCanvas *C2  = new TCanvas("C2", "C2", 1400, 400);
    C2->Divide(3, 1);
    C2->cd(1);
    hThetaPhiCoinc->SetMinimum(0.0);
    hThetaPhiCoinc->Draw("COLZ");
    C2->cd(2);
    hThetaPhiLeft->SetMinimum(0.0);
    hThetaPhiLeft->Draw("COLZ");
    C2->cd(3);
    hThetaPhiRight->SetMinimum(0.0);
    hThetaPhiRight->Draw("COLZ");

    TCanvas *C3 = new TCanvas("C3", "C3", 1400, 400);
    C3->Divide(3, 1);
    C3->cd(1);
    hFluxPLeft->SetMinimum(0.0);
    hFluxPLeft->Draw("HISTE");
    MolPolSetStatsBoxNDC(hFluxPLeft, statX1, statY1, statX2, statY2);
    C3->cd(2);
    hFluxPRight->SetMinimum(0.0);
    hFluxPRight->Draw("HISTE");
    MolPolSetStatsBoxNDC(hFluxPRight, statX1, statY1, statX2, statY2);
    C3->cd(3);
    hFluxPSum->SetMinimum(0.0);
    hFluxPSum->Draw("HISTE");
    MolPolSetStatsBoxNDC(hFluxPSum, statX1, statY1, statX2, statY2);

    TCanvas *C4 = new TCanvas("C4", "C4", 600, 550);
    hFluxPLR->SetMinimum(0.0);
    hFluxPLR->Draw("COLZ");

    TH2F *hHitLxLyETotalSinglesPct = (TH2F *)hHitLxLyETotalSinglesMap->Clone("hHitLxLyETotalSinglesPct");
    TH2F *hHitLxLyETotalCoincPct = (TH2F *)hHitLxLyETotalCoincMap->Clone("hHitLxLyETotalCoincPct");
    const Double_t intETotalSingles = hHitLxLyETotalSinglesPct->Integral();
    if (intETotalSingles > 0.0) hHitLxLyETotalSinglesPct->Scale(100.0 / intETotalSingles);
    const Double_t intETotalCoinc = hHitLxLyETotalCoincPct->Integral();
    if (intETotalCoinc > 0.0) hHitLxLyETotalCoincPct->Scale(100.0 / intETotalCoinc);
    hHitLxLyETotalCoincPct->SetMaximum(hHitLxLyETotalSinglesPct->GetMaximum());

    TCanvas *C5 = new TCanvas("C5", "C5", 1600, 800);
    C5->Divide(4, 1);
    C5->cd(1);
    gPad->SetRightMargin(0.15);
    hHitLxLyETotalSinglesPct->SetStats(kFALSE);
    hHitLxLyETotalSinglesPct->SetMinimum(0.0);
    hHitLxLyETotalSinglesPct->Draw("COLZ");
    gPad->SetGridx(kTRUE);
    gPad->SetGridy(kTRUE);
    C5->cd(2);
    gPad->SetRightMargin(0.15);
    hHitLxLyETotalCoincPct->SetStats(kFALSE);
    hHitLxLyETotalCoincPct->SetMinimum(0.0);
    hHitLxLyETotalCoincPct->Draw("COLZ");
    gPad->SetGridx(kTRUE);
    gPad->SetGridy(kTRUE);
    C5->cd(3);
    gPad->SetRightMargin(0.15);
    hHitLxLyEBarSinglesMap->SetStats(kFALSE);
    hHitLxLyEBarSinglesMap->SetMinimum(0.0);
    hHitLxLyEBarSinglesMap->Draw("COLZ");
    gPad->SetGridx(kTRUE);
    gPad->SetGridy(kTRUE);
    C5->cd(4);
    gPad->SetRightMargin(0.15);
    hHitLxLyEBarCoincMap->SetStats(kFALSE);
    hHitLxLyEBarCoincMap->SetMinimum(0.0);
    hHitLxLyEBarCoincMap->SetMaximum( hHitLxLyEBarSinglesMap->GetMaximum() );
    hHitLxLyEBarCoincMap->Draw("COLZ");
    gPad->SetGridx(kTRUE);
    gPad->SetGridy(kTRUE);

    // mollerCrossSection.C: Fe rate using MC theta range (track-qualified sumEvP uses hitFluxTrackCoinc elsewhere)
    Double_t molrate = RunMollerCrossSection(kMollerRateFe, minThcom, maxThcom, fabs(maxPhcom - minPhcom), maxSumEvP, 1.0, kTRUE);

    printf("\n========= Kinematics =========\n");
    printf("Beam energy (maxSumEvP): %.3f GeV\n", maxSumEvP);
    printf("Theta_CM range: [%.3f, %.3f] deg\n", minThcom, maxThcom);
    printf("Phi_CM range:   [%.3f, %.3f] deg\n", minPhcom, maxPhcom);
    printf("Moller rate:    %.3f events/s/muA/um\n", molrate);

    // Print the summary of asymmetry analysis results for each channel
    printf("\n========= Summary =========\n");
    printf("Entries: %lld\n", nEntries);

    for (Int_t ch = 0; ch < kMolPolAsymN; ch++) {
        printf("--- %s ---\n", MolPolAsymChannelName(ch));

        Double_t acceptedEventFrac = 0.0;
        if (nEntries > 0) acceptedEventFrac = (Double_t)accum[ch].nPass / (Double_t)nEntries;

        Double_t weightedEventFrac = 0.0;
        if (weightTotalThrown > 0.0) weightedEventFrac = accum[ch].weightAccepted / weightTotalThrown;

        Double_t mollerWtdEventFrac = 0.0;
        if (totalMollerEvXs > 0.0) mollerWtdEventFrac = accum[ch].sumMollerEvXs / totalMollerEvXs;

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
            polarizationFacErr = sqrt(fabs(varAzz) / (Double_t)accum[ch].nPass) / fabs(meanAzz) * fabs(polarizationFactor);
            meanAnalyzingPower = 1.0 / polarizationFactor / targetPolarization;
            meanAnalyzingPowerErr = polarizationFacErr / polarizationFactor / polarizationFactor / targetPolarization;
        }

        printf("    Polarization factor: %.5f +/- %.5f\n",
               polarizationFactor, polarizationFacErr);
        printf("   Mean analyzing power: %.7f +/- %.7f\n",
               meanAnalyzingPower, meanAnalyzingPowerErr);
    }

    delete T;
}