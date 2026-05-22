#ifndef MOLPOL_ANALYSIS_H
#define MOLPOL_ANALYSIS_H

#include "MolPolData.h"
#include "TH1.h"
#include "TH2.h"

///////////////////////////////////////////////////////////////
// MolPolAnalysis.h
// Analysis helpers for MolPol ROOT trees (flux gates, PMT rows,
// event-level asymmetry, histogram axis helpers). Requires MolPolData.h (branch buffers
// filled after SetupMolpolBranches + GetEntry).
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
// MolPolEventAsymmetry()
// Polarized-weight asymmetry for the current entry (after GetEntry):
//   (evPolPlusWghtZ - evPolMinusWghtZ) / (evPolPlusWghtZ + evPolMinusWghtZ).
// Returns 0 if the denominator is zero.
///////////////////////////////////////////////////////////////

inline Double_t MolPolEventAsymmetry() {
    Double_t numer = evPolPlusWghtZ - evPolMinusWghtZ;
    Double_t denom = evPolPlusWghtZ + evPolMinusWghtZ;
    if (denom == 0.0) return 0.0;
    return numer / denom;
}

///////////////////////////////////////////////////////////////
// PMT row active state
// Detector 9 is segmented into 4 rows by hitLy.
// Rows are always active in L/R pairs (R/L looking upstream):
//   Row 1 (PMT 1/5):  hitLy >= 0.075
//   Row 2 (PMT 2/6):  0.0 <= hitLy < 0.075
//   Row 3 (PMT 3/7): -0.075 <= hitLy < 0.0
//   Row 4 (PMT 4/8):  hitLy < -0.075
///////////////////////////////////////////////////////////////

Bool_t pmtRow1Active = true;
Bool_t pmtRow2Active = true;
Bool_t pmtRow3Active = true;
Bool_t pmtRow4Active = true;

inline void SetPmtRowActive(Int_t row) {
    switch (row) {
        case 1: pmtRow1Active = true; break;
        case 2: pmtRow2Active = true; break;
        case 3: pmtRow3Active = true; break;
        case 4: pmtRow4Active = true; break;
    }
}

inline void SetPmtRowInactive(Int_t row) {
    switch (row) {
        case 1: pmtRow1Active = false; break;
        case 2: pmtRow2Active = false; break;
        case 3: pmtRow3Active = false; break;
        case 4: pmtRow4Active = false; break;
    }
}

///////////////////////////////////////////////////////////////
// IsActivePmtRow()
// Given a hit array index, determines which PMT row the hit
// falls in based on hitLy and returns whether that row is active.
///////////////////////////////////////////////////////////////

inline Bool_t IsActivePmtRow(Int_t i) {
    Double_t ly = hitLy[i];
    if (ly >= 0.075)             return pmtRow1Active;
    if (ly >= 0.0)               return pmtRow2Active;
    if (ly >= -0.075)            return pmtRow3Active;
                                 return pmtRow4Active;
}

///////////////////////////////////////////////////////////////
// hitFluxPLeftSum() / hitFluxPRightSum()
// Summed hit |p| on detector 9, excluding neutrons (pid != 2112),
// gated on active PMT rows.
// Left/Right defined looking upstream: left = hitLx >= 0,
// right = hitLx < 0.
///////////////////////////////////////////////////////////////

inline Double_t hitFluxPLeftSum() {
    Double_t sum = 0.0;
    for (Int_t i = 0; i < hitN; i++) {
        if (hitDet[i] == 9 && hitPid[i] < 2112 && hitLx[i] >= 0.0 && IsActivePmtRow(i))
            sum += hitP[i];
    }
    return sum;
}

inline Double_t hitFluxPRightSum() {
    Double_t sum = 0.0;
    for (Int_t i = 0; i < hitN; i++) {
        if (hitDet[i] == 9 && hitPid[i] < 2112 && hitLx[i] < 0.0 && IsActivePmtRow(i))
            sum += hitP[i];
    }
    return sum;
}

///////////////////////////////////////////////////////////////
// hitFluxPLeft(pMin) / hitFluxPRight(pMin)
// True if summed flux on that arm meets the threshold (uses
// hitFluxPLeftSum / hitFluxPRightSum). For pMin <= 0, require sum > 0.
// For pMin > 0, require sum >= pMin.
///////////////////////////////////////////////////////////////

inline Bool_t hitFluxPLeft(Double_t pMin = 0.0) {
    Double_t pSum = hitFluxPLeftSum();
    if (pMin <= 0.0) return (pSum > 0.0);
    return (pSum >= pMin);
}

inline Bool_t hitFluxPRight(Double_t pMin = 0.0) {
    Double_t pSum = hitFluxPRightSum();
    if (pMin <= 0.0) return (pSum > 0.0);
    return (pSum >= pMin);
}

///////////////////////////////////////////////////////////////
// hitFluxPCoinc(cut)
// Momentum coincidence on both arms: delegates threshold test to
// hitFluxPLeft(cut) && hitFluxPRight(cut).
///////////////////////////////////////////////////////////////

inline Bool_t hitFluxPCoinc(Double_t pMin = 0.0) {
    return hitFluxPLeft(pMin) && hitFluxPRight(pMin);
}

///////////////////////////////////////////////////////////////
// hitFluxTrackLeft() / hitFluxTrackRight()
// True if detector 9 has at least one electron (pid == 11) on that
// arm (hitLx) with hitTrid == 1 or 2, on an active PMT row.
///////////////////////////////////////////////////////////////

inline Bool_t hitFluxTrackLeft() {
    for (Int_t i = 0; i < hitN; i++) {
        if (hitDet[i] == 9 && hitPid[i] == 11 && IsActivePmtRow(i) && hitLx[i] >= 0.0) {
            if (hitTrid[i] == 1 || hitTrid[i] == 2) return true;
        }
    }
    return false;
}

inline Bool_t hitFluxTrackRight() {
    for (Int_t i = 0; i < hitN; i++) {
        if (hitDet[i] == 9 && hitPid[i] == 11 && IsActivePmtRow(i) && hitLx[i] < 0.0) {
            if (hitTrid[i] == 1 || hitTrid[i] == 2) return true;
        }
    }
    return false;
}

///////////////////////////////////////////////////////////////
// hitFluxTrackCoinc()
// True if both spectrometer arms see track-qualified electron hits:
// hitFluxTrackLeft() && hitFluxTrackRight().
///////////////////////////////////////////////////////////////

inline Bool_t hitFluxTrackCoinc() {
    return hitFluxTrackLeft() && hitFluxTrackRight();
}

///////////////////////////////////////////////////////////////
// MolPolGateMode — select track-qualified vs momentum-cut gates
///////////////////////////////////////////////////////////////

enum MolPolGateMode {
    kMolPolGateTrack = 0,
    kMolPolGateMomentum
};

inline const char *MolPolGateModeName(MolPolGateMode mode) {
    if (mode == kMolPolGateTrack) return "track-qualified";
    return "momentum-cut";
}

inline Bool_t MolPolPassCoinc(MolPolGateMode mode, Double_t pMin = 0.0) {
    if (mode == kMolPolGateTrack) return hitFluxTrackCoinc();
    return hitFluxPCoinc(pMin);
}

inline Bool_t MolPolPassLeft(MolPolGateMode mode, Double_t pMin = 0.0) {
    if (mode == kMolPolGateTrack) return hitFluxTrackLeft();
    return hitFluxPLeft(pMin);
}

inline Bool_t MolPolPassRight(MolPolGateMode mode, Double_t pMin = 0.0) {
    if (mode == kMolPolGateTrack) return hitFluxTrackRight();
    return hitFluxPRight(pMin);
}

///////////////////////////////////////////////////////////////
// MolPolSetXAxisFromContent1D()
// Set x-axis range from five bins before the first non-zero bin to
// five bins after the last non-zero bin (clamped to histogram edges).
// No-op if every bin is empty or if the range is invalid.
///////////////////////////////////////////////////////////////

inline void MolPolSetXAxisFromContent1D(TH1 *h) {
    if (!h) return;
    Int_t nb = h->GetNbinsX();
    Int_t first = -1;
    Int_t last  = -1;
    for (Int_t b = 1; b <= nb; b++) {
        if (h->GetBinContent(b) > 0.0) {
            if (first < 0) first = b;
            last = b;
        }
    }
    if (first < 0 || last < 0) return;

    const Int_t margin = 5;
    Int_t lowBin  = first - margin;
    if (lowBin < 1) lowBin = 1;
    Int_t highBin = last + margin;
    if (highBin > nb) highBin = nb;

    Double_t xmin = h->GetXaxis()->GetBinLowEdge(lowBin);
    Double_t xmax = h->GetXaxis()->GetBinUpEdge(highBin);
    if (xmin >= xmax) return;
    h->GetXaxis()->SetRangeUser(xmin, xmax);
}

///////////////////////////////////////////////////////////////
// MolPolSetAxesFromContent2D()
// For each axis, set range from five bins before the first bin that
// has non-zero content (in any slice along the other axis) to five bins
// after the last such bin. Clamped to histogram edges. No-op if the
// histogram is empty or a range is invalid.
///////////////////////////////////////////////////////////////

inline void MolPolSetAxesFromContent2D(TH2 *h) {
    if (!h) return;
    const Int_t margin = 5;
    Int_t nbx = h->GetNbinsX();
    Int_t nby = h->GetNbinsY();

    Int_t firstX = -1, lastX = -1;
    for (Int_t bx = 1; bx <= nbx; bx++) {
        Bool_t colHas = false;
        for (Int_t by = 1; by <= nby; by++) {
            if (h->GetBinContent(bx, by) > 0.0) {
                colHas = true;
                break;
            }
        }
        if (colHas) {
            if (firstX < 0) firstX = bx;
            lastX = bx;
        }
    }
    if (firstX < 0 || lastX < 0) return;

    Int_t firstY = -1, lastY = -1;
    for (Int_t by = 1; by <= nby; by++) {
        Bool_t rowHas = false;
        for (Int_t bx = 1; bx <= nbx; bx++) {
            if (h->GetBinContent(bx, by) > 0.0) {
                rowHas = true;
                break;
            }
        }
        if (rowHas) {
            if (firstY < 0) firstY = by;
            lastY = by;
        }
    }
    if (firstY < 0 || lastY < 0) return;

    Int_t lowX  = firstX - margin;
    if (lowX < 1) lowX = 1;
    Int_t highX = lastX + margin;
    if (highX > nbx) highX = nbx;

    Int_t lowY  = firstY - margin;
    if (lowY < 1) lowY = 1;
    Int_t highY = lastY + margin;
    if (highY > nby) highY = nby;

    Double_t xmin = h->GetXaxis()->GetBinLowEdge(lowX);
    Double_t xmax = h->GetXaxis()->GetBinUpEdge(highX);
    Double_t ymin = h->GetYaxis()->GetBinLowEdge(lowY);
    Double_t ymax = h->GetYaxis()->GetBinUpEdge(highY);
    if (xmin >= xmax || ymin >= ymax) return;
    h->GetXaxis()->SetRangeUser(xmin, xmax);
    h->GetYaxis()->SetRangeUser(ymin, ymax);
}

#endif // MOLPOL_ANALYSIS_H
