///////////////////////////////////////////////////////////////
// MolPolDetectorResponseAnalysis.h
// Eric King 06/02/2026
//
// Helper types, functions, and globals for the detector
// response analysis script.
///////////////////////////////////////////////////////////////

#ifndef MOLPOL_DETECTOR_RESPONSE_ANALYSIS_H
#define MOLPOL_DETECTOR_RESPONSE_ANALYSIS_H

#include "MolPolDetectorResponse.h"

#include "TString.h"
#include "TMath.h"

#include <cmath>

///////////////////////////////////////////////////////////////
// Active PMT row parsing
// Input string contains digits 1-4 for active rows.
// E.g., "1234" = all rows, "23" = inner rows, "14" = outer.
// Returns a boolean array activeRow[4] (0-indexed: row 1 = [0]).
///////////////////////////////////////////////////////////////

inline void ParseActiveRows(const char *rowStr, Bool_t activeRow[4]) {
    for (Int_t r = 0; r < 4; r++) activeRow[r] = false;
    TString s(rowStr);
    for (Int_t i = 0; i < s.Length(); i++) {
        Int_t row = s[i] - '1';
        if (row >= 0 && row < 4) activeRow[row] = true;
    }
}

///////////////////////////////////////////////////////////////
// Asymmetry accumulator
///////////////////////////////////////////////////////////////

struct DetRespAsymAccum {
    Int_t    nPass;
    Double_t weightAccepted;
    Double_t sumMollerEvXs;
    Double_t summedAsym;
    Double_t summedAsym2;
    Double_t sumPolWtPos;
    Double_t sumPolWtNeg;
};

inline void DetRespAccumInit(DetRespAsymAccum &a) {
    a.nPass = 0;
    a.weightAccepted = 0.0;
    a.sumMollerEvXs = 0.0;
    a.summedAsym = 0.0;
    a.summedAsym2 = 0.0;
    a.sumPolWtPos = 0.0;
    a.sumPolWtNeg = 0.0;
}

inline void DetRespAccumulate(DetRespAsymAccum &a, Double_t eventAzz) {
    a.nPass++;
    a.weightAccepted += dr_evUnpolWght;
    a.sumMollerEvXs  += dr_evXs;
    a.summedAsym  += eventAzz * (dr_evPolPlusWghtZ + dr_evPolMinusWghtZ);
    a.summedAsym2 += eventAzz * eventAzz * (dr_evPolPlusWghtZ + dr_evPolMinusWghtZ);
    a.sumPolWtPos += dr_evPolPlusWghtZ;
    a.sumPolWtNeg += dr_evPolMinusWghtZ;
}

///////////////////////////////////////////////////////////////
// Channel enumeration
///////////////////////////////////////////////////////////////

enum DetRespAsymChannel {
    kDetRespAsymCoinc = 0,
    kDetRespAsymLeft,
    kDetRespAsymRight,
    kDetRespAsymN
};

inline const char *DetRespAsymChannelName(Int_t ch) {
    switch (ch) {
        case kDetRespAsymCoinc: return "Coincidence";
        case kDetRespAsymLeft:  return "Left singles";
        case kDetRespAsymRight: return "Right singles";
        default:                return "?";
    }
}

///////////////////////////////////////////////////////////////
// Solve90PctRadius()
// Bisection solver for the radius containing 90% of the
// two-component radial energy distribution at a given depth.
//
// CDF: F(r) = p * r^2/(r^2+Rc^2) + (1-p) * r^2/(r^2+Rt^2)
// Solves F(r) = 0.9 for r in units of R_M.
///////////////////////////////////////////////////////////////

inline Double_t Solve90PctRadius(Double_t p, Double_t Rc, Double_t Rt) {
    const Double_t target = 0.90;
    Double_t rLo = 0.0;
    Double_t rHi = 20.0;
    for (Int_t iter = 0; iter < 100; iter++) {
        Double_t rMid = 0.5 * (rLo + rHi);
        Double_t r2 = rMid * rMid;
        Double_t F = p * r2 / (r2 + Rc * Rc)
                   + (1.0 - p) * r2 / (r2 + Rt * Rt);
        if (F < target) rLo = rMid;
        else            rHi = rMid;
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

#endif // MOLPOL_DETECTOR_RESPONSE_ANALYSIS_H