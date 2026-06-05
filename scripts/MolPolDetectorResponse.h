#ifndef MOLPOL_DETECTOR_RESPONSE_H
#define MOLPOL_DETECTOR_RESPONSE_H

///////////////////////////////////////////////////////////////
// MolPolDetectorResponse.h
// Detector response model for the MolPol lead-fiber calorimeter.
// Implements electromagnetic shower development in Pb using the
// Grindhammer-Peters parameterization and distributes deposited
// energy across the 2x4 PMT segment grid.
//
// Contents:
//   PbConst     – Physical constants for pure lead
//   DetGeom     – Calorimeter segment geometry and helpers
//   GPShower    – Grindhammer-Peters longitudinal and radial
//                 shower parameterization (Appendix A.1)
//   SegIntegral – 2D Gauss-Legendre integration of the radial
//                 profile over rectangular PMT segments
//   ComputePmtFractions() – PMT segment energy fractions
//   ProcessHitShower()    – Full shower processing per hit
//   SetupDetRespChain()   – TChain loader for _detresponse.root
//   SetupDetRespBranches()– Branch address wiring for TDetResp
//
// Shower parameterization reference:
//   G. Grindhammer and S. Peters,
//   "The Parameterized Simulation of Electromagnetic Showers
//    in Homogeneous and Sampling Calorimeters",
//   arXiv:hep-ex/0001020 (2000).
//   https://arxiv.org/abs/hep-ex/0001020
//
// Requires: ROOT (TTree, TChain, TMath).
///////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"

#include <cmath>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <string>

///////////////////////////////////////////////////////////////
// Physical constants for pure lead (Pb, Z=82, A=207.2)
///////////////////////////////////////////////////////////////

namespace PbConst {
    const Double_t Z      = 82.0;
    const Double_t A      = 207.2;
    const Double_t rho    = 11.35;    // density [g/cm^3]
    const Double_t X0     = 0.5612;   // radiation length [cm] (NIST)
    const Double_t X0g    = 6.3696;   // radiation length [g/cm^2] = X0 * rho
    const Double_t Es     = 21.2;     // scale energy [MeV]
    // Critical energy [MeV] (Eq. 5): Ec = 2.66 * (X0[g/cm^2] * Z / A)^1.1
    // Note: Eq. 5 uses X0 in g/cm^2, not cm.
    // 2.66 * pow(6.3696 * 82.0 / 207.2, 1.1) = 2.66 * pow(2.5208, 1.1)
    //                                         = 2.66 * 2.7651
    //                                         = 7.3548 MeV
    const Double_t Ec     = 7.3548;   // MeV (precomputed to avoid Cling static init issues)
    // Moliere radius: RM = Es * X0[g/cm^2] / Ec  gives RM in g/cm^2,
    // then divide by rho to get cm.
    // RM = 21.2 * 6.3696 / 7.3548 / 11.35 = 1.6176 cm
    const Double_t RM     = 1.6176;   // Moliere radius [cm]
}

///////////////////////////////////////////////////////////////
// Detector geometry constants
// Calorimeter face: 30 cm tall, segmented into 2 (L/R) x 4 rows.
// Hit coordinates from GEANT4 are in meters; converted to cm
// internally by the calling code.
//
// PMT segment indexing:
//   L1=0 (top-left,    hitLy >= +7.5 cm,  hitLx >= 0)
//   L2=1 (mid-up-left, 0 <= hitLy < +7.5, hitLx >= 0)
//   L3=2 (mid-dn-left, -7.5 <= hitLy < 0, hitLx >= 0)
//   L4=3 (bot-left,    hitLy < -7.5,       hitLx >= 0)
//   R1=4, R2=5, R3=6, R4=7 (same rows, hitLx < 0)
///////////////////////////////////////////////////////////////

namespace DetGeom {
    const Int_t    NPMT     = 8;
    const Double_t blockHalfY = 15.0;   // [cm] half-height of block
    const Double_t blockHalfX = 9.0;    // [cm] half-width of block
    const Double_t detAngleDeg = 7.5;          // [deg] detector face angle from beam axis
    const Double_t detAngleRad = detAngleDeg * 3.14159265358979 / 180.0;  // [rad]
    const Double_t fiberDiam   = 0.1;          // [cm] scintillating fiber diameter (1 mm)
    const Double_t fiberX0     = 45.0;         // [cm] radiation length of fiber (polystyrene)
    const Double_t fiberVolFrac = 0.20;        // fiber-to-total volume fraction (1/5 fiber, 4/5 lead)
    // Y boundaries of the 4 rows (top to bottom), in cm
    const Double_t rowYBound[5] = { 15.0, 7.5, 0.0, -7.5, -15.0 };

    // Calorimeter longitudinal depth
    const Double_t caloDepth   = 30.0;  // [cm] physical depth of Pb block
    // Depth in radiation lengths, computed from the material X0
    const Double_t caloDepthX0 = caloDepth / PbConst::X0;  // = 53.46

    // Segment name strings for labels and histograms
    const char *PmtSegName[8] = {
        "L1", "L2", "L3", "L4",
        "R1", "R2", "R3", "R4"
    };
}

///////////////////////////////////////////////////////////////
// Calorimeter energy resolution
// σ/E = 12.9%/√E + 1.2%  (E in GeV)
//
// From: R. Livan, M. Vercesi, R. Wigmans,
//   "Scintillating-Fibre Calorimetry" (1995).
//
// This parameterization models the combined stochastic and
// constant resolution terms for a lead-fiber sampling
// calorimeter. Applied per-PMT segment to model independent
// sampling fluctuations and photostatistics in each segment.
///////////////////////////////////////////////////////////////

namespace CaloResolution {
    // Stochastic term: 12.9%
    const Double_t a = 0.129;
    // Constant term: 1.2%
    const Double_t b = 0.012;

    // Returns sigma in GeV for a given energy in GeV
    inline Double_t Sigma(Double_t E_GeV) {
        if (E_GeV <= 0.0) return 0.0;
        return E_GeV * (a / sqrt(E_GeV) + b);
    }
}

///////////////////////////////////////////////////////////////
// Grindhammer-Peters shower parameterization
// Appendix A.1 (homogeneous media) of arXiv:hep-ex/0001020
///////////////////////////////////////////////////////////////

namespace GPShower {

    //--- A.1.1: Average longitudinal profile parameters ---
    // T_hom = ln(y) - 0.858                       (Eq. 7, App A.1.1)
    // alpha_hom = 0.21 + (0.492 + 2.38/Z) * ln(y) (Eq. 8, App A.1.1)
    // where y = E / Ec  (E in MeV, Ec in MeV)

    inline Double_t Thom(Double_t lny) {
        return lny - 0.858;
    }

    inline Double_t AlphaHom(Double_t lny, Double_t Z) {
        return 0.21 + (0.492 + 2.38 / Z) * lny;
    }

    //--- Gamma distribution PDF for longitudinal profile (Eq. 2) ---
    // f(t) = (beta*t)^(alpha-1) * beta * exp(-beta*t) / Gamma(alpha)
    // beta = (alpha - 1) / T

    inline Double_t GammaPDF(Double_t t, Double_t alpha, Double_t beta) {
        if (t <= 0.0 || alpha <= 1.0) return 0.0;
        Double_t bt = beta * t;
        return TMath::Power(bt, alpha - 1.0) * beta * TMath::Exp(-bt) / TMath::Gamma(alpha);
    }

    //--- Incomplete gamma integral over [t1, t2] via numerical quadrature ---
    // Returns integral of f(t) dt over the interval, using Simpson's rule.
    inline Double_t GammaSliceIntegral(Double_t t1, Double_t t2,
                                       Double_t alpha, Double_t beta,
                                       Int_t nSteps = 20) {
        Double_t h = (t2 - t1) / nSteps;
        Double_t sum = GammaPDF(t1, alpha, beta) + GammaPDF(t2, alpha, beta);
        for (Int_t i = 1; i < nSteps; i++) {
            Double_t t = t1 + i * h;
            Double_t w = (i % 2 == 0) ? 2.0 : 4.0;
            sum += w * GammaPDF(t, alpha, beta);
        }
        return sum * h / 3.0;
    }

    //--- A.1.3: Average radial profile parameters (Eqs. 24-26) ---
    // All in Moliere radius units.
    // lnE = ln(E) where E is in MeV.

    // RC(tau) = z1 + z2 * tau                          (Eq. 24)
    inline Double_t RC(Double_t tau, Double_t lnE, Double_t Z) {
        Double_t z1 = 0.0251 + 0.00319 * lnE;
        Double_t z2 = 0.1162 - 0.000381 * Z;
        return z1 + z2 * tau;
    }

    // RT(tau) = k1 * { exp(k3*(tau-k2)) + exp(k4*(tau-k2)) }  (Eq. 25)
    inline Double_t RT(Double_t tau, Double_t lnE, Double_t Z) {
        Double_t k1 = 0.659 - 0.00309 * Z;
        Double_t k2 = 0.645;
        Double_t k3 = -2.59;
        Double_t k4 = 0.3585 + 0.0421 * lnE;
        Double_t dk = tau - k2;
        return k1 * (TMath::Exp(k3 * dk) + TMath::Exp(k4 * dk));
    }

    // p(tau) = p1 * exp( (p2-tau)/p3 - exp((p2-tau)/p3) )     (Eq. 26)
    inline Double_t CoreWeight(Double_t tau, Double_t lnE, Double_t Z) {
        Double_t p1 = 2.632 - 0.00094 * Z;
        Double_t p2 = 0.401 + 0.00187 * Z;
        Double_t p3 = 1.313 - 0.0686 * lnE;
        Double_t arg = (p2 - tau) / p3;
        Double_t val = p1 * TMath::Exp(arg - TMath::Exp(arg));
        // Clamp to [0, 1]
        if (val < 0.0) return 0.0;
        if (val > 1.0) return 1.0;
        return val;
    }

    //--- Radial CDF for a single component (Eq. 27) ---
    // F(r) = r^2 / (r^2 + R^2)
    // where r and R are both in Moliere units.
    // R stands for RC or RT as appropriate.
    inline Double_t RadialCDF(Double_t r, Double_t R) {
        if (r <= 0.0) return 0.0;
        return (r * r) / (r * r + R * R);
    }

    //--- Two-component radial CDF (Eq. 23 integrated) ---
    // F_total(r) = p * r^2/(r^2+RC^2) + (1-p) * r^2/(r^2+RT^2)
    inline Double_t RadialCDFTotal(Double_t r, Double_t p,
                                   Double_t Rc, Double_t Rt) {
        return p * RadialCDF(r, Rc) + (1.0 - p) * RadialCDF(r, Rt);
    }

} // namespace GPShower

///////////////////////////////////////////////////////////////
// Segment integration over finite block dimensions
//
// Computes the fraction of shower energy deposited in each
// PMT segment using a semi-analytic method over the physical
// block boundaries (±9 cm in x, row boundaries in y).
//
// For each segment rectangle [x1,x2] × [y1,y2]:
//   - Inner integral (y): analytic antiderivative
//   - Outer integral (x): 32-point Gauss-Legendre quadrature
//     over the finite segment width
//
// The core component (R_C ~ 0.05-0.3 RM) is too narrow for
// GL quadrature in x, so all core energy is assigned to the
// hit's segment. The tail component (R_T ~ 0.4-2.4 RM) is
// smooth enough for 32-point GL over the finite x-range.
//
// Reference: Grindhammer & Peters, arXiv:hep-ex/0001020
///////////////////////////////////////////////////////////////

namespace SegIntegral {

    //--- Analytic y-integral for one component ---
    // Evaluates the definite integral of
    //   R^2 / (pi * (x^2 + y^2 + R^2)^2)
    // over y from y1 to y2, at a fixed x.
    // All coordinates in RM units.
    inline Double_t AnalyticYIntegral(Double_t x, Double_t y1, Double_t y2,
                                      Double_t R) {
        Double_t a2 = x * x + R * R;
        Double_t a  = sqrt(a2);
        Double_t R2 = R * R;
        // Antiderivative: (R^2/pi) * [y/(2*a^2*(a^2+y^2)) + atan(y/a)/(2*a^3)]
        Double_t F2 = (R2 / TMath::Pi()) *
            (y2 / (2.0 * a2 * (a2 + y2 * y2)) +
             atan2(y2, a) / (2.0 * a * a2));
        Double_t F1 = (R2 / TMath::Pi()) *
            (y1 / (2.0 * a2 * (a2 + y1 * y1)) +
             atan2(y1, a) / (2.0 * a * a2));
        return F2 - F1;
    }

    //--- 32-point Gauss-Legendre nodes and weights ---
    const Int_t nGL = 32;
    const Double_t glX[32] = {
        -0.9972638618494816, -0.9856115115452684,
        -0.9647622555875064, -0.9349060759377397,
        -0.8963211557660521, -0.8493676137325700,
        -0.7944837959679424, -0.7321821187402897,
        -0.6630442669302152, -0.5877157572407623,
        -0.5068999089322294, -0.4213512761306353,
        -0.3318686022821276, -0.2392873622521371,
        -0.1444719615827965, -0.0483076656877383,
         0.0483076656877383,  0.1444719615827965,
         0.2392873622521371,  0.3318686022821276,
         0.4213512761306353,  0.5068999089322294,
         0.5877157572407623,  0.6630442669302152,
         0.7321821187402897,  0.7944837959679424,
         0.8493676137325700,  0.8963211557660521,
         0.9349060759377397,  0.9647622555875064,
         0.9856115115452684,  0.9972638618494816
    };
    const Double_t glW[32] = {
         0.0070186100094701,  0.0162743947309057,
         0.0253920653092621,  0.0342738629130214,
         0.0428358980222267,  0.0509980592623762,
         0.0586840934785355,  0.0658222227763618,
         0.0723457941088485,  0.0781938957870703,
         0.0833119242269467,  0.0876520930044038,
         0.0911738786957639,  0.0938443990808046,
         0.0956387200792749,  0.0965400885147278,
         0.0965400885147278,  0.0956387200792749,
         0.0938443990808046,  0.0911738786957639,
         0.0876520930044038,  0.0833119242269467,
         0.0781938957870703,  0.0723457941088485,
         0.0658222227763618,  0.0586840934785355,
         0.0509980592623762,  0.0428358980222267,
         0.0342738629130214,  0.0253920653092621,
         0.0162743947309057,  0.0070186100094701
    };

    //--- Semi-analytic integration over a finite rectangle ---
    // Integrates the density for one component over [x1,x2] × [y1,y2]
    // in hit-relative RM units. Analytic in y, 32-pt GL in x.
    inline Double_t IntegrateSegment(Double_t x1, Double_t x2,
                                     Double_t y1, Double_t y2,
                                     Double_t R) {
        Double_t mx = 0.5 * (x2 - x1);
        Double_t cx = 0.5 * (x2 + x1);
        Double_t sum = 0.0;
        for (Int_t i = 0; i < nGL; i++) {
            Double_t xi = mx * glX[i] + cx;
            sum += glW[i] * AnalyticYIntegral(xi, y1, y2, R);
        }
        return sum * mx;
    }

} // namespace SegIntegral

///////////////////////////////////////////////////////////////
// ComputePmtFractions()
// For a single hit entering at (hitXcm, hitYcm) on the
// calorimeter face, compute the fractional energy deposited
// in each of the 8 PMT segments at a given shower depth.
//
// Uses the finite block dimensions (±blockHalfX in x, row
// boundaries in y). Energy falling outside the block is not
// counted (lateral and vertical leakage).
//
// Strategy:
//   - Core: all energy assigned to the hit's segment (core is
//     narrow enough that <1% crosses any boundary)
//   - Tail: semi-analytic integration over each segment's
//     finite rectangle (analytic in y, 32-pt GL in x)
//
// hitXcm, hitYcm: entry point in cm
// p, Rc, Rt: Grindhammer-Peters radial parameters (in RM units)
// frac[8]: output array of energy fractions per PMT segment
///////////////////////////////////////////////////////////////

inline void ComputePmtFractions(Double_t hitXcm, Double_t hitYcm,
                                Double_t p, Double_t Rc, Double_t Rt,
                                Double_t frac[8]) {

    const Double_t RM = PbConst::RM;

    // Hit position in RM units
    Double_t hx = hitXcm / RM;
    Double_t hy = hitYcm / RM;

    // Block x-boundaries in hit-relative RM units
    Double_t blockXmin = -DetGeom::blockHalfX / RM - hx;  // left edge of block
    Double_t blockXmax =  DetGeom::blockHalfX / RM - hx;  // right edge of block

    // L/R boundary in hit-relative RM units
    Double_t lrBound = -hx;  // x=0 in detector coords

    // Determine which column the hit is in
    Int_t hitCol = (hitXcm >= 0.0) ? 0 : 1;

    for (Int_t row = 0; row < 4; row++) {
        // Hit-relative y boundaries in RM units
        Double_t ry1 = DetGeom::rowYBound[row + 1] / RM - hy;  // lower
        Double_t ry2 = DetGeom::rowYBound[row] / RM - hy;      // upper

        // --- Core component: assign all to hit's segment ---
        // Core is so narrow (<0.3 RM) that essentially all energy
        // stays in the hit's segment. Compute the total core
        // fraction in this row using the analytic y-integral
        // evaluated at x=0 (hit location), integrated over a
        // narrow x-range. For simplicity, use the infinite-strip
        // result since core leakage beyond ±9 cm is negligible.
        Double_t coreRow = ry2 / (2.0 * sqrt(ry2 * ry2 + Rc * Rc))
                         - ry1 / (2.0 * sqrt(ry1 * ry1 + Rc * Rc));
        Double_t coreLeft  = (hitCol == 0) ? coreRow : 0.0;
        Double_t coreRight = (hitCol == 1) ? coreRow : 0.0;

        // --- Tail component: semi-analytic over finite segments ---
        // Left column: x ∈ [max(blockXmin, lrBound), blockXmax]
        //   but left means x > 0 in detector coords, so x > lrBound
        // Right column: x ∈ [blockXmin, min(blockXmax, lrBound)]
        Double_t leftX1  = (lrBound > blockXmin) ? lrBound : blockXmin;
        Double_t leftX2  = blockXmax;
        Double_t rightX1 = blockXmin;
        Double_t rightX2 = (lrBound < blockXmax) ? lrBound : blockXmax;

        Double_t tailLeft  = 0.0;
        Double_t tailRight = 0.0;

        if (leftX2 > leftX1)
            tailLeft = SegIntegral::IntegrateSegment(leftX1, leftX2,
                                                     ry1, ry2, Rt);
        if (rightX2 > rightX1)
            tailRight = SegIntegral::IntegrateSegment(rightX1, rightX2,
                                                      ry1, ry2, Rt);

        // Clamp to non-negative
        if (tailLeft  < 0.0) tailLeft  = 0.0;
        if (tailRight < 0.0) tailRight = 0.0;

        // --- Combine core and tail ---
        frac[0 * 4 + row] = p * coreLeft  + (1.0 - p) * tailLeft;   // Left
        frac[1 * 4 + row] = p * coreRight + (1.0 - p) * tailRight;  // Right
    }
}

///////////////////////////////////////////////////////////////
// Shower processing — two-stage approach
//
// Stage 1: ComputeShowerProfile()
//   Computes the full shower shape (slice energies, PMT fractions,
//   radial parameters) at each depth step. This is independent of
//   where the shower starts in the block.
//
// Stage 2: AccumulateShowerInBlock()
//   Shifts the profile by a start depth offset and accumulates
//   energy into pmtE[8], stopping at the block back face
//   (longitudinal leakage).
//
// The separation allows multiple MC throws of fiber penetration
// depth to reuse a single shower profile computation.
///////////////////////////////////////////////////////////////

const Int_t MAX_DEPTH_STEPS = 320;
const Int_t DR_MAXNHIT      = 20;    // max EM hits per event on det 9
const Int_t DR_NTHROWS       = 20;   // MC throws for fiber penetration

// Minimum energy for shower development [GeV].
// Below this, particles are deposited as point-like.
const Double_t MIN_SHOWER_E_GEV = 0.050;  // 50 MeV
const Double_t SLICE_X0         = 0.25;   // slice thickness [X0]

///////////////////////////////////////////////////////////////
// ShowerProfile — pre-computed shower shape for a single hit.
///////////////////////////////////////////////////////////////

struct ShowerProfile {
    Int_t    nSteps;                       // valid depth steps
    Double_t hitEGeV;                      // hit energy [GeV]
    Double_t hitXcm, hitYcm;              // hit entry point [cm]
    Bool_t   belowThreshold;               // true if E < minShowerE
    Int_t    pointSeg;                     // segment for below-threshold
    Double_t sliceFrac[MAX_DEPTH_STEPS];   // fractional energy per slice
    Double_t frac[MAX_DEPTH_STEPS * 8];    // PMT fractions [step*8+seg]
    Double_t depthT[MAX_DEPTH_STEPS];      // depth midpoint [X0]
    Double_t rc[MAX_DEPTH_STEPS];          // core radius [RM]
    Double_t rt[MAX_DEPTH_STEPS];          // tail radius [RM]
    Double_t pw[MAX_DEPTH_STEPS];          // core weight
};

///////////////////////////////////////////////////////////////
// ComputeShowerProfile()
///////////////////////////////////////////////////////////////

inline void ComputeShowerProfile(Double_t hitEGeV, Double_t hitXcm,
                                  Double_t hitYcm, ShowerProfile &prof) {
    prof.hitEGeV = hitEGeV;
    prof.hitXcm  = hitXcm;
    prof.hitYcm  = hitYcm;
    prof.nSteps  = 0;
    prof.belowThreshold = false;

    Double_t hitEMeV = hitEGeV * 1000.0;
    Double_t y = hitEMeV / PbConst::Ec;

    if (y <= 1.0 || hitEGeV < MIN_SHOWER_E_GEV) {
        prof.belowThreshold = true;
        Int_t col = (hitXcm >= 0.0) ? 0 : 1;
        Int_t row = 3;
        for (Int_t r = 0; r < 4; r++) {
            if (hitYcm >= DetGeom::rowYBound[r + 1]) { row = r; break; }
        }
        prof.pointSeg = col * 4 + row;
        return;
    }

    Double_t lny = TMath::Log(y);
    Double_t lnE = TMath::Log(hitEMeV);

    Double_t T     = GPShower::Thom(lny);
    Double_t alpha = GPShower::AlphaHom(lny, PbConst::Z);
    if (alpha <= 1.0) alpha = 1.01;
    Double_t beta  = (alpha - 1.0) / T;

    // Full shower extent (block depth cap applied during accumulation)
    Double_t tMax = std::max(3.5 * T, 20.0);
    if (tMax > (Double_t)MAX_DEPTH_STEPS * SLICE_X0)
        tMax = (Double_t)MAX_DEPTH_STEPS * SLICE_X0;
    Int_t nSteps = (Int_t)std::ceil(tMax / SLICE_X0);
    if (nSteps > MAX_DEPTH_STEPS) nSteps = MAX_DEPTH_STEPS;

    Double_t cumFrac = 0.0;
    Int_t nValid = 0;

    for (Int_t iStep = 0; iStep < nSteps; iStep++) {
        Double_t t1 = iStep * SLICE_X0;
        Double_t t2 = t1 + SLICE_X0;
        Double_t tMid = t1 + SLICE_X0 / 2.0;

        Double_t sf = GPShower::GammaSliceIntegral(t1, t2, alpha, beta);
        cumFrac += sf;

        if (sf < 1.0e-8) {
            if (cumFrac > 0.999) break;
            continue;
        }

        Double_t tau = tMid / T;
        Double_t rc  = GPShower::RC(tau, lnE, PbConst::Z);
        Double_t rt  = GPShower::RT(tau, lnE, PbConst::Z);
        Double_t p   = GPShower::CoreWeight(tau, lnE, PbConst::Z);

        Double_t stepFrac[8];
        ComputePmtFractions(hitXcm, hitYcm, p, rc, rt, stepFrac);

        prof.sliceFrac[nValid] = sf;
        prof.depthT[nValid]    = tMid;
        prof.rc[nValid]        = rc;
        prof.rt[nValid]        = rt;
        prof.pw[nValid]        = p;
        for (Int_t s = 0; s < 8; s++)
            prof.frac[nValid * 8 + s] = stepFrac[s];

        nValid++;
        if (cumFrac > 0.999) break;
    }

    prof.nSteps = nValid;
}

///////////////////////////////////////////////////////////////
// AccumulateShowerInBlock()
// Shifts profile by startDepth_cm and accumulates into pmtE[8].
// Slices beyond the calorimeter depth are skipped (leakage).
///////////////////////////////////////////////////////////////

inline void AccumulateShowerInBlock(const ShowerProfile &prof,
                                     Double_t startDepth_cm,
                                     Double_t pmtE[8]) {
    if (prof.belowThreshold) {
        if (startDepth_cm < DetGeom::caloDepth &&
            prof.hitYcm >= DetGeom::rowYBound[4] &&
            prof.hitYcm <= DetGeom::rowYBound[0]) {
            pmtE[prof.pointSeg] += prof.hitEGeV;
        }
        return;
    }

    Double_t startDepthX0 = startDepth_cm / PbConst::X0;

    for (Int_t iStep = 0; iStep < prof.nSteps; iStep++) {
        // Slice boundaries in X0 from shower start
        Double_t sliceStartX0 = prof.depthT[iStep] - SLICE_X0 / 2.0;
        Double_t sliceEndX0   = prof.depthT[iStep] + SLICE_X0 / 2.0;

        // Physical position in block [X0]
        Double_t physStartX0 = startDepthX0 + sliceStartX0;
        Double_t physEndX0   = startDepthX0 + sliceEndX0;

        if (physStartX0 >= DetGeom::caloDepthX0) break;  // entire slice beyond block

        // Fractional last slice: if slice partially extends beyond block
        Double_t fracUsed = 1.0;
        if (physEndX0 > DetGeom::caloDepthX0) {
            fracUsed = (DetGeom::caloDepthX0 - physStartX0) / SLICE_X0;
        }

        Double_t sliceE = prof.hitEGeV * prof.sliceFrac[iStep] * fracUsed;
        for (Int_t s = 0; s < 8; s++)
            pmtE[s] += sliceE * prof.frac[iStep * 8 + s];

        if (fracUsed < 1.0) break;  // partial slice was the last one
    }
}

///////////////////////////////////////////////////////////////
// DetermineStartMaterial()
// Flat MC to determine if the particle enters lead or fiber.
// Returns true if fiber, false if lead.
///////////////////////////////////////////////////////////////

inline Bool_t DetermineStartMaterial(TRandom3 *rng) {
    return (rng->Uniform() < DetGeom::fiberVolFrac);
}

///////////////////////////////////////////////////////////////
// SampleShowerStartDepth()
// Samples the depth at which the shower initiates, drawn from
// an exponential distribution with mean = X0_start.
//
// PDF:  P(z) = (1/X0) * exp(-z/X0)
// CDF:  F(z) = 1 - exp(-z/X0)
// Inverse CDF: z = -X0 * ln(1-u), u ~ Uniform(0,1)
///////////////////////////////////////////////////////////////

inline Double_t SampleShowerStartDepth(Double_t X0_start, TRandom3 *rng) {
    return rng->Exp(X0_start);
}

///////////////////////////////////////////////////////////////
// SampleFiberTransverseDistance()
// Samples the transverse distance (chord length) a particle
// traverses through a fiber, using the inverse CDF method.
//
// CDF: f(u) = 1 - sqrt(1 - u^2),  u ~ Uniform(0,1)
// Returns a chord length as a fraction of the fiber diameter,
// scaled to cm.
///////////////////////////////////////////////////////////////

inline Double_t SampleFiberTransverseDistance(TRandom3 *rng) {
    Double_t u = rng->Uniform();
    Double_t frac = sqrt(1.0 - u * u);
    return frac * DetGeom::fiberDiam;
}

///////////////////////////////////////////////////////////////
// ComputeFiberPenetration()
// Monte Carlo determination of the shower start depth [cm]
// in the calorimeter block, accounting for fiber penetration.
//
// Process:
//   1. Determine if the particle enters lead or fiber.
//   2. Sample z_start from exponential(X0_start).
//   3. If fiber: compute z_available from the incidence angle
//      and fiber transverse distance. If z_start > z_available,
//      the particle exits the fiber at z_available and enters
//      lead (effective start depth = z_available).
//   4. If lead: z_start is the effective start depth directly.
//
// Returns the effective shower start depth [cm] in the block.
// Returns -1.0 if the particle passes through the entire block
// without showering (complete leakage).
// If z_avail_out is non-null, stores the available fiber depth
// (only meaningful when the particle hits a fiber).
// If z_start_raw_out is non-null, stores the raw exponential
// draw before clamping to z_available (for diagnostics).
// If d_trans_out is non-null, stores the sampled fiber transverse
// distance (only meaningful when the particle hits a fiber).
// If verbose is true, prints detailed MC calculation steps.
///////////////////////////////////////////////////////////////

inline Double_t ComputeFiberPenetration(Double_t hitPx, Double_t hitPy,
                                         Double_t hitPz, TRandom3 *rng,
                                         Double_t *z_avail_out = 0,
                                         Double_t *z_start_raw_out = 0,
                                         Double_t *d_trans_out = 0,
                                         Bool_t verbose = false) {

    if (z_avail_out) *z_avail_out = -1.0;
    if (z_start_raw_out) *z_start_raw_out = -1.0;
    if (d_trans_out) *d_trans_out = -1.0;

    Bool_t isFiber = DetermineStartMaterial(rng);
    Double_t X0_start = isFiber ? DetGeom::fiberX0 : PbConst::X0;

    // Sample shower start depth from exponential distribution
    Double_t z_start = SampleShowerStartDepth(X0_start, rng);
    if (z_start_raw_out) *z_start_raw_out = z_start;

    if (verbose) {
        printf("      Material: %s (X0 = %.4f cm)\n",
               isFiber ? "FIBER" : "LEAD", X0_start);
        printf("      z_start (raw): %.4f cm\n", z_start);
        printf("      hitP: (%.4f, %.4f, %.4f) GeV\n", hitPx, hitPy, hitPz);
    }

    if (isFiber) {
        // Compute transverse distance available in the fiber
        Double_t d_trans = SampleFiberTransverseDistance(rng);

        if (d_trans_out) *d_trans_out = d_trans;

        // Compute incidence angle relative to fiber axis.
        Double_t pT = sqrt(hitPx * hitPx + hitPy * hitPy);
        Double_t theta_particle = atan2(pT, hitPz);
        Double_t theta_incident = fabs(theta_particle - DetGeom::detAngleRad);

        if (verbose) {
            printf("      d_trans: %.4f cm\n", d_trans);
            printf("      theta_particle: %.4f rad (%.2f deg)\n",
                   theta_particle, theta_particle * 180.0 / 3.14159265358979);
            printf("      theta_incident: %.6f rad (%.4f deg)\n",
                   theta_incident, theta_incident * 180.0 / 3.14159265358979);
        }

        // Protect against near-zero incidence (nearly parallel to fiber)
        if (theta_incident < 1.0e-6) theta_incident = 1.0e-6;

        // Available depth in the fiber before exiting into lead
        Double_t z_available = d_trans / sin(theta_incident);

        if (z_avail_out) *z_avail_out = z_available;

        if (verbose) {
            printf("      z_available: %.4f cm\n", z_available);
            printf("      z_start %s z_available => %s\n",
                   z_start > z_available ? ">" : "<=",
                   z_start > z_available ? "exits fiber, clamped" : "showers in fiber");
        }

        // If z_start exceeds z_available, particle exits fiber
        // at z_available and enters lead
        if (z_start > z_available) {
            z_start = z_available;
        }
    }

    // Check for complete pass-through (beyond block depth)
    if (z_start > DetGeom::caloDepth) {
        return -1.0;
    }

    if (verbose) {
        printf("      effective start depth: %.4f cm\n", z_start);
    }

    return z_start;
}

///////////////////////////////////////////////////////////////
// FillShowerDiag()
// Fills the TShower diagnostic tree from a ShowerProfile.
///////////////////////////////////////////////////////////////

Int_t   diag_eventIdx;
Int_t   diag_hitIdx;
Float_t diag_hitX;
Float_t diag_hitY;
Float_t diag_hitE;
Int_t   diag_nDepthSteps;
Float_t diag_depthT[MAX_DEPTH_STEPS];
Float_t diag_sliceE[MAX_DEPTH_STEPS];
Float_t diag_RC[MAX_DEPTH_STEPS];
Float_t diag_RT[MAX_DEPTH_STEPS];
Float_t diag_p[MAX_DEPTH_STEPS];
Float_t diag_pmtFrac[MAX_DEPTH_STEPS * 8];
Float_t diag_hitZ[DR_NTHROWS];  // shower start depth per throw [cm]

inline void FillShowerDiag(const ShowerProfile &prof,
                            TTree *diagTree, Int_t eventIdx, Int_t hitIdx,
                            const Float_t hitZ[DR_NTHROWS] = 0) {
    diag_eventIdx = eventIdx;
    diag_hitIdx   = hitIdx;
    diag_hitX     = (Float_t)prof.hitXcm;
    diag_hitY     = (Float_t)prof.hitYcm;
    diag_hitE     = (Float_t)prof.hitEGeV;
    diag_nDepthSteps = prof.nSteps;

    for (Int_t k = 0; k < MAX_DEPTH_STEPS; k++) {
        diag_depthT[k] = 0.0f;
        diag_sliceE[k] = 0.0f;
        diag_RC[k]     = 0.0f;
        diag_RT[k]     = 0.0f;
        diag_p[k]      = 0.0f;
    }
    for (Int_t k = 0; k < MAX_DEPTH_STEPS * 8; k++)
        diag_pmtFrac[k] = 0.0f;

    for (Int_t iStep = 0; iStep < prof.nSteps; iStep++) {
        diag_depthT[iStep] = (Float_t)prof.depthT[iStep];
        diag_sliceE[iStep] = (Float_t)(prof.hitEGeV * prof.sliceFrac[iStep]);
        diag_RC[iStep]     = (Float_t)prof.rc[iStep];
        diag_RT[iStep]     = (Float_t)prof.rt[iStep];
        diag_p[iStep]      = (Float_t)prof.pw[iStep];
        for (Int_t s = 0; s < 8; s++)
            diag_pmtFrac[iStep * 8 + s] = (Float_t)prof.frac[iStep * 8 + s];
    }

    // Copy hitZ (start depths per throw)
    for (Int_t k = 0; k < DR_NTHROWS; k++)
        diag_hitZ[k] = (hitZ != 0) ? hitZ[k] : 0.0f;

    diagTree->Fill();
}


///////////////////////////////////////////////////////////////
// DetResp output tree branch buffers and reader functions
// (for use by downstream plotting/analysis scripts)
///////////////////////////////////////////////////////////////

// Branch buffers for TDetResp tree
Double_t dr_evXs, dr_evAsym, dr_evUnpolWght;
Double_t dr_evPolPlusWghtZ;
Double_t dr_evPolMinusWghtZ;
Float_t  dr_pmtE[DR_NTHROWS * 8];       // [throw*8 + seg]
Float_t  dr_pmtETotal[DR_NTHROWS];
Float_t  dr_pmtELeft[DR_NTHROWS];
Float_t  dr_pmtERight[DR_NTHROWS];
Float_t  dr_hitZ[DR_NTHROWS];           // shower start depth per throw [cm]
Double_t dr_hitPx[DR_MAXNHIT];
Double_t dr_hitPy[DR_MAXNHIT];
Double_t dr_hitPz[DR_MAXNHIT];
Double_t dr_hitE[DR_MAXNHIT];
Int_t    dr_nHitsDet9;

// Branch pointers
TBranch *b_dr_evXs, *b_dr_evAsym, *b_dr_evUnpolWght;
TBranch *b_dr_evPolPlusWghtZ;
TBranch *b_dr_evPolMinusWghtZ;
TBranch *b_dr_pmtE, *b_dr_pmtETotal, *b_dr_pmtELeft, *b_dr_pmtERight;
TBranch *b_dr_hitZ;
TBranch *b_dr_hitPx, *b_dr_hitPy, *b_dr_hitPz, *b_dr_hitE;
TBranch *b_dr_nHitsDet9;

///////////////////////////////////////////////////////////////
// SetupDetRespChain()
///////////////////////////////////////////////////////////////

inline TChain *SetupDetRespChain(const char *fileList) {
    TChain *chain = new TChain("TDetResp");
    std::string input(fileList);

    if (input.find(".root") != std::string::npos) {
        chain->Add(fileList);
        printf("DetResp: added single file: %s\n", fileList);
    } else {
        std::ifstream ifile(fileList);
        if (!ifile.is_open()) {
            printf("Error: cannot open file list %s\n", fileList);
            delete chain;
            return 0;
        }
        std::string line;
        Int_t nFiles = 0;
        while (ifile >> line) {
            chain->Add(line.c_str());
            nFiles++;
        }
        ifile.close();
        printf("DetResp: added %d files from %s\n", nFiles, fileList);
    }

    if (chain->GetEntries() == 0) {
        printf("Error: no entries found in TDetResp chain\n");
        delete chain;
        return 0;
    }
    return chain;
}

///////////////////////////////////////////////////////////////
// SetupDetRespBranches()
///////////////////////////////////////////////////////////////

inline void SetupDetRespBranches(TTree *tree) {
    tree->SetBranchAddress("evXs",            &dr_evXs,            &b_dr_evXs);
    tree->SetBranchAddress("evAsym",          &dr_evAsym,          &b_dr_evAsym);
    tree->SetBranchAddress("evUnpolWght",     &dr_evUnpolWght,     &b_dr_evUnpolWght);
    tree->SetBranchAddress("evPolPlusWghtZ",  &dr_evPolPlusWghtZ,  &b_dr_evPolPlusWghtZ);
    tree->SetBranchAddress("evPolMinusWghtZ", &dr_evPolMinusWghtZ, &b_dr_evPolMinusWghtZ);
    tree->SetBranchAddress("pmtE",            dr_pmtE,             &b_dr_pmtE);
    tree->SetBranchAddress("pmtETotal",       dr_pmtETotal,        &b_dr_pmtETotal);
    tree->SetBranchAddress("pmtELeft",        dr_pmtELeft,         &b_dr_pmtELeft);
    tree->SetBranchAddress("pmtERight",       dr_pmtERight,        &b_dr_pmtERight);
    tree->SetBranchAddress("hitZ",            dr_hitZ,             &b_dr_hitZ);
    tree->SetBranchAddress("hitPx",           dr_hitPx,            &b_dr_hitPx);
    tree->SetBranchAddress("hitPy",           dr_hitPy,            &b_dr_hitPy);
    tree->SetBranchAddress("hitPz",           dr_hitPz,            &b_dr_hitPz);
    tree->SetBranchAddress("hitE",            dr_hitE,             &b_dr_hitE);
    tree->SetBranchAddress("nHitsDet9",       &dr_nHitsDet9,       &b_dr_nHitsDet9);
}

///////////////////////////////////////////////////////////////
// DetResp event asymmetry (Z-polarization)
///////////////////////////////////////////////////////////////

inline Double_t DetRespEventAsymmetry() {
    Double_t numer = dr_evPolPlusWghtZ - dr_evPolMinusWghtZ;
    Double_t denom = dr_evPolPlusWghtZ + dr_evPolMinusWghtZ;
    if (denom == 0.0) return 0.0;
    return numer / denom;
}

#endif // MOLPOL_DETECTOR_RESPONSE_H