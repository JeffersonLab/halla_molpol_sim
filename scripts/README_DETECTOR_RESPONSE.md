# MolPol Detector Response Model

Parameterized electromagnetic shower model for the MolPol lead-fiber calorimeter, applied to flux-plane (detector 9) simulation data. Replaces the "perfect detector" approximation with realistic energy deposition across the 2×4 PMT segment grid, enabling studies of shower leakage effects on measured asymmetry and rate as a function of active PMT configuration.

## Shower Parameterization

The shower model follows the Grindhammer-Peters parameterization for homogeneous media (Appendix A.1):

> G. Grindhammer and S. Peters,
> "The Parameterized Simulation of Electromagnetic Showers in Homogeneous and Sampling Calorimeters,"
> arXiv:hep-ex/0001020 (2000).
> [https://arxiv.org/abs/hep-ex/0001020](https://arxiv.org/abs/hep-ex/0001020)

The longitudinal shower profile uses a gamma distribution (Eq. 2) with parameters T and α determined by the incident energy (Eqs. 7–8). The lateral profile at each depth slice uses a two-component radial distribution (Eq. 23) with depth-dependent core radius R_C (Eq. 24), tail radius R_T (Eq. 25), and core weight p (Eq. 26). The two components are a mathematical decomposition (narrow + broad terms) chosen to reproduce the radial profile shape — they do not represent physically distinct shower processes. R_C and R_T are the medians of each component, and p controls the relative weight. Each component has the form `f(r) = 2rR² / (r² + R²)²`, which integrates to the CDF `F(r) = r² / (r² + R²)` (Eq. 27). Note: the paper uses unsubscripted R in Eqs. 27–28 as a generic stand-in for either R_C or R_T; this is not explicitly stated.

Material constants use pure lead values (precomputed, NIST):

| Constant | Value | Description |
|----------|-------|-------------|
| Z | 82.0 | Atomic number |
| A | 207.2 | Atomic mass |
| ρ | 11.35 g/cm³ | Density |
| X₀ | 0.5612 cm (6.37 g/cm²) | Radiation length |
| E_c | 7.3548 MeV | Critical energy (Eq. 5, using X₀ in g/cm²) |
| R_M | 1.6176 cm | Molière radius |
| E_s | 21.2 MeV | Scale energy |
| Depth | 30.0 cm (53.5 X₀) | Calorimeter longitudinal depth |

### Segment Integration Method

Energy fractions per PMT segment are computed using a hybrid analytic/numerical method:

1. **Row fractions (y-direction)**: Fully analytic. Each row is treated as an infinite strip in x. The density integrates analytically over all x, and the resulting y-integral has the closed form: `F_strip = y₂/(2√(y₂²+R²)) - y₁/(2√(y₁²+R²))`.

2. **L/R split (x-direction)**: The left-side fraction is decomposed as `F_left = F_symmetric + F_correction`, where the symmetric half is analytic (`y/(4√(y²+R²))`) and the shift correction accounts for the hit's displacement from the L/R boundary using 32-point Gauss-Legendre quadrature in y.

3. **Core/tail separation**: The core component (R_C ~ 0.05–0.3 R_M) is too narrow for the GL correction to resolve accurately, so all core energy is assigned to the hit's column. The tail component (R_T ~ 0.4–2.4 R_M) is smooth enough for the 32-point GL.

The infinite-strip approximation in x is valid because 99% of shower energy is contained within ~4.6 cm of the hit, while the detector face is 18 cm wide. Even for hits at the detector opening edge (±6.3 cm), less than 1% of shower energy leaks off the detector face laterally.

### Assumptions and Limitations

- **EM particles only**: Only e⁻ (PID 11), e⁺ (PID -11), and γ (PID 22) are processed. All other particles are excluded (see Things Encountered below).
- **Minimum shower energy**: Particles below 50 MeV are deposited as point-like in the containing segment rather than showered. Below this energy, the Grindhammer-Peters parameterization produces unphysical radial parameters (see Things Encountered below).
- **Normal incidence**: All particles are assumed to enter the calorimeter face head-on. The shower axis is not tilted along the hit momentum vector.
- **Pure lead**: The Molière radius uses homogeneous Pb values rather than an effective lead-fiber sampling value.
- **Linear superposition**: Multiple hits within an event are showered independently and their PMT deposits are summed.

### Longitudinal Leakage

The calorimeter has a physical depth of 30 cm, corresponding to 53.5 X₀ in pure lead. The shower depth integration is capped at this value — any shower energy developing beyond the back face of the block is lost. This models longitudinal leakage physically: the gamma distribution determines how much energy would develop at each depth, and slices beyond the calorimeter depth are simply not accumulated. The calorimeter depth in radiation lengths is computed from `DetGeom::caloDepth / PbConst::X0`, so changing the material or depth automatically adjusts the leakage.

### Energy Resolution Smearing

After shower energy is distributed across the 8 PMT segments, each segment's energy is independently smeared with a Gaussian to model sampling fluctuations and photostatistics. The resolution parameterization follows:

> R. Livan, M. Vercesi, R. Wigmans,
> "Scintillating-Fibre Calorimetry" (1995).

$$\sigma/E = 12.9\%/\sqrt{E} + 1.2\%$$

where E is in GeV. Representative values:

| E (GeV) | σ/E | σ (MeV) |
|---------|-----|---------|
| 1.0 | 14.1% | 141 |
| 2.0 | 10.3% | 206 |
| 6.5 | 6.3% | 407 |
| 11.0 | 5.1% | 560 |

The smearing is applied per-PMT (not per-hit) because the physical source of the resolution — sampling fluctuations in the lead-fiber stack and photostatistics at each PMT — acts independently on each segment's collected light. Negative values after smearing are clamped to zero. The energy conservation checks run before smearing, so they validate the shower model itself without confusion from the stochastic component. Arm sums (`pmtELeft`, `pmtERight`, `pmtETotal`) are recomputed from the smeared segment values.

The smearing uses ROOT's `TRandom3` generator (Mersenne Twister).

## Files

| File | Description |
|------|-------------|
| `MolPolDetectorResponse.h` | Shower model implementation: material constants (`PbConst`), detector geometry (`DetGeom`), Grindhammer-Peters parameterization (`GPShower`), segment integration (`SegIntegral`), `ComputePmtFractions()`, `ProcessHitShower()`. Also provides `SetupDetRespChain()` / `SetupDetRespBranches()` for reading output files, and `DetRespEventAsymmetry()`. |
| `MolPolDetectorResponse.C` | ROOT macro driver: beam energy determination, event loop with shower processing, output tree creation, optional diagnostic tree (capped at 10,000 events), and energy conservation checks. |
| `DetectorResponsePlots.C` | Plotting macro for detector response output files. Gate mode (coincidence/left/right) and energy threshold selection. Produces PMT segment energy spectra, left/right/total energy sums, left-vs-right heatmap, and interactive 3D shower visualization. |

### Dependencies

- `MolPolData.h` — branch buffer definitions and `SetupMolpolChain()` / `SetupMolpolBranches()`.
- `MolPolAnalysis.h` — `hitFluxTrackCoinc()` for beam energy determination.
- `MolPolDetectorResponse.h` — required by both `.C` scripts.
- ROOT (TTree, TChain, TFile, TMath, TGraph2D, TPolyLine3D, TCanvas, TH1F, TH2F).

## Usage

### Running the Detector Response Model

```bash
# Basic run
root 'MolPolDetectorResponse.C("input.root")'

# With shower diagnostic output (limited to 10,000 events)
root 'MolPolDetectorResponse.C("input.root", true)'

# From a file list
root 'MolPolDetectorResponse.C("filelist.txt")'
```

The output file is named `<input>_detresponse.root` (or `detresponse_output.root` if the input is a file list).

### Plotting

```bash
# Default: coincidence mode, no energy cut
root 'DetectorResponsePlots.C("input_detresponse.root")'

# Coincidence with 0.5 GeV cut
root 'DetectorResponsePlots.C("input_detresponse.root", "coincidence", 0.5)'

# Left singles with 2.0 GeV cut
root 'DetectorResponsePlots.C("input_detresponse.root", "left", 2.0)'

# Right singles with 1.0 GeV cut (short form)
root 'DetectorResponsePlots.C("input_detresponse.root", "r", 1.0)'
```

Gate modes (case-insensitive): `"coincidence"` / `"c"`, `"left"` / `"l"`, `"right"` / `"r"`. Partial matches on first letter are accepted with a warning.

Gate definitions: coincidence requires `pmtELeft >= cut AND pmtERight >= cut`. Left/right requires the respective arm total to exceed the cut. Individual PMT segment values within passing events may be below the cut threshold (energy is distributed across segments by the shower model).

#### Shower Visualization (interactive)

```cpp
root [0] .L DetectorResponsePlots.C
root [1] TChain *t = new TChain("TShower");
root [2] t->Add("input_detresponse.root");
root [3] PlotEventShower(t, 1);   // 1st unique event (1-based indexing)
root [4] PlotEventShower(t, 42);  // 42nd unique event
```

#### Plot Canvases

**cPmtSeg** (4 rows × 2 columns) — Energy spectrum for each of the 8 PMT segments. Canvas layout mirrors the physical detector face: left column = L1–L4, right column = R1–R4, top to bottom.

**cPmtSums** (3 rows × 1 column) — Energy spectra for `pmtELeft` (L1+L2+L3+L4), `pmtERight` (R1+R2+R3+R4), and `pmtETotal` (all 8 segments).

**cPmtLR** (single pad) — 2D heatmap (COLZ) of `pmtELeft` vs `pmtERight`, showing left-right arm correlation.

**cShower_evN** (single pad, 3D) — Shower development for event N. Each depth slice is drawn as a TGraph2D circle at (hitX, hitY, depth). Disc radius represents the 90% energy containment radius computed from the two-component CDF. Marker size scales with slice energy. Multiple hits in the same event are distinguished by color. PMT segment grid boundaries are drawn at z = 0.

## Output ROOT File

### TDetResp Tree (per event)

Event-level physics branches, copied through unchanged from the input:

| Branch | Type | Description |
|--------|------|-------------|
| `evXs` | `Double_t` | Moller cross section for this event |
| `evAsym` | `Double_t` | Event asymmetry |
| `evUnpolWght` | `Double_t` | Unpolarized event weight |
| `evPolPlusWghtX` | `Double_t` | Positive helicity weight, X polarization |
| `evPolPlusWghtY` | `Double_t` | Positive helicity weight, Y polarization |
| `evPolPlusWghtZ` | `Double_t` | Positive helicity weight, Z (longitudinal) polarization |
| `evPolMinusWghtX` | `Double_t` | Negative helicity weight, X polarization |
| `evPolMinusWghtY` | `Double_t` | Negative helicity weight, Y polarization |
| `evPolMinusWghtZ` | `Double_t` | Negative helicity weight, Z (longitudinal) polarization |

Detector response branches:

| Branch | Type | Description |
|--------|------|-------------|
| `pmtE[8]` | `Double_t[8]` | Deposited energy per PMT segment [GeV]. Summed over all EM hits and all shower depth slices. See segment indexing below. |
| `pmtETotal` | `Double_t` | Sum of all 8 segments — total detected energy [GeV]. |
| `pmtELeft` | `Double_t` | Sum of L1+L2+L3+L4 (segments 0–3) [GeV]. |
| `pmtERight` | `Double_t` | Sum of R1+R2+R3+R4 (segments 4–7) [GeV]. |
| `nHitsDet9` | `Int_t` | Number of EM hits (e⁻, e⁺, γ) on detector 9 for this event. |

#### PMT Segment Indexing

The calorimeter face is 30 cm tall, divided into a 2 (left/right) × 4 (vertical rows) grid. Left/right is defined by the sign of `hitLx`: left (≥ 0) and right (< 0).

```
             hitLx >= 0        hitLx < 0
            (Left side)       (Right side)
          +--------------+--------------+
  +15 cm  |     L1 (0)   |     R1 (4)   |   hitLy >= +7.5 cm
          +--------------+--------------+
  +7.5 cm |     L2 (1)   |     R2 (5)   |   0 <= hitLy < +7.5 cm
          +--------------+--------------+
    0 cm  |     L3 (2)   |     R3 (6)   |   -7.5 cm <= hitLy < 0
          +--------------+--------------+
  -7.5 cm |     L4 (3)   |     R4 (7)   |   hitLy < -7.5 cm
          +--------------+--------------+
  -15 cm
```

Array index: `pmtE[col * 4 + row]` where `col` = 0 (left) or 1 (right), `row` = 0 (top) to 3 (bottom).

### TShower Tree (per hit, optional)

Created only when `writeShowerDiag = true`. Contains one entry per EM hit on detector 9, limited to the first 10,000 events with hits (controlled by `MAX_DIAG_EVENTS`). All diagnostic arrays use `Float_t` to reduce file size. Arrays are dimensioned to `MAX_DEPTH_STEPS = 80`; only indices `0` through `nDepthSteps - 1` contain valid data. Entries beyond that index are zeroed.

| Branch | Type | Description |
|--------|------|-------------|
| `eventIdx` | `Int_t` | Index of the parent event in the input TChain. Cross-references to `TDetResp` entry number. |
| `hitIdx` | `Int_t` | Sequential index of this hit within detector 9 for its event (after EM filtering). Distinguishes multiple hits in the same event. |
| `hitX` | `Float_t` | Hit entry point, horizontal [cm]. |
| `hitY` | `Float_t` | Hit entry point, vertical [cm]. |
| `hitE` | `Float_t` | Incident hit energy [GeV]. |
| `nDepthSteps` | `Int_t` | Number of valid depth slices for this hit. Varies with energy. |
| `depthT[80]` | `Float_t[80]` | Depth of each slice midpoint [X₀]. |
| `sliceE[80]` | `Float_t[80]` | Energy deposited in each depth slice [GeV]. |
| `RC[80]` | `Float_t[80]` | Core radius at each depth [R_M units]. |
| `RT[80]` | `Float_t[80]` | Tail radius at each depth [R_M units]. |
| `p[80]` | `Float_t[80]` | Core weight fraction at each depth [dimensionless, 0–1]. |
| `pmtFrac[640]` | `Float_t[640]` | Fractional energy per PMT segment at each depth. Indexed as `pmtFrac[step * 8 + segment]`. |

## Algorithm Overview

The macro performs two passes over the input data:

**Pass 1 — Beam energy determination**: Scans all events using `hitFluxTrackCoinc()` to find the maximum `evP[0] + evP[1]` for track-qualified coincidences.

**Pass 2 — Shower processing**: For each event:

1. Zero the `pmtE[8]` accumulator.
2. Loop over hits on detector 9. Accept only EM particles (e⁻, e⁺, γ); skip all others.
3. For each EM hit, convert position from meters to cm.
4. If `hitE < E_c` (no shower development), deposit all energy in the segment containing the hit.
5. Otherwise, compute longitudinal profile parameters (T, α, β) from the hit energy.
6. Step through the shower in 1 X₀ depth slices, up to the calorimeter depth (30 cm = 53.5 X₀). Energy beyond this depth is longitudinal leakage and is not deposited:
   - Compute the fractional energy in the slice from the gamma distribution.
   - Compute τ = t/T and evaluate the depth-dependent radial parameters R_C(τ), R_T(τ), p(τ).
   - Compute row fractions using analytic strip integration.
   - Split L/R using core-to-hit-column assignment and tail correction with 32-point GL quadrature.
   - Add the slice energy × segment fraction to the PMT accumulator.
   - Stop when cumulative energy fraction exceeds 99.9% or calorimeter depth is reached.
7. Sum segments into `pmtELeft`, `pmtERight`, `pmtETotal`.
8. Run energy conservation checks (pmtETotal vs sumHitE, and vs beam energy).
9. Apply per-PMT Gaussian energy smearing (σ/E = 12.9%/√E + 1.2%). Clamp negatives to zero.
10. Recompute `pmtELeft`, `pmtERight`, `pmtETotal` from smeared values and fill the output tree.

## Things Encountered

Issues identified and resolved during development, documented here for reference.

### Hadronic particles in the EM shower model

**Symptom**: Events with `pmtETotal` significantly exceeding the beam energy (e.g., 14 GeV from a 10.7 GeV beam), despite the integration passing energy conservation checks.

**Cause**: The flux plane records all particles crossing it, including hadrons (protons, pions, etc.) produced by upstream interactions. Two problems arise from processing these through the EM shower model: (1) the Grindhammer-Peters parameterization describes electromagnetic showers only — hadronic showers have different lateral and longitudinal profiles; (2) `hitE` includes rest mass energy, so a proton with 0.5 GeV kinetic energy has `hitE` = 1.44 GeV (0.938 GeV rest mass + kinetic). When multiple massive secondaries reach the detector, the summed `hitE` can exceed the beam energy from rest masses alone.

**Resolution**: Only EM particles (e⁻, e⁺, γ) are accepted. All other particle types are excluded and counted as `nHitsSkippedNonEM` in the processing summary.

### 2D Gauss-Legendre quadrature failure on sharp peaks

**Symptom**: Energy conservation violations where `pmtETotal` exceeded `sumHitE` by factors of 2–16×. The inflated energy always appeared in the PMT segment containing the hit. Occurred after correcting the Molière radius from 23 cm to 1.6 cm.

**Cause**: The original implementation integrated the 2D radial density `R²/(π(x²+y²+R²)²)` over each PMT segment rectangle using 8×8 Gauss-Legendre quadrature. With the corrected R_M = 1.6 cm, the core component (R_C ~ 0.13 cm) produces a very sharp peak in a rectangle that spans ~7.5 × 28 cm. The GL nodes cannot resolve a feature this narrow relative to the integration domain. When a node lands near the peak, the large Jacobian (proportional to rectangle area) amplifies the sampled value, producing integrals exceeding 1.0 — physically impossible for a normalized density function, but a well-known failure mode of fixed-order quadrature on under-resolved integrands. The bug was invisible at R_M = 23 cm because the integrand was smooth on that scale.

**Resolution**: Replaced the 2D quadrature with a hybrid analytic/numerical method. Row fractions use an exact strip formula. The L/R split uses a decomposition into a symmetric analytic half plus a 1D correction integrated with 32-point GL quadrature. The core component (too narrow for even 1D GL) is assigned entirely to the hit's column, which is >99% accurate since 99% of core energy is within ~4.6 cm and hits are typically 5–7 cm from the L/R boundary.

### X₀ units in critical energy formula (Eq. 5)

**Symptom**: Computed E_c ≈ 0.52 MeV and R_M ≈ 23 cm, versus PDG values of E_c ≈ 7.4 MeV and R_M ≈ 1.6 cm.

**Cause**: Eq. 5 of Grindhammer-Peters uses X₀ in g/cm² (radiation mass thickness), not X₀ in cm (linear radiation length). The NIST value X₀ = 0.5612 cm was passed directly, but the formula requires X₀ × ρ = 6.37 g/cm².

**Resolution**: Precomputed E_c and R_M using the correct units. Added both `X0` (cm) and `X0g` (g/cm²) to `PbConst` namespace with comments explaining the distinction.

### ROOT Cling JIT static initializer failure

**Symptom**: All symbols from the `.h` file unresolved when running the macro interpreted. Error manifested as a massive `Failed to materialize symbols` dump.

**Cause**: The `PbConst` namespace contained `const Double_t` variables initialized with `TMath::Power()` calls. Cling's JIT compiler cannot resolve function calls in namespace-scope static initializers — the guard variables fail to materialize.

**Resolution**: Replaced computed initializers with precomputed literal constants. All `PbConst` values are now hardcoded with the computation shown in comments.

### Cling `std::_Identity` symbol failure

**Symptom**: `symbol '_ZNKSt9_IdentityIiEclERKi' unresolved while linking` when loading `DetectorResponsePlots.C`, even if `PlotEventShower()` was not called.

**Cause**: `std::set<Int_t>` in `PlotEventShower()` triggers Cling to instantiate `std::_Identity<int>`, a GCC internal template that Cling cannot resolve. The error occurs at load time because Cling compiles all functions in the file.

**Resolution**: Replaced `std::set<Int_t>` with linear scan on `std::vector<Int_t>`. Moved the `HitSlice` struct from function-local scope to file scope (Cling also has issues with `std::vector` of locally-defined structs).

### Grindhammer-Peters parameterization breakdown at low energy

**Symptom**: Shower visualization showed enormous rings (18–32 cm radius) at early depth slices for a 23 MeV hit. Debug output revealed R_T growing exponentially to hundreds of R_M, with p collapsing to zero.

**Cause**: The depth-dependent radial parameter R_T (Eq. 25) contains an `exp(k4 × (τ - k2))` term that grows exponentially for τ >> 1. For low-energy particles (E ~ E_c), the shower maximum T is very small (< 1 X₀), so τ = t/T exceeds the parameterization's fitted range by the second depth step. The core weight p drops to zero, leaving the unphysically large tail component to dominate. At E = 23 MeV, by depth step 3: R_T = 122 R_M (197 cm), p = 0.0001. The parameterization was fitted to GeV-scale shower data and is not valid in this regime.

**Resolution**: Added a minimum shower energy threshold of 50 MeV. Particles below this energy are deposited as point-like in the segment containing the hit, bypassing the shower model entirely. This is physically reasonable: a 50 MeV electron in lead has a range of only a few X₀ and produces a shower too compact to spread meaningfully across 7.5 cm PMT segments. The energy carried by sub-threshold particles is typically < 1% of the event total.