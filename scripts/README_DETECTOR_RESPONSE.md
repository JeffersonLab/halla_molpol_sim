# MolPol Detector Response Model

This document describes the detector response simulation chain for the MolPol
Møller polarimeter lead–fiber calorimeter. The chain takes flux-plane hit data
from the MolPol Geant4 simulation (detector 9), models the electromagnetic
shower development and fiber penetration in the calorimeter, and produces a new
ROOT file containing the per-PMT energy response. A set of analysis and
inspection scripts then operate on that output.

The code is organized into two stages:

1. **`MolPolDetectorResponse.C` / `.h`** — reads the raw simulation output and
   writes a `*_detresponse.root` file containing the modeled PMT response.
2. **Analysis scripts** — `MolPolDetectorResponseAnalysis.C` (asymmetry and
   spectrum analysis) and `MolPolDetectorResponseShowerInspection.C`
   (interactive per-event shower inspection), both of which read the
   `*_detresponse.root` file.

---

# Part 1 — `MolPolDetectorResponse.C`

## Purpose

`MolPolDetectorResponse.C` converts raw flux-plane hits on detector 9 into a
modeled calorimeter response. The flux plane in the simulation is a 1 nm thick
surface; it records particles crossing it but does not itself model energy
deposition. This script supplies the missing physics: it develops each
electromagnetic particle into a full shower, distributes the deposited energy
across the 2×4 grid of PMT segments, models where in the block the shower
starts (fiber penetration), and applies calorimeter energy resolution.

## Invocation

```
root -l 'MolPolDetectorResponse.C("input.root")'
root -l 'MolPolDetectorResponse.C("input.root", true)'           // write shower diagnostics
root -l 'MolPolDetectorResponse.C("input.root", false, 12345)'   // custom RNG seed
root -l 'MolPolDetectorResponse.C("filelist.txt")'               // text file list of .root files
```

Arguments:

| Argument | Type | Default | Meaning |
|----------|------|---------|---------|
| `fileList` | `const char*` | — | Single `.root` file, or a text file listing `.root` files (one per line) |
| `writeShowerDiag` | `Bool_t` | `false` | If true, writes the `TShower` per-hit diagnostic tree (capped at 10,000 events) |
| `seed` | `UInt_t` | `1800707365` | Seed for the `TRandom3` generator (fiber MC and energy smearing). Fixed default for reproducibility |

The output file name is the input name with `.root` replaced by
`_detresponse.root`.

## Physics Models

### Electromagnetic Shower Parameterization

Showers are modeled using the Grindhammer–Peters parameterization
(G. Grindhammer and S. Peters, arXiv:hep-ex/0001020, 2000), which describes the
average longitudinal and radial energy profile of an electromagnetic shower in
a homogeneous medium. The shower is treated in two stages so that the expensive
profile computation is performed once per hit and reused across all Monte Carlo
throws.

**Stage 1 — `ComputeShowerProfile()`** computes the full shower shape,
independent of where the shower starts in the block:

- The longitudinal profile follows a gamma distribution, with the shower
  maximum `T`, shape parameter `α`, and rate `β = (α−1)/T` taken from the
  homogeneous-medium formulas in the paper (functions of `ln y`, where
  `y = E/E_c`).
- The profile is integrated in 1 X₀ slices to give the fractional energy
  deposited in each slice (`sliceFrac`).
- At each depth the radial profile is described by a core and a tail term,
  with the core radius `RC`, tail radius `RT`, and core weight `p` evaluated
  as functions of the scaled depth `τ = t/T`.
- For each slice the radial profile is integrated over the finite block
  boundaries to give the fraction of that slice's energy collected by each of
  the 8 PMT segments (`frac[step*8 + seg]`). See **Segment Integration Method**
  below.

The profile is computed out to the natural shower extent (up to `3.5·T` or
20 X₀), *without* a block-depth cap; the cap is applied during accumulation.

**Stage 2 — `AccumulateShowerInBlock()`** shifts a computed profile by a start
depth offset and accumulates the energy into the 8 PMT segments. Slices whose
physical depth (start depth + slice depth) exceeds the calorimeter depth
(`caloDepthX0`) are dropped — this is longitudinal leakage.

Particles below `MIN_SHOWER_E_GEV = 50 MeV` are not showered: they deposit
their energy point-like in the segment containing the hit. This avoids the
regime where the parameterization is unreliable (the physical critical energy
for lead is `E_c = 7.35 MeV`, but the radial parameters become unphysical well
above that, so a conservative 50 MeV floor is used).

### Segment Integration Method

The PMT segment energy fractions are computed by a semi-analytic integration of
the radial profile over each finite segment rectangle:

- The transverse integral over the *y* direction (the segment row boundaries)
  is done analytically.
- The transverse integral over the *x* direction (the finite block half-width
  ±9 cm) is done with 32-point Gauss–Legendre quadrature.

This hybrid approach was adopted because a pure 2D Gauss–Legendre quadrature
failed to capture the sharp core peak of the radial profile (the core radius
is small compared with the segment size).

### Fiber Penetration Model

The calorimeter is a lead–fiber sandwich: roughly 4/5 lead and 1/5
scintillating fiber by volume on the face. Because the radiation length of the
fiber (`fiberX0 = 45 cm`) is far longer than that of lead (`X₀ = 0.56 cm`), a
particle that enters a fiber can travel a significant distance before
showering. This produces a distribution of shower start depths, which in turn
broadens and adds a low-energy tail to the PMT response — the effect this model
is designed to capture.

`ComputeFiberPenetration()` returns the effective shower start depth [cm] for a
single particle, via the following Monte Carlo procedure:

1. **Material selection** (`DetermineStartMaterial`): a flat random draw selects
   fiber with probability `fiberVolFrac = 0.20`, otherwise lead.

2. **Start depth sampling** (`SampleShowerStartDepth`): the depth at which the
   shower initiates is drawn from an exponential distribution with mean equal
   to the radiation length of the chosen material:
   `P(z) = (1/X₀) exp(−z/X₀)`, sampled with ROOT's `rng->Exp(X₀)`.

3. **Fiber geometry** (fiber hits only):
   - **Transverse distance** (`SampleFiberTransverseDistance`): the chord
     length through the fiber, sampled via the inverse CDF
     `d_trans = D · √(1 − u²)`, `u ~ Uniform(0,1)`. This peaks near the full
     diameter `D = fiberDiam` and falls toward zero at the edge.
   - **Incidence angle**: the angle of the particle relative to the fiber axis,
     `θ_incident = |atan2(p_T, p_z) − detAngleRad|`, where the detector face is
     tilted at `detAngleDeg = 7.5°` from the beam axis.
   - **Available depth**: `z_available = d_trans / sin(θ_incident)` — how far
     the particle can travel inside the fiber before exiting its side into
     lead.
   - **Clamping**: if the sampled `z_start` exceeds `z_available`, the particle
     exits the fiber and showers in lead, so the effective start depth is set
     to `z_available`.

4. **Pass-through**: if the effective start depth exceeds the block depth
   (`caloDepth = 30 cm`), the particle traverses the entire block without
   showering and the function returns `−1.0`. With the default geometry this
   requires both `z_start > 30` and `z_available > 30` (a near-parallel fiber
   trajectory), and is rare.

Each helper is its own function so the physics can be refined independently.
A `// TODO` flag marks the chord-length model as a candidate for further
refinement.

### Monte Carlo Throws

For every event the script performs `DR_NTHROWS = 20` independent throws. Each
throw samples a fresh set of fiber-penetration start depths (one per hit) and a
fresh set of per-PMT energy smearings, producing 20 statistically independent
realizations of the PMT response for the same underlying event. The shower
*profiles* are computed once per hit and shared across throws; only the start
depth and smearing change. This 20-throw structure is what produces the
asymmetric low-energy tails that should be compared against FADC pulse-height
data.

### Energy Resolution Smearing

After accumulation, each PMT segment energy is independently smeared with a
Gaussian whose width follows the scintillating-fiber calorimeter resolution
(R. Livan, M. Vercesi, R. Wigmans, 1995):

```
σ/E = 12.9% / √E ⊕ 1.2%
```

Smearing is applied per segment, per throw. Negative smeared values are clamped
to zero. The Left/Right/Total sums are recomputed from the smeared segment
values.

### Energy Conservation and Leakage

On throw 0 (before smearing), the script checks that the total deposited energy
does not exceed the sum of incident hit energies (an "energy violation" would
indicate a bug) and tracks how much energy is lost to longitudinal leakage. The
**energy leakage fraction** reported in the summary is energy-weighted:

```
leakage fraction = 1 − (total deposited energy) / (total incident energy)
```

summed over all events. It is a fraction of energy, not a fraction of events.

## Material and Geometry Constants

Defined in the `PbConst` and `DetGeom` namespaces in
`MolPolDetectorResponse.h`. All are precomputed literals (rather than evaluated
via `TMath` functions) to avoid ROOT Cling static-initializer failures.

| Constant | Value | Meaning |
|----------|-------|---------|
| `PbConst::Z` | 82 | Lead atomic number |
| `PbConst::A` | 207.2 | Lead mass number |
| `PbConst::rho` | 11.35 g/cm³ | Lead density |
| `PbConst::X0` | 0.5612 cm | Lead radiation length |
| `PbConst::Ec` | 7.3548 MeV | Lead critical energy |
| `PbConst::RM` | 1.6176 cm | Lead Molière radius |
| `DetGeom::caloDepth` | 30.0 cm | Physical block depth (≈53.5 X₀) |
| `DetGeom::blockHalfX` | 9.0 cm | Block half-width |
| `DetGeom::blockHalfY` | 15.0 cm | Block half-height |
| `DetGeom::detAngleDeg` | 7.5° | Detector face angle from beam axis |
| `DetGeom::fiberDiam` | 0.1 cm | Fiber diameter (1 mm) |
| `DetGeom::fiberX0` | 45.0 cm | Fiber radiation length |
| `DetGeom::fiberVolFrac` | 0.20 | Fiber volume fraction |
| `MAX_DEPTH_STEPS` | 80 | Max longitudinal slices |
| `DR_MAXNHIT` | 20 | Max EM hits per event |
| `DR_NTHROWS` | 20 | MC throws per event |
| `MIN_SHOWER_E_GEV` | 0.050 | Showering threshold [GeV] |

## Diagnostic Histograms

After processing, the script draws a six-panel diagnostic canvas (`cFiberDiag`,
palette 55) and writes the histograms to the output file:

1. **z_start** — shower start depth, total (black) / lead (red) / fiber (blue),
   log scale.
2. **z_available** — available fiber depth, log scale.
3. **z_{start,raw}^{MC} vs z_avail** — 2D distribution (COLZ) with a magenta
   dashed L-corner marking the pass-through boundary at (30, 30).
4. **d_trans** — fiber transverse (chord) distance.
5. **θ_incident by particle type** — All (black) / Gen Møllers, `hitTrid` 1 or 2
   (red) / Other e⁻ (blue) / γ (green) / proton (magenta), log scale.
6. **Shower vs Pass-through scatter** — showered hits (black) vs pass-through
   hits (red) in the z_avail / z_start_raw plane.

Log-scale panels have their minimum fixed to 0.9 so single-count bins are
visible. Histograms are detached from the output file (`SetDirectory(0)`) before
the file is closed so they survive for drawing.

---

# Part 2 — Analysis Scripts

Both analysis scripts read the `*_detresponse.root` file produced by Part 1 and
use the shared reader functions in `MolPolDetectorResponse.h`
(`SetupDetRespChain`, `SetupDetRespBranches`, and the `dr_*` global buffers).

## 2a — `MolPolDetectorResponseAnalysis.C`

### Purpose

Computes the longitudinal analyzing power (A_zz) and PMT energy spectra from the
modeled detector response, with selectable active PMT rows and an energy
threshold. It is the detector-response analogue of the geometric `MolPolAnalysis`
(see the separate `README_MOLPOL_ANALYSIS.md`), and the two can be compared
apples-to-apples when run with comparable gates.

### Invocation

```
root -l 'MolPolDetectorResponseAnalysis.C("input_detresponse.root")'
root -l 'MolPolDetectorResponseAnalysis.C("input_detresponse.root", 0.5)'          // 0.5 GeV cut
root -l 'MolPolDetectorResponseAnalysis.C("input_detresponse.root", 0.5, "12")'    // rows 1,2 only
```

| Argument | Type | Default | Meaning |
|----------|------|---------|---------|
| `fileList` | `const char*` | — | A `*_detresponse.root` file or text list |
| `energyCut` | `Float_t` | `0.0` | Per-arm energy threshold [GeV]; gate uses strict `>` |
| `activeRowStr` | `const char*` | `"1234"` | Which PMT rows are active (e.g. `"12"` = rows 1,2) |

The analysis currently reads **throw 0** (`dr_pmtE[0..7]`) of the 20-throw
output. Selecting other throws, or averaging across all 20, is a possible
extension for FADC tail comparisons.

### Channels and Gates

Three asymmetry channels are computed in a single event loop, gated on the
active-row energy sums against `energyCut` (strict `>`):

- **Coincidence**: both arms above the cut.
- **Left singles**: left arm above the cut.
- **Right singles**: right arm above the cut.

For each event the longitudinal asymmetry is computed from the polarized
weights (`evPolPlusWghtZ`, `evPolMinusWghtZ`); only polarized events contribute
to the A_zz histograms (unpolarized events have zero polarized weight and yield
NaN, which is rejected).

### Output Canvases

- **cPmtSeg** (4×2): individual PMT segment energy spectra; inactive rows shaded.
- **cPmtSums** (3×1 linear) and **cPmtSumsLog** (3×1 log, y-min fixed to 0.1):
  Left / Right / Total energy sums, with Active PMT (blue solid) overlaid on
  All PMT (red dashed).
- **cPmtLR**: Left-vs-Right energy heatmap for active rows.
- **cAzz** (3×1): A_zz spectra for the three channels, range [0, 7/9], with the
  x-axis auto-ranged to the populated region. The Møller analyzing power peaks
  at 7/9 ≈ 0.778 (θ_CM = 90°).

The helper functions (`ParseActiveRows`, the `DetRespAsymAccum` accumulator, the
`DetRespAsymChannel` enum, `Solve90PctRadius`, and the `HitSlice` struct) live
in `MolPolDetectorResponseAnalysis.h`.

## 2b — `MolPolDetectorResponseShowerInspection.C`

### Purpose

Interactive, per-event inspection of the shower diagnostics. Requires that the
input file was produced with `writeShowerDiag = true` so that the `TShower`
tree is present.

### Invocation

```
root -l 'MolPolDetectorResponseShowerInspection.C("input_detresponse.root")'
```

On load it opens the file, wires up the `TDetResp` branches via
`SetupDetRespBranches`, builds an index that maps each unique `eventIdx` in
`TShower` to its arm energies (throw 0) in `TDetResp`, and prints the available
interactive commands.

### Interactive Commands

| Command | Description |
|---------|-------------|
| `ListShowerEvents(N, M, energyCut)` | List `N` events starting at position `M` (1-based), labeling each by gate type (Coincidence / Left / Right / None) against `energyCut` |
| `PlotEventShower(eventNum)` | 3D visualization of the shower(s) in one event |
| `PlotShowerParams(eventNum)` | Grindhammer–Peters parameter plots (slice energy, RC, RT, p vs depth) for the event |

`eventNum` is 1-based and indexes into the ordered list of unique events in
`TShower`. Because the diagnostic tree stores one entry per hit (not per throw),
multi-hit events appear as consecutive `TShower` entries sharing the same
`eventIdx` but with different `hitIdx`.

---

# Output ROOT File Structure (`*_detresponse.root`)

The output file contains one main tree (`TDetResp`), an optional diagnostic
tree (`TShower`, only when `writeShowerDiag = true`), and the diagnostic
histograms listed above.

## `TDetResp` tree — one entry per event

Throws are indexed `k = 0 … 19`; PMT segments are indexed `s = 0 … 7`.

| Branch | Type | Size | Description |
|--------|------|------|-------------|
| `evXs` | `Double_t` | scalar | Event cross-section weight (passed through) |
| `evAsym` | `Double_t` | scalar | Event asymmetry (passed through) |
| `evUnpolWght` | `Double_t` | scalar | Unpolarized weight |
| `evPolPlusWghtX/Y/Z` | `Double_t` | scalar | Polarized (+) weights, X/Y/Z |
| `evPolMinusWghtX/Y/Z` | `Double_t` | scalar | Polarized (−) weights, X/Y/Z |
| `pmtE` | `Float_t` | `[160]` | Per-segment deposited energy [GeV], post-smear, indexed `[k*8 + s]` |
| `pmtETotal` | `Float_t` | `[20]` | Total energy per throw [GeV] |
| `pmtELeft` | `Float_t` | `[20]` | Left-arm sum (segments 0–3) per throw [GeV] |
| `pmtERight` | `Float_t` | `[20]` | Right-arm sum (segments 4–7) per throw [GeV] |
| `hitZ` | `Float_t` | `[20]` | Representative shower start depth per throw [cm] |
| `nHitsDet9` | `Int_t` | scalar | Number of EM hits on detector 9 in this event |
| `hitPx` | `Double_t` | `[nHitsDet9]` | Per-hit momentum x [GeV] |
| `hitPy` | `Double_t` | `[nHitsDet9]` | Per-hit momentum y [GeV] |
| `hitPz` | `Double_t` | `[nHitsDet9]` | Per-hit momentum z [GeV] |

### PMT Segment Indexing

The 8 PMT segments form a 2×4 grid (2 columns × 4 rows). Segment index
`s = col*4 + row`, with column 0 = Left (segments 0–3) and column 1 = Right
(segments 4–7). Rows are ordered from top to bottom by the `rowYBound`
boundaries. Within a throw `k`, segment `s` is stored at `pmtE[k*8 + s]`.

## `TShower` tree — one entry per hit (optional)

Present only when `writeShowerDiag = true`, capped at 10,000 showering events.
One entry per EM hit; the shower profile is the same for all throws, so the
profile is stored once and the 20 per-throw start depths are stored in `hitZ`.

| Branch | Type | Size | Description |
|--------|------|------|-------------|
| `eventIdx` | `Int_t` | scalar | Source event index (== `TDetResp` entry number) |
| `hitIdx` | `Int_t` | scalar | Hit index within the event (0-based) |
| `hitX`, `hitY` | `Float_t` | scalar | Hit entry point on the face [cm] |
| `hitE` | `Float_t` | scalar | Hit energy [GeV] |
| `nDepthSteps` | `Int_t` | scalar | Number of valid longitudinal slices |
| `depthT` | `Float_t` | `[80]` | Slice depth midpoint [X₀] |
| `sliceE` | `Float_t` | `[80]` | Energy deposited in each slice [GeV] |
| `RC` | `Float_t` | `[80]` | Core radius per slice [Molière radii] |
| `RT` | `Float_t` | `[80]` | Tail radius per slice [Molière radii] |
| `p` | `Float_t` | `[80]` | Core weight per slice |
| `pmtFrac` | `Float_t` | `[640]` | Per-slice PMT fractions, indexed `[step*8 + seg]` |
| `hitZ` | `Float_t` | `[20]` | Shower start depth for each throw [cm] (−1 = pass-through) |

---

# Dependencies

- ROOT (TTree, TChain, TBranch, TRandom3, TH1F/TH2F, TGraph, TCanvas)
- `MolPolData.h` — raw simulation tree reader (input side)
- `MolPolAnalysis.h` — active-row helpers shared with the geometric analysis
- `MolPolDetectorResponse.h` — shower model, fiber model, and `TDetResp` reader

# References

- G. Grindhammer and S. Peters, *The Parameterized Simulation of
  Electromagnetic Showers in Homogeneous and Sampling Calorimeters*,
  arXiv:hep-ex/0001020 (2000).
- R. Livan, M. Vercesi, R. Wigmans, *Scintillating-Fibre Calorimetry* (1995).