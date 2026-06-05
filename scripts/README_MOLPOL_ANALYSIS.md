# MolPol Analysis

`MolPolAnalysis.C` computes the longitudinal analyzing power (A_zz) and event
acceptance of the MolPol Møller polarimeter directly from the raw Geant4
simulation output, using a **geometric** flux-plane cut on detector 9. It is
the counterpart to `MolPolDetectorResponseAnalysis.C`, which performs the
equivalent analysis on the *modeled* calorimeter response; the two agree when
run with comparable gates.

This analysis works on the flux-plane hits themselves (momentum and track
information), without modeling shower development or PMT response.

## Invocation

```
root -l 'MolPolAnalysis.C("input.root")'
root -l 'MolPolAnalysis.C("input.root", "12")'                  // rows 1,2 active
root -l 'MolPolAnalysis.C("input.root", "1234", 0.25, kFALSE)'  // momentum-cut gates, 0.25 GeV
root -l 'MolPolAnalysis.C("filelist.txt")'                      // text list of .root files
```

| Argument | Type | Default | Meaning |
|----------|------|---------|---------|
| `fileList` | `const char*` | — | A `.root` file or a text file listing `.root` files |
| `activeRowStr` | `const char*` | `"1234"` | Which PMT rows are active (e.g. `"12"` = rows 1 and 2) |
| `momentumCut` | `Double_t` | `0.0` | Momentum threshold [GeV] for the momentum-cut gates |
| `useTrackGates` | `Bool_t` | `kTRUE` | `kTRUE` = track-coincidence gates; `kFALSE` = momentum-cut gates |

## Gate Modes

The flux plane (detector 9) is segmented into four PMT rows by `hitLy`
boundaries, always read out in Left/Right pairs. The active rows are selected
with `activeRowStr`; only hits in active rows contribute. Two gate modes are
available:

- **Track gates** (`useTrackGates = kTRUE`): an arm fires if it has a hit from a
  generated Møller track (`hitTrid` 1 or 2). This is the cleanest selection of
  primary Møller electrons. Helpers: `hitFluxTrackLeft`, `hitFluxTrackRight`,
  `hitFluxTrackCoinc`.

- **Momentum-cut gates** (`useTrackGates = kFALSE`): an arm fires if the summed
  non-neutron momentum on that side exceeds `momentumCut`. Helpers:
  `hitFluxPLeft`, `hitFluxPRight`, `hitFluxPCoinc`. The gate uses a strict `>`
  comparison, so `momentumCut = 0.0` requires nonzero momentum.

Left and Right are defined looking upstream (toward the target).

## Channels

Three channels are accumulated per event:

- **Coincidence**: both arms fire.
- **Left singles**: left arm fires.
- **Right singles**: right arm fires.

For each passing event the longitudinal asymmetry is computed from the polarized
weights and accumulated. The analyzing power is then derived from the
polarization factor and the known target polarization
(`targetPolarization = 0.08014`).

## Output

Printed to the terminal:

- **Header**: file, gate mode, momentum cut (momentum mode only), active rows
  and segments, number of entries.
- **Kinematics**: beam energy (max summed track momentum), θ_CM and φ_CM
  ranges, and the Møller scattering rate.
- **Summary** per channel: raw event fraction, Møller-weighted and
  weighted event fractions, polarization factor, and mean analyzing power with
  uncertainty.

The Møller analyzing power peaks at 7/9 ≈ 0.778 at θ_CM = 90°, which is the
natural sanity check for the reported mean analyzing power.

## Active-Row Helpers (in `MolPolAnalysis.h`)

| Function | Purpose |
|----------|---------|
| `SetActiveRowsFromStr(str)` | Set the active rows from a string such as `"1234"` |
| `SetPmtRowActive(row)` / `SetPmtRowInactive(row)` | Toggle a single row (explicit set/unset for clarity) |
| `IsActivePmtRow(i)` | Whether hit `i` lands in an active row |
| `GetActiveRowStr()` / `PrintActiveRows()` | Report the active rows |

The flux gate helpers (`hitFluxPLeft/Right/Coinc`,
`hitFluxTrackLeft/Right/Coinc`) and the `MolPolGateMode` enum also live in
`MolPolAnalysis.h`, shared with the detector-response analysis.

## Dependencies

- ROOT (TTree, TChain, TH1F, TCanvas)
- `MolPolData.h` — raw simulation tree reader
- `MolPolAnalysis.h` — gate helpers, active-row helpers, accumulators
- `MollerCrossSection.C` — Møller rate calculation (see its README)