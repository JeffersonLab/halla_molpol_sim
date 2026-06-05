# Møller Cross Section

`MollerCrossSection.C` is a small ROOT macro providing Møller scattering cross
section and event-rate calculations. It is used standalone for quick estimates
and is called by `MolPolAnalysis.C` to report the Møller rate.

## Invocation

```
root -l -q MollerCrossSection.C
root -l -q 'MollerCrossSection.C(kMollerCx, 75, 105, 40, 11)'
root -l -q 'MollerCrossSection.C(kMollerRateFe, 75, 105, 40, 11, 1.0)'
```

Quiet mode (returns the value, no printout) by passing `kTRUE` as the last
argument:

```
RunMollerCrossSection(kMollerRateFe, 75, 105, 40, 11, 1.0, kTRUE);
```

## Calculation Selector

The first argument selects which calculation to run:

| Enum | Value | Calculation |
|------|-------|-------------|
| `kMollerCx` | 0 | Integrated cross section [mb] |
| `kMollerCxBD` | 1 | Bjorken–Drell form of the integrated cross section [mb] |
| `kMollerRateFe` | 2 | Scattering rate for an Fe foil [events/s/µA] |
| `kMollerRateH` | 3 | Scattering rate for a 1 m LH₂ target [events/s/µA] |
| `kMollerDiffBD` | 4 | Differential cross section dσ/dΩ [barns] (uses only `theta1` and `T`) |

## Parameters

| Argument | Default | Meaning |
|----------|---------|---------|
| `theta1` | 75 | Lower CM scattering angle [deg] |
| `theta2` | 105 | Upper CM scattering angle [deg] |
| `dphi` | 40 | Azimuthal acceptance [deg] |
| `T` | 11 | Beam kinetic energy [GeV] |
| `th` | 1.0 | Target thickness [µm] (Fe rate only) |
| `quiet` | `kFALSE` | Suppress printout, return value only |

For `kMollerDiffBD`, only `theta1` (the angle) and `T` are used; the other
geometric arguments are ignored.

## Physics

The cross section uses the standard Møller formula in the center-of-mass frame.
`Ecom(T)` computes the CM energy from the beam kinetic energy. Constants:
electron mass `me`, fine structure constant `alpha`, and the GeV²→mb conversion
`GeV2mb = 0.3894`. The rate calculations fold in target density, Avogadro's
number, atomic mass, and atomic number for iron (`kMollerRateFe`) or hydrogen
(`kMollerRateH`).

## Dependencies

- ROOT (TF1, TMath, TString)