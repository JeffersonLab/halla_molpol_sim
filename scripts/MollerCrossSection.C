#ifndef MOLLER_CROSS_SECTION_C
#define MOLLER_CROSS_SECTION_C

///////////////////////////////////////////////////////////////
// mollerCrossSection.C
// Moller cross section and event-rate helpers (ROOT macro).
//
// Standalone (from project root, verbose by default):
//   root -l -q scripts/mollerCrossSection.C
//   root -l -q 'scripts/mollerCrossSection.C(kMollerCx, 75, 105, 40, 11)'
//   root -l -q 'scripts/mollerCrossSection.C(kMollerRateFe, 75, 105, 40, 11, 1.0)'
//
// Quiet (return value only, last argument kTRUE):
//   RunMollerCrossSection(kMollerRateFe, 75, 105, 40, 11, 1.0, kTRUE);
//
// Calc selector: 0=kMollerCx, 1=kMollerCxBD, 2=kMollerRateFe,
//   3=kMollerRateH, 4=kMollerDiffBD (only theta1 and T used).
///////////////////////////////////////////////////////////////

#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include "TString.h"
#include "Rtypes.h"

enum MollerCxCalc {
    kMollerCx = 0,
    kMollerCxBD,
    kMollerRateFe,
    kMollerRateH,
    kMollerDiffBD
};

const char *MollerCxCalcName(Int_t calc) {
    switch (calc) {
        case kMollerCx:     return "molcx (integrated sigma)";
        case kMollerCxBD:   return "molcxBD (Bjorken-Drell sigma)";
        case kMollerRateFe: return "molrateFe (Fe foil rate)";
        case kMollerRateH:  return "molrateH (1 m LH2 rate)";
        case kMollerDiffBD: return "diffmolcxBD (d sigma/d Omega)";
        default:            return "unknown";
    }
}

const double me = 0.0005109989461;//mass of electron in GeV
const double alpha = 0.0072973525664;//fine structure constant
const double GeV2mb = 0.3894;

double Ecom(double T){
  double E = T + me;
  double Ecm = 2*sqrt(me*(E + me)/2.0);
  return Ecm;
}

double molcx(double theta1=75, double theta2=105, double dphi=40, double T=11,
           Bool_t quiet=kFALSE){
  double low = cos(theta2*TMath::Pi()/180.0);
  double high = cos(theta1*TMath::Pi()/180.0);
  double Ecm = Ecom(T);
  double coeff = 2*TMath::Pi()*GeV2mb*pow(alpha/Ecm,2)*dphi/180.0;
  coeff *= 0.5;
  if (!quiet) {
    cout<<"Center of mass frame energy: "<<Ecm<<" GeV"<<endl;
    cout<<"Coefficient: "<<coeff<<endl;
  }
  TF1 *f = new TF1("f","pow((x*x+3)/(x*x-1),2)", low, high);
  double cx = coeff;
  cx *= f->Integral(low,high);
  if (!quiet)
    printf("Moller scattering cross section for these parameters: %e mb\n", cx);
  return cx;
}

double molcxBD(double theta1=75, double theta2=105, double dphi=40, double T=11,
               Bool_t quiet=kFALSE){
  double low = theta1*TMath::Pi()/180.0;
  double high = theta2*TMath::Pi()/180.0;
  double Ecm = Ecom(T);
  Ecm *= 0.5;
  double coeff = pow(alpha/Ecm,2)/8.0*GeV2mb*2*TMath::Pi()*dphi/180.0;
  coeff *= 0.5;
  if (!quiet) {
    cout<<"Center of mass frame energy: "<<Ecm<<" GeV"<<endl;
    cout<<"Coefficient: "<<coeff<<endl;
  }
  TString fn = Form("sin(x)*((1+pow(cos(x/2.0),4))/pow(sin(x/2.0),4) + 2/pow(sin(x/2)*cos(x/2.0),2)+(1+pow(sin(x/2.0),4))/pow(cos(x/2.0),4))");
  if (!quiet) cout<<fn.Data()<<endl;
  TF1 *f = new TF1("f",fn.Data(), low, high);
  double cx = coeff;
  cx *= f->Integral(low, high);
  if (!quiet)
    printf("Moller scattering cross section for these parameters: %e mb\n", cx);
  return cx;
}

double molrateFe(double theta1=75, double theta2=105, double dphi=40, double T=11,
                 double th=1, Bool_t quiet=kFALSE){
  const double rhoFe = 7874;
  const double avogadro = 6.0221409e26;
  const double atmFe = 55.845;
  const double qe = 1.60217662e-13;
  const double Z_Fe = 26;
  double cx = molcx(theta1, theta2, dphi, T, quiet);
  double mrate = 1e-31 * cx * th * 1e-6 / qe * rhoFe / atmFe * avogadro * Z_Fe;
  if (!quiet)
    printf("Moller scattering rate per microamp for %0.2f micron Fe foil: %i events per second\n",
           th, int(mrate));
  return mrate;
}

double molrateH(double theta1=75, double theta2=105, double dphi=40, double T=11,
                Bool_t quiet=kFALSE){
  const double rhoH = 70.8;
  const double avogadro = 6.0221409e26;
  const double atmH = 1.008;
  const double qe = 1.60217662e-13;
  const double Z_H = 1;
  double cx = molcx(theta1, theta2, dphi, T, quiet);
  double mrate = 1e-31 * cx / qe * rhoH / atmH * avogadro * Z_H;
  if (!quiet)
    printf("Moller scattering rate per microamp for 1 m LH2 taret: %i events per second\n",
           int(mrate));
  return mrate;
}

double diffmolcxBD(double theta=70.9, double T=0.0157, Bool_t quiet=kFALSE){
  double angle = theta * TMath::Pi() / 180.0;
  double Ecm = Ecom(T);
  Ecm *= 0.5;
  if (!quiet) cout<<"Center of mass energy: "<<Ecm<<endl;
  double coeff = pow(alpha/Ecm,2)/8.0*GeV2mb;
  TString fn = Form("((1+pow(cos(x/2.0),4))/pow(sin(x/2.0),4) + 2.0/pow(sin(x/2.0)*cos(x/2.0),2)+(1+pow(sin(x/2.0),4))/pow(cos(x/2.0),4))");
  if (!quiet) cout<<fn.Data()<<endl;
  TF1 *f = new TF1("f",fn.Data(), 0,1);
  double diffcx = coeff;
  diffcx *= f->Eval(angle)*1e-3;
  if (!quiet)
    printf("Differential Moller scattering cross section for these parameters: %e barns\n",
           diffcx);
  return diffcx;
}

double RunMollerCrossSection(Int_t calc,
                             double theta1 = 75, double theta2 = 105,
                             double dphi = 40, double T = 11, double th = 1.0,
                             Bool_t quiet = kFALSE) {
    switch (calc) {
        case kMollerCx:     return molcx(theta1, theta2, dphi, T, quiet);
        case kMollerCxBD:   return molcxBD(theta1, theta2, dphi, T, quiet);
        case kMollerRateFe: return molrateFe(theta1, theta2, dphi, T, th, quiet);
        case kMollerRateH:  return molrateH(theta1, theta2, dphi, T, quiet);
        case kMollerDiffBD: return diffmolcxBD(theta1, T, quiet);
        default:
            printf("RunMollerCrossSection: unknown calc %d\n", calc);
            return 0.0;
    }
}

void MollerCrossSection(Int_t calc = kMollerCx,
                        double theta1 = 75, double theta2 = 105,
                        double dphi = 40, double T = 11, double th = 1.0,
                        Bool_t quiet = kFALSE) {
    if (!quiet) {
        printf("=== %s ===\n", MollerCxCalcName(calc));
        if (calc == kMollerDiffBD)
            printf("  (using theta1=%g deg, T=%g GeV; theta2, dphi, th ignored)\n",
                   theta1, T);
        else if (calc == kMollerRateFe)
            printf("  theta1=%g theta2=%g dphi=%g T=%g GeV th=%g um\n",
                   theta1, theta2, dphi, T, th);
        else if (calc != kMollerDiffBD)
            printf("  theta1=%g theta2=%g dphi=%g T=%g GeV\n",
                   theta1, theta2, dphi, T);
    }

    double result = RunMollerCrossSection(calc, theta1, theta2, dphi, T, th, quiet);

    if (!quiet)
        printf("=== returned value: %e ===\n", result);
}

#endif // MOLLER_CROSS_SECTION_C
