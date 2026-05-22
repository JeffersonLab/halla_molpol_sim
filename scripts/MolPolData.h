#ifndef MOLPOL_DATA_H
#define MOLPOL_DATA_H

#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"

#include <cstdio>
#include <fstream>
#include <string>

///////////////////////////////////////////////////////////////
// MolPolData.h
// MolPol simulation ROOT I/O: branch buffers, TChain setup, and
// branch address wiring. Analysis helpers live in MolPolAnalysis.h.
//
// Usage:
//   #include "MolPolData.h"
//   ...
//   TChain *T = SetupMolpolChain("path/to/file.root");  // or file list .txt
//   SetupMolpolBranches(T);
//   // or: TTree *T = (TTree*)f->Get("T"); SetupMolpolBranches(T);
//
// Note: Variable-length arrays (ev.* indexed by ev.npart,
// hit.* indexed by hit.n) are read into static arrays.
// Adjust MAXNPART and MAXNHIT if your simulation exceeds
// these limits.
///////////////////////////////////////////////////////////////

static const Int_t MAXNPART = 10;
static const Int_t MAXNHIT  = 1000;

// --- Event-level scalar branches ---
Int_t    evNpart;
Double_t evPhcom;
Double_t evThcom;
Double_t evXs;
Double_t evAsym;
Double_t evUnpolWght;
Double_t evPolPlusWghtX;
Double_t evPolPlusWghtY;
Double_t evPolPlusWghtZ;
Double_t evPolMinusWghtX;
Double_t evPolMinusWghtY;
Double_t evPolMinusWghtZ;
Double_t evTargMom;

// --- Event-level array branches (indexed by ev.npart) ---
Int_t    evPid[MAXNPART];
Double_t evVx[MAXNPART];
Double_t evVy[MAXNPART];
Double_t evVz[MAXNPART];
Double_t evP[MAXNPART];
Double_t evPx[MAXNPART];
Double_t evPy[MAXNPART];
Double_t evPz[MAXNPART];
Double_t evTh[MAXNPART];
Double_t evPh[MAXNPART];

// --- Hit-level branches ---
Int_t    hitN;
Int_t    hitDet[MAXNHIT];
Int_t    hitVid[MAXNHIT];
Int_t    hitPid[MAXNHIT];
Int_t    hitTrid[MAXNHIT];
Int_t    hitMtrid[MAXNHIT];
Int_t    hitGen[MAXNHIT];
Double_t hitX[MAXNHIT];
Double_t hitY[MAXNHIT];
Double_t hitZ[MAXNHIT];
Double_t hitLx[MAXNHIT];
Double_t hitLy[MAXNHIT];
Double_t hitLz[MAXNHIT];
Double_t hitPx[MAXNHIT];
Double_t hitPy[MAXNHIT];
Double_t hitPz[MAXNHIT];
Double_t hitVx[MAXNHIT];
Double_t hitVy[MAXNHIT];
Double_t hitVz[MAXNHIT];
Double_t hitVdx[MAXNHIT];
Double_t hitVdy[MAXNHIT];
Double_t hitVdz[MAXNHIT];
Double_t hitP[MAXNHIT];
Double_t hitE[MAXNHIT];
Double_t hitM[MAXNHIT];

// --- TBranch pointers: event scalars ---
TBranch *b_evNpart;
TBranch *b_evPhcom;
TBranch *b_evThcom;
TBranch *b_evXs;
TBranch *b_evAsym;
TBranch *b_evUnpolWght;
TBranch *b_evPolPlusWghtX;
TBranch *b_evPolPlusWghtY;
TBranch *b_evPolPlusWghtZ;
TBranch *b_evPolMinusWghtX;
TBranch *b_evPolMinusWghtY;
TBranch *b_evPolMinusWghtZ;
TBranch *b_evTargMom;

// --- TBranch pointers: event arrays ---
TBranch *b_evPid;
TBranch *b_evVx;
TBranch *b_evVy;
TBranch *b_evVz;
TBranch *b_evP;
TBranch *b_evPx;
TBranch *b_evPy;
TBranch *b_evPz;
TBranch *b_evTh;
TBranch *b_evPh;

// --- TBranch pointers: hit branches ---
TBranch *b_hitN;
TBranch *b_hitDet;
TBranch *b_hitVid;
TBranch *b_hitPid;
TBranch *b_hitTrid;
TBranch *b_hitMtrid;
TBranch *b_hitGen;
TBranch *b_hitX;
TBranch *b_hitY;
TBranch *b_hitZ;
TBranch *b_hitLx;
TBranch *b_hitLy;
TBranch *b_hitLz;
TBranch *b_hitPx;
TBranch *b_hitPy;
TBranch *b_hitPz;
TBranch *b_hitVx;
TBranch *b_hitVy;
TBranch *b_hitVz;
TBranch *b_hitVdx;
TBranch *b_hitVdy;
TBranch *b_hitVdz;
TBranch *b_hitP;
TBranch *b_hitE;
TBranch *b_hitM;

///////////////////////////////////////////////////////////////
// SetupMolpolChain()
// Creates TChain("T") and adds files from fileList:
//   - If fileList contains ".root", it is treated as a single ROOT file path.
//   - Otherwise it is treated as a text file with one ROOT path per line.
// Returns NULL if the file list cannot be opened or the chain has no entries.
// Caller owns the returned TChain* (delete when finished).
///////////////////////////////////////////////////////////////

inline TChain *SetupMolpolChain(const char *fileList) {
    TChain *chain = new TChain("T");
    std::string input(fileList);

    if (input.find(".root") != std::string::npos) {
        chain->Add(fileList);
        printf("Added single file: %s\n", fileList);
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
        printf("Added %d files from %s\n", nFiles, fileList);
    }

    if (chain->GetEntries() == 0) {
        printf("Error: no entries found in chain\n");
        delete chain;
        return 0;
    }
    return chain;
}

///////////////////////////////////////////////////////////////
// SetupMolpolBranches()
// Sets all branch addresses for the given TTree pointer.
// Call once after opening the ROOT file and retrieving tree T.
///////////////////////////////////////////////////////////////

inline void SetupMolpolBranches(TTree *tree) {

    // --- Event-level scalar branches ---
    tree->SetBranchAddress("evNpart",         &evNpart,         &b_evNpart);
    tree->SetBranchAddress("evPhcom",         &evPhcom,         &b_evPhcom);
    tree->SetBranchAddress("evThcom",         &evThcom,         &b_evThcom);
    tree->SetBranchAddress("evXs",            &evXs,            &b_evXs);
    tree->SetBranchAddress("evAsym",          &evAsym,          &b_evAsym);
    tree->SetBranchAddress("evUnpolWght",     &evUnpolWght,     &b_evUnpolWght);
    tree->SetBranchAddress("evPolPlusWghtX",  &evPolPlusWghtX,  &b_evPolPlusWghtX);
    tree->SetBranchAddress("evPolPlusWghtY",  &evPolPlusWghtY,  &b_evPolPlusWghtY);
    tree->SetBranchAddress("evPolPlusWghtZ",  &evPolPlusWghtZ,  &b_evPolPlusWghtZ);
    tree->SetBranchAddress("evPolMinusWghtX", &evPolMinusWghtX, &b_evPolMinusWghtX);
    tree->SetBranchAddress("evPolMinusWghtY", &evPolMinusWghtY, &b_evPolMinusWghtY);
    tree->SetBranchAddress("evPolMinusWghtZ", &evPolMinusWghtZ, &b_evPolMinusWghtZ);
    tree->SetBranchAddress("evTargMom",       &evTargMom,       &b_evTargMom);

    // --- Event-level array branches (indexed by ev.npart) ---
    tree->SetBranchAddress("evPid",  evPid,  &b_evPid);
    tree->SetBranchAddress("evVx",   evVx,   &b_evVx);
    tree->SetBranchAddress("evVy",   evVy,   &b_evVy);
    tree->SetBranchAddress("evVz",   evVz,   &b_evVz);
    tree->SetBranchAddress("evP",    evP,    &b_evP);
    tree->SetBranchAddress("evPx",   evPx,   &b_evPx);
    tree->SetBranchAddress("evPy",   evPy,   &b_evPy);
    tree->SetBranchAddress("evPz",   evPz,   &b_evPz);
    tree->SetBranchAddress("evTh",   evTh,   &b_evTh);
    tree->SetBranchAddress("evPh",   evPh,   &b_evPh);

    // --- Hit-level branches ---
    tree->SetBranchAddress("hitN",     &hitN,    &b_hitN);
    tree->SetBranchAddress("hitDet",   hitDet,   &b_hitDet);
    tree->SetBranchAddress("hitVid",   hitVid,   &b_hitVid);
    tree->SetBranchAddress("hitPid",   hitPid,   &b_hitPid);
    tree->SetBranchAddress("hitTrid",  hitTrid,  &b_hitTrid);
    tree->SetBranchAddress("hitMtrid", hitMtrid, &b_hitMtrid);
    tree->SetBranchAddress("hitGen",   hitGen,   &b_hitGen);
    tree->SetBranchAddress("hitX",     hitX,     &b_hitX);
    tree->SetBranchAddress("hitY",     hitY,     &b_hitY);
    tree->SetBranchAddress("hitZ",     hitZ,     &b_hitZ);
    tree->SetBranchAddress("hitLx",    hitLx,    &b_hitLx);
    tree->SetBranchAddress("hitLy",    hitLy,    &b_hitLy);
    tree->SetBranchAddress("hitLz",    hitLz,    &b_hitLz);
    tree->SetBranchAddress("hitPx",    hitPx,    &b_hitPx);
    tree->SetBranchAddress("hitPy",    hitPy,    &b_hitPy);
    tree->SetBranchAddress("hitPz",    hitPz,    &b_hitPz);
    tree->SetBranchAddress("hitVx",    hitVx,    &b_hitVx);
    tree->SetBranchAddress("hitVy",    hitVy,    &b_hitVy);
    tree->SetBranchAddress("hitVz",    hitVz,    &b_hitVz);
    tree->SetBranchAddress("hitVdx",   hitVdx,   &b_hitVdx);
    tree->SetBranchAddress("hitVdy",   hitVdy,   &b_hitVdy);
    tree->SetBranchAddress("hitVdz",   hitVdz,   &b_hitVdz);
    tree->SetBranchAddress("hitP",     hitP,     &b_hitP);
    tree->SetBranchAddress("hitE",     hitE,     &b_hitE);
    tree->SetBranchAddress("hitM",     hitM,     &b_hitM);
}

#endif // MOLPOL_DATA_H