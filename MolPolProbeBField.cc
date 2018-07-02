#include <iostream>
#include <fstream>
#include <vector>
#include "include/MolPolEMField.hh"
#include "G4SystemOfUnits.hh"
#include <string>
#include <cmath>

using namespace std;

int main(){

  std::vector<G4String> fileNames;
  std::vector<G4double> fileScales;
  std::vector<G4double> fileOffsets;

  fileNames.push_back("../TOSCA/q1_6.47kG.table");
  //fileNames.push_back("../TOSCA/q2_6.145kG.table");
  //fileNames.push_back("../TOSCA/q1_6.47kG.table");
  //fileNames.push_back("../TOSCA/q1_6.47kG.table");
  //fileNames.push_back("../TOSCA/lilly_119kG.table");
  fileScales.push_back(1.);//-6.5/5.1
  //fileScales.push_back(-6.0/5.1);
  //fileScales.push_back(2.50/5.1);
  //fileScales.push_back(1.);
  //fileScales.push_back(1.);
  fileOffsets.push_back(75.19*cm);
  //fileOffsets.push_back(140.46*cm);
  //fileOffsets.push_back(209.59*cm);
  //fileOffsets.push_back(274.59*cm);
  //fileOffsets.push_back(422.80*cm);

  MolPolEMField * globalField;

  globalField = new MolPolEMField(fileNames,fileScales,fileOffsets);

  G4bool keepGoing = true;

  G4double px;
  G4double py;
  G4double pz;

  G4cout << " >>> QUICK AND DIRTY FIELD PROBER. CTRL-C TO EXIT. <<< " << G4endl << G4endl;

  while(keepGoing){
    cout << "Enter global x-coordinate: ";
    cin  >> px;
    cout << "Enter global y-coordinate: ";
    cin  >> py;
    cout << "Enter global z-coordinate: ";
    cin  >> pz;

    px *= cm;
    py *= cm;
    pz *= cm;

    //Must use 4-position
    G4double pt[4] = {px,py,pz,0};
    //MolPol EM Field is B & E ... B: 0,1,2 -and- E: 3,4,5
    G4double B[6] = {0.,0.,0.,0.,0.,0.};

    globalField->GetFieldValue(pt,B);

    cout << "Returned field value @ Point ( "
         << px / cm << " , " << py / cm << " , " << pz / cm << " , 0 "
         << " ) cm is ( " << B[0] / tesla << " , " << B[1] / tesla << " , "
         << B[2] / tesla << " ) Teslas" << endl;

    G4cout << G4endl;
  }

}
