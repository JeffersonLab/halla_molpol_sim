
// ********************************************************************
//
// $Id: MolPolEMField.hh,v 1.0, 2010/12/26   MolPol Exp $
// GEANT4 tag $Name: geant4-09-04 $
//

#ifndef MolPolEMField_H
#define MolPolEMField_H

#include "G4ThreeVector.hh"
#include "G4ElectroMagneticField.hh"
#include "G4SystemOfUnits.hh"
#include "MolPolTOSCAField.hh"

class MolPolEMField : public G4ElectroMagneticField
{
public:

	MolPolEMField() ;
	MolPolEMField( std::vector<G4String> files , std::vector<G4double> scales , std::vector<G4double> beamlineOffsets);
	~MolPolEMField() ;

	inline void GetFieldValue(const G4double Point[4], G4double *Bfield ) const;
	//  Point[4] x,y,z,time
	//  Return as Bfield[0], [1], [2] the magnetic field x, y & z components
	//   and   as Bfield[3], [4], [5] the electric field x, y & z components

	G4bool DoesFieldChangeEnergy() const { return true; }
	//  For field with an electric component this should be true
	//  For pure magnetic field this should be false
	//  Alternative: default safe implementation { return true; }

	inline void SetErDC(G4double val) { ErDC = val; }
	inline G4double GetErDC() const { return ErDC; }

	inline void SetErInner(G4double val) { ErInner = val; }
	inline G4double GetErInner() const { return ErInner; }

	inline void SetEField3V(G4ThreeVector v) { EField3V = v; bUseUniformEField=true;}
	inline G4ThreeVector GetEField3V() const { return EField3V; }

	inline void SetBField3V(G4ThreeVector v) { BField3V = v; bUseUniformBField=true;}
	inline G4ThreeVector GetBField3V() const { return BField3V; }

	inline void clearToscaFields(){
		fFields.clear();
	}

private:

	std::vector<MolPolTOSCAField*> fFields;

	void AddNewField(G4String& name,G4double scale,G4double zOffset);

	bool bUseUniformEField;
	bool bUseUniformBField;
	bool bUseBFieldMaps;

	G4double ErDC;
	G4double ErInner;
	G4ThreeVector EField3V;
	G4ThreeVector BField3V;

};

#endif
