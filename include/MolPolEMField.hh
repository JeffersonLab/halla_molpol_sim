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
#include "G4MagneticField.hh"

#include "MolPolSolenoid.hh"
#include "MolPolDipole.hh"
#include "MolPolTOSCAField.hh"
#include "MolPolQuad.hh"

class MolPolEMField : public G4ElectroMagneticField
{
public:

	MolPolEMField() ;
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

private:

	bool bUseUniformEField;
	bool bUseUniformBField;

	G4double ErDC;
	G4double ErInner;
	G4ThreeVector EField3V;
	G4ThreeVector BField3V;

	bool                bBfieldInUse[6] = {false,false,false,false,false,false};
  MolPolSolenoid*     cSolenoidField;
	G4MagneticField*    cDipoleField;
	G4MagneticField*    cQuad1Field;
	G4MagneticField*    cQuad2Field;
	G4MagneticField*    cQuad3Field;
	G4MagneticField*    cQuad4Field;

public:
	// Solenoid
	void setSolenoidStatus(G4bool state)       {bBfieldInUse[0] = state;G4cout << "  Solenoid status set to: " << bBfieldInUse[0] << G4endl;}
	void setSolenoidObject(MolPolSolenoid* sol){cSolenoidField = sol;}
	void setSolenoidStrength(G4double strength){cSolenoidField->UpdateSolenoid( strength );}

  // Dipole
	void setDipoleStatus(G4bool state){bBfieldInUse[5] = state;G4cout << "  Dipole status set to: " << bBfieldInUse[5] << G4endl;}
	void setDipoleObject(G4int type, G4MagneticField* dip){cDipoleField = dip;}//Probably don't need type here, sub-class should be polymorphic
	void setDipoleStrength(G4int type, G4double strength){
		if(type==1) ((MolPolDipole*)cDipoleField)->setDipoleStrength( strength );
		if(type==2) ((MolPolTOSCAField*)cDipoleField)->setFieldScale( strength );
	}

  // Quad 1
	void setQuad1Status(G4bool state){bBfieldInUse[1] = state;G4cout << "  Quad1 status set to: " << bBfieldInUse[1] << G4endl;}
	void setQuad1Object(G4int type, G4MagneticField* q1){cQuad1Field = q1;}//Probably don't need type here, sub-class should be polymorphic
	void setQuad1Strength(G4int type, G4double strength){
		if(type==1) ((MolPolQuad*)cQuad1Field)->updateQuad( strength );
		if(type==2) ((MolPolTOSCAField*)cQuad1Field)->setFieldScale( strength );
	}

	// Quad 2
	void setQuad2Status(G4bool state){bBfieldInUse[2] = state;G4cout << "  Quad2 status set to: " << bBfieldInUse[2] << G4endl;}
	void setQuad2Object(G4int type, G4MagneticField* q2){cQuad2Field = q2;}//Probably don't need type here, sub-class should be polymorphic
	void setQuad2Strength(G4int type, G4double strength){
		if(type==1) ((MolPolQuad*)cQuad2Field)->updateQuad( strength );
		if(type==2) ((MolPolTOSCAField*)cQuad2Field)->setFieldScale( strength );
	}

	//Quad 3
	void setQuad3Status(G4bool state){bBfieldInUse[3] = state;G4cout << "  Quad3 status set to: " << bBfieldInUse[3] << G4endl;}
	void setQuad3Object(G4int type, G4MagneticField* q3){cQuad3Field = q3;}//Probably don't need type here, sub-class should be polymorphic
	void setQuad3Strength(G4int type, G4double strength){
	  if(type==1) ((MolPolQuad*)cQuad3Field)->updateQuad( strength );
		if(type==2) ((MolPolTOSCAField*)cQuad3Field)->setFieldScale( strength );
	}

	// Quad 4
	void setQuad4Status(G4bool state){bBfieldInUse[4] = state;G4cout << "  Quad4 status set to: " << bBfieldInUse[4] << G4endl;}
	void setQuad4Object(G4int type, G4MagneticField* q4){cQuad4Field = q4;}//Probably don't need type here, sub-class should be polymorphic
	void setQuad4Strength(G4int type, G4double strength){
	  if(type==1) ((MolPolQuad*)cQuad4Field)->updateQuad( strength );
		if(type==2) ((MolPolTOSCAField*)cQuad4Field)->setFieldScale( strength );
	}



};

#endif
