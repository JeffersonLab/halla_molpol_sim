// ********************************************************************
//
// $Id: MolPolEMField.hh,v 1.0, 2010/12/26   MolPol Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//   User Field class Setup implementation.
//
//
#include "MolPolEMField.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

MolPolEMField::MolPolEMField()
{

  EField3V.set(0,0,0);
  BField3V.set(0,0,0);
}

//////////////////////////////////////////////////////////////////////////
//
//  Deconstructors:
MolPolEMField::~MolPolEMField()
{

}


////////////////////////////////////////////////////////////////////////////
//input Point[4] (x,y,z,t)
//
inline void MolPolEMField::GetFieldValue(const G4double Point[4],G4double *Bfield) const
{
  //////////////////////////////////////////////////////////
  //get BField
  if(this->bUseUniformBField)
    {
      Bfield[0]=BField3V.x();
      Bfield[1]=BField3V.y();
      Bfield[2]=BField3V.z();
    }
  else
    {
      double pB[3],pPos[3]={Point[0]/cm,Point[1]/cm,Point[2]/cm};  //turn into cm
      for(int i=0;i<3;i++) Bfield[i]=0.0;  //reset

    }

  //////////////////////////////////////////////////////////
  //get EFiled,

  double  *Efield=&Bfield[3];
  if(this->bUseUniformEField)
    {
      Efield[0]=EField3V.x();
      Efield[1]=EField3V.y();
      Efield[2]=EField3V.z();
    }
  else
    {
      for(int i=0;i<3;i++) Efield[i]=0.;
    }
}
