#include "MolPolEMField.hh"
#include "MolPolTOSCAField.hh"

// Left continued usage of bUseUniformBField & EField and added bool for BFieldMaps.
// This is probably not wholly necessary.

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// EMField Constructor :: No Global Field.
MolPolEMField::MolPolEMField()
{
  bUseUniformBField = true;
  bUseUniformEField = true;
	bUseBFieldMaps = false;

  EField3V.set(0,0,0);
  BField3V.set(0,0,0);
}

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// EMField Constructor :: Using TOSCA maps for global fields.
MolPolEMField::MolPolEMField( std::vector<G4String> files , std::vector<G4double> scales , std::vector<G4double> beamlineOffsets )
{
  bUseUniformBField = false;
  bUseUniformEField = true;
	bUseBFieldMaps = true;

  for(int i = 0; i < abs(files.size()); i++){
    G4cout << "file: " << files[i] << " with scale " << scales[i] << " and offset " << beamlineOffsets[i] << G4endl;
  }
  for(G4int i = 0; i < abs(files.size()); i++){
    AddNewField(files[i],scales[i],beamlineOffsets[i]);
  }
}

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// EMField Destructor
MolPolEMField::~MolPolEMField()
{

}


//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Get Field Value
inline void MolPolEMField::GetFieldValue(const G4double Point[4],G4double *Bfield) const
{
  // Member func is 'const' so everything used within must be declared within.
  G4double Bsum[3]  = {0.,0.,0.};
  G4double thisB[3] = {0.,0.,0.};

  //////////////////////////////////////////////////////////
  //get BField
  if(this->bUseUniformBField){
    Bfield[0]=BField3V.x();
    Bfield[1]=BField3V.y();
    Bfield[2]=BField3V.z();
  } else if(this->bUseBFieldMaps){
    std::vector<MolPolTOSCAField*>::const_iterator it = fFields.begin();
    for (it = fFields.begin(); it != fFields.end(); it++){
      (*it)->GetFieldValue(Point, thisB);
      for (int i = 0; i < 3; i++) Bsum[i] += thisB[i];
    }
  }
  for (int i = 0; i < 3; i++) Bfield[i] = Bsum[i];

  //////////////////////////////////////////////////////////
  //get EField,
  double  *Efield=&Bfield[3];
  if(this->bUseUniformEField || this->bUseBFieldMaps)
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

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Constucts new TOSCAField object for each field currently given and adds those
// objects to vector fFields.
void MolPolEMField::AddNewField(G4String& name,G4double scale,G4double zOffset)
{
    MolPolTOSCAField *newToscaField = new MolPolTOSCAField(name,scale,zOffset);
    if (newToscaField->IsInit()) {
        fFields.push_back(newToscaField);
        G4cout << __FUNCTION__ << ": field " << name << " was added." << G4endl;
    }
}
