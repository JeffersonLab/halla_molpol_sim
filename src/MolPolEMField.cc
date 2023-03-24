#include "MolPolEMField.hh"
#include "MolPolTOSCAField.hh"
#include <iomanip>

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// EMField Constructor :: No Global Field.
MolPolEMField::MolPolEMField()
{
  EField3V.set(0,0,0);
  BField3V.set(0,0,0);
  bUseBFieldMaps = false;
}

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// EMField Constructor :: Using TOSCA maps for global fields.
MolPolEMField::MolPolEMField( std::vector<G4String> files , std::vector<G4double> scales , std::vector<G4double> beamlineOffsets )
{
  EField3V.set(0,0,0);
  BField3V.set(0,0,0);
	bUseBFieldMaps = true;
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

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // B-Field
  Bfield[0]=BField3V.x();
  Bfield[1]=BField3V.y();
  Bfield[2]=BField3V.z();
  if(this->bUseBFieldMaps){
    for (G4int i = 0; i < abs(fFields.size()); i++ ){
      fFields[i]->GetFieldValue(Point, thisB);
      for (G4int j = 0; j < 3; j++) Bsum[j] += thisB[j];
    }
  }
  for (int i = 0; i < 3; i++) Bfield[i] = Bsum[i];

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // E-Field
  G4double  *Efield=&Bfield[3];
  Efield[0]=EField3V.x();
  Efield[1]=EField3V.y();
  Efield[2]=EField3V.z();

  //G4cout << std::setprecision(6) << std::scientific << "MolPolEMField Sending B-field: ( " << Bfield[0] /tesla << " , " << Bfield[1] / tesla << " , " << Bfield[2] / tesla << " ) Tesla"
  //       << std::setprecision(6) << std::fixed << " at Point (" << Point[0] / cm << "," << Point[1] /cm << "," << Point[2] / cm << ") cm" << G4endl << "-------------------------------------------------------------------------------------------------------" << G4endl;

}

//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// Clear TOSCA field vector
void MolPolEMField::clearToscaFields(){
  fFields.clear();
  bUseBFieldMaps = false;
}


//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
// SetTOSCAFields
void MolPolEMField::SetTOSCAFields( std::vector<G4String> files , std::vector<G4double> scales , std::vector<G4double> beamlineOffsets ){
  clearToscaFields();
  bUseBFieldMaps = true;
  for(G4int i = 0; i < abs(files.size()); i++){
    AddNewField(files[i],scales[i],beamlineOffsets[i]);
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
        G4cout << __FUNCTION__ << ": TOSCA field " << name << " was added." << G4endl << G4endl;
    }
}
