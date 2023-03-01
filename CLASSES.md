# `MolPol Class Descriptions`
Brief description of the MolPol classes.

## `MolPolDetector()`

Custom class that inherits from **G4VSensitiveDetector()**. Takes a G4String for the name of the detector and a integer assignment. Detectors are registered to the SensitiveDetector manager in **MolPolDetectorConstruction()::Construct()**.

## `MolPolDetectorConstruction()`
Polarimeter detector construction. Major methods of the class described below:

### DefineGeometryCommands()
Generic messenger to pass geometry specific macro commands for simulations. See [macro summary](https://github.com/JeffersonLab/halla_molpol_sim#molpol-macro-usage) in the main README file under "Geometry Modifiables".

### ConstructMaterials()
Inserted as a method to compartmentalize the building of materials for the spectrometer. Materials are constructed here and then neatly called at the beginning of **Construct()**. Allows the whole thing to be much more readable.

### Construct()
Basic Geant4 construction method from inherited class. Returns the world volume.

## `MolPolDetectorHit()`
Inherits from **G4Hit()** and describes the information which will be collected for each detector hit. 

## `MolPolDipole()`
Creates idealized dipole field within assigned volume.  Class utilized within **MolPolEMFieldSetup()** (see below).

*Note: There exists commented out code in the MolPolDipole class for snaked fringe fields. I don't know who commented it out or why. Perhaps this is something that we should look at in the future.*

## `MolPolEMFieldMessenger()`
Messenger class for magnetic fields.  Please see the tabled list  found here under 'Field Information' for a list of options. This is straight forward. 

## `MolPolEMFieldSetup()`
Sets up the fields for the quads, dipole and holding field. 

## `MolPolEMField()`
I'm unsure what the purpose of this class is. At first glance it appears to be unused as all other active fields have their own classes Quads, Dipole and 'Solenoid'. This is where information is passed to from the **MolPolEMFieldMessenger**.

## `MolPolEvent()`
Custom class which holds the information for each event created in **MolPolPrimaryActionGenerator()**.

## `MolPolEventAction()`
Inherits from the standard Geant4 class for Event Action. The important method in this class is the **EndOfEventAction()** where the Hits Collection is sorted through and MolPolIO is instructed to store that information in its members.

## `MolPolIO()`
Creates or re-creates (if the same name is specified) the ROOT file for simulation output.  

#### InitializeTree()
Builds the tree structure of the ROOT file.

#### AddDetectorHit() 
Takes a pointer to a **MolPolDetectorHit** object and stores the data in the **MolPolIO** members.

#### SetEventData()
Takes the event data from **MolPolEvent** object.

#### FillTree()
Standard ROOT command, initiates the writing of a tree branch from objects specified during branch construction in **InitializeTree()**.

#### Flush()
This is somewhat misrepresented. It does not flush the member arrays of **MolPolIO** but rather sets the number of hits sent from the Hits Collection (fNDetHit) recorded during **MolPolEventAction::EndOfEventAction()** back to zero. It is this value, fNDetHit, that is used to specify the (portion) size of the buffered member arrays of **MolPolIO** which will be recorded.

> fTree->Branch("hitN", &fNDetHit, "hit.n/I");
> 
> fTree->Branch("hitDet",> &fDetHit_det, "hit.det[hit.n]/I");
>  
> fTree->Branch("hitVid",> &fDetHit_id, "hit.vid[hit.n]/I");

As can be seen in the snippet the value stored at  **&fNDetHit** is stored in the root branch as **hit.n** which is then used as the array size for subsequent branch entries. This is how root knows how many entries from the **MolPolIO** hit member arrays to save to the tree. 


## `MolPolMessenger()`
Messenger class for non-field non-geometry macros.

## `MolPolPrimaryGeneratorAction()`
This is the event generator.

|Generator| Brief Description |
|----|------|
| Moller | Generates Moller Pairs, general documentation on the computational methods used can be found in the Swartz paper [Observation of Target Electron Momentum Effects in Single-Arm Møller Polarimetry](https://www.sciencedirect.com/science/article/abs/pii/0168900295003843), in the J.P. Alexander paper [Radiative corrections to the  Z0  resonance](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.37.56), with an excellent reference on the Hydrogenlike Wavefunctions used for the iron target specified in the Swartz paper written by Don Jones in the Github Repository found [here](https://github.com/jonesdc76/MollerPolarimetry).    |
| Beam   | Generates single electrons using relevant options such as ***fx***, ***fy***, ***fz***, ***beamRot__*** and ***beamE***. Can be useful if we're ever interested in what Geant4 says should be happening with beam on target .|
| LUND   | Generates individual events based on specifications in LUND file.|
| Inelastics (Coming soon) | Inelastic generator modeled after remoll inelastic generator. 

### Levchuk Effect 
Has been converted from old hydrogenic models to Hartree-Fock calculated models. Details of this are covered in [NIMA](https://www.sciencedirect.com/science/article/abs/pii/S0168900222007987) or, alternatively, [arXiv](https://arxiv.org/abs/2207.02150).

## `MolPolQuad()`
Used to assign a quadrupole field in a specific volume.  Self explanatory. Takes the *field gradient* from pole tip to opposing pole tip, the *origin of the field* (where should the center of the field be), a *rotation matrix* (if necessary) and the radius of the quadrupole defined from center to one of the pole tips. 

*Note: The offsets for the field entered in the simulation macro are passed to* ***MolPolEMFieldSetup***  *and then passed to* ***MolPolQuad*** *in the origin of the field*.

## `MolPolRunAction()`
Basic functionality of a RunAction class in Geant4. Instructs **MolPolIO** to initialize the root tree, times the overall simulation, and instructs **MolPolIO** to write the tree.


## `MolPolSolenoid()`
Creates instance of the holding field. Loads solenoid.map which is normalized to 1 Tesla. The scaling factor for the map is passed by the macro */field/setQ6T*.

The class returns a linearly interpolated value of the field from the field map. 


## `MolPolSteppingAction()`
Basic SteppingAction functionality. Two macros are associated with this. Information can be found in the 'Stepping Action Options' in the table found [here](https://github.com/JeffersonLab/halla_molpol_sim/blob/master/README.md#molpol-macro-usage).
