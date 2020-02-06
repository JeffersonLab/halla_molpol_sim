# MolPol ROOT Output Tree Structure

The output trees are organized into two main parts. Information about the generated events and information collected at sensitive detectors/flux planes.

## Generated Events Information

**_evNpart_** will be a scaler that defines the size of each of the subsequent arrays with the exception of **_evAsym_** which will is the analyzing power for the moller scattering event. _Note: In the case of a **beam** simulation **evAsym** would be recorded as zero._

To extract the information about the Z position of the vertex for the who particles we would access evVz[0] or evVz[1].

| Branch Name | Description |
|:--------|:---------|
| evNpart |        number of generated particles	   |
| evPid		| Geant4 particle type	   |
| evVx	   | creation vertex, lab frame [m]	   |
| evVy		   |creation vertex, lab frame [m]	   |
| evVz	   |	creation vertex, lab frame [m]	   |
| evP	   |		Particle initial momentum [GeV]	   |
| evPx	   |	Particle initial momentum components, lab frame [GeV]	   |
| evPy	   |	Particle initial momentum components, lab frame [GeV]	   |
| evPz	   |	Particle initial momentum components, lab frame [GeV]	   |
| evTh	   |		Particle initial polar angle [rad]	   |
| evPh	   |		Particle initial azimuthal angle [rad]	   |
| evThcom	   |	Particle initial polar angle in CM frame	   |
| evPhcom	   |	Particle initial azimuthal angle in CM frame	   |
| evAsym   	   |      Azz calucated 	   |

## Hit Data

**_hitN_** is a scalar that defines (to you) the size of each of the arrays that contain information about a particle as it passed through a flux plane or hit a detector. Let's suppose two particles passed through 14 planes and each hit the detector; in this case, there would be 15 entries for each particle so **_hitN_** would be 30. Each recorded hit will tell you which detector **_hitDet_** and the information about the particle as it passed through that detector or flux plane.

| Branch Name | Description |
|:--------|:---------|
| hitN	| 	Number of hits for the event, number of electrons emitted | 
| hitDet	| 	Detector number| 
| hitVid	| 	Volume ID number (not yet implemented)| 
| hitPid	| 	Geant4 particle type| 
| hitTrid | 	Geant4 track ID number (1 = first particle created)| 
| hitMtrid	| Geant4 mother track ID number (0 = particle from gun)| 
| hitX | 	Hit X coordinate, lab frame [m]| 
| hitY | 	Hit Y coordinate, lab frame [m]| 
| hitZ | 	Hit Z coordinate, lab frame [m]| 
| hitLx |  	Hit Local X coordinate, lab frame [m]| 
| hitLy | Hit Local X coordinate, lab frame [m]| 
| hitLz  |   Hit Local X coordinate, lab frame [m]| 
| hitP	| Momentum magnitude of particle [GeV]| 
| hitPx	| Momentum X-component of particle, lab frame [GeV]| 
| hitPy	| Momentum Y-component of particle, lab frame [GeV]| 
| hitPz	| Momentum Z-component of particle, lab frame [GeV]| 
| hitVx | 	Creation X-coordinate of vertex of particles| 
| hitVx | 	Creation Y-coordinate of vertex of particles| 
| hitVx | 	Creation Z-coordinate of vertex of particles| 
| hitE	| 	Energy of particle [GeV]| 
| hitM	| 	Mass of particle [GeV]| 
