#include <cstdio> 

void GetG3Kin(string fin)
{

  // Read comgeant nt and extract kinematic variables
  // Text output format
  // : BEAME, THETACOM, PHICOM, vtx[3], KE1, KE2, p1[3], p2[3], ANPOW

  // First remove the existing file
  int dummy = remove( "g3input.dat" );

  // output file
  ofstream ofstr("g3input.dat", ios::app);

  TFile* file = new TFile(fin.c_str());
  TTree* tree = (TTree*)file->Get("h2");

  int NMAXTR = 10;       // usuallly there are only 3 tracks (beam particle, two moller electrons)
  int NVERT  = 10;       // usuallly there are only 1 nver
  int ntra;              // number of tracks
  int nver;              // number of vtx
  float bpara[6];        // beam parameters (0-2: PBEAM[3], 3:meh or PBEAM_x/PBEAM_z 4:meh or PBEAM_y/PBEAM_z 5:|PBEAM|)
  float ptra[NMAXTR][3]; // three momentum (0: beam particle)
  float vert[NVERT][3];  // vertex positions 
  float thetcm;          
  float phicm;
  float anpow;

  tree->SetBranchAddress("ntra",   &ntra);
  tree->SetBranchAddress("nver",   &nver);
  tree->SetBranchAddress("bpara",  bpara);
  tree->SetBranchAddress("ptra",   ptra);
  tree->SetBranchAddress("vert",   vert);
  tree->SetBranchAddress("thetcm", &thetcm);
  tree->SetBranchAddress("phicm",  &phicm);
  tree->SetBranchAddress("anpower",&anpow);

  int nentries = tree->GetEntries();
  //  for(int ientry=0; ientry<nentries; ientry++)
  for(int ientry=0; ientry<10000; ientry++)
    {

      tree->GetEntry(ientry);
      if(ntra>3)
        {
          cout << "Number of tracks > 3... skip the event:\t" << ientry << endl;
          continue;
        }

      if(nver>1)
        {
          cout << "Number of vertex > 1... skip the event:\t" << ientry << endl;
          continue;
        }

      ofstr << bpara[5] << " \t"
	    << thetcm   << " \t"
	    << phicm    << " \t"
	    << vert[0][0] << " \t"
	    << vert[0][1] << " \t"
	    << vert[0][2] << " \t"
	    << ptra[1][0] << " \t"
	    << ptra[1][1] << " \t"
	    << ptra[1][2] << " \t"
	    << ptra[2][0] << " \t"
	    << ptra[2][1] << " \t"
	    << ptra[2][2] << " \t"
	    << anpow << endl;     
    }// end of event loop

  ofstr.close();

  file->Close();

  return;

}
