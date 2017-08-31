// Generate a sample of events, starting from 
// a sample of LHE and write the output in a ROOT file

#include "Pythia8/Pythia.h"
//#include "Pythia8/Pythia8ToHepMC.h"
#include "Pythia8Plugins/HepMC2.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/Units.h"

#include <TFile.h>
#include <TTree.h>
#include <GenTree.hh>
#include <GenCandidateFiller.hh>

#include <stdio.h>
#include <string.h>

using namespace Pythia8; 

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc != 4) {
    cerr << " Unexpected number of command-line arguments. \n You are"
         << " expected to provide the name of an input card, an input LHE file, and an output file. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Check that the provided input names corresponds to existing files
  ifstream is1(argv[1]);  
  if (!is1) {
    cerr << " The requested PYTHIA card " << argv[1] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }
  ifstream is2(argv[2]);  
  if (!is2) {
    cerr << " The requested LHE event file " << argv[2] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Confirm that external files will be used for input and output.
  cout << " PYTHIA settings will be read from file " << argv[1] << endl;

  // Interface for conversion from Pythia8::Event to HepMC one. 
  HepMC::Pythia8ToHepMC ToHepMC;
  // Switch off warnings for parton-level events.
  //  ToHepMC.set_print_inconsistency(false);
  //  ToHepMC.set_free_parton_warnings(false);

  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(argv[3], std::ios::out);

  // Generator. 
  Pythia pythia;

  // Read in commands from external file.
  //pythia.readFile(argv[1]);    
  
  // Read the LHE file
  char name[256];
  pythia.readString("Beams:frameType = 4");
  sprintf(name,"Beams:LHEF = %s", argv[2]);
  cout << name << endl;
  pythia.readString(name);

  // set seed
  int jobpid = getpid();
  TDatime *now = new TDatime();
  int today = now->GetDate();
  int clock = now->GetTime();
  int myseed = today+clock+jobpid*1000;
  if(myseed>900000000) myseed = myseed - 900000000;
  pythia.readString("Random:setSeed=on");
  char command[512];
  sprintf(command,"Random:seed=%i",myseed);
  pythia.readString(command);

  // Initialize Les Houches Event File run. List initialization information.
  pythia.init();

  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nList    = pythia.mode("Main:numberToList");
  int nShow    = pythia.mode("Main:timesToShow");
  int nAbort   = pythia.mode("Main:timesAllowErrors"); 
  bool showCS  = pythia.flag("Main:showChangedSettings");
  bool showAS  = pythia.flag("Main:showAllSettings");
  bool showCPD = pythia.flag("Main:showChangedParticleData");
  bool showAPD = pythia.flag("Main:showAllParticleData");
  

  // List settings.
  if (showCS) pythia.settings.listChanged();
  if (showAS) pythia.settings.listAll();

  // List particle data.  
  if (showCPD) pythia.particleData.listChanged();
  if (showAPD) pythia.particleData.listAll();

  // Begin event loop.
  int nPace = max(1, nEvent / max(1, nShow) ); 
  int iAbort = 0;
  for (int iEvent = 0; ; ++iEvent) {
    if (nShow > 0 && iEvent%nPace == 0) 
      cout << " Now begin event " << iEvent << endl;
    
    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break; 
      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }
 
    // List first few events.
    if (iEvent < nList) { 
      pythia.LHAeventList();               
      pythia.info.list();
      pythia.process.list();          
      pythia.event.list();
    }

    // Construct new empty HepMC event. 
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    // Fill HepMC event, including PDF info.
    ToHepMC.fill_next_event( pythia, hepmcevt );

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;

    // End of event loop.
  }

  // Give statistics. 
  pythia.stat();
 
  // Done.
  return 0;
}
