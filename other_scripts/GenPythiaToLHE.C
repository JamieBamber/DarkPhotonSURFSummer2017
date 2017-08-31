#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

//#include "HepMC/GenEvent.h"
//#include "HepMC/IO_GenEvent.h"
//#include "HepMC/Units.h"

#include <TDatime.h>

#include <stdio.h>
#include <string.h>
#include <unistd.h>

using namespace Pythia8; 

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc != 3) {
    cerr << " To run the code provide the name of the input pythia card and the output LHE file. \n"
	 << " example: ./GenPythiaToHepMC data/pythiaCards/EXO/RSGraviton_gg_EXAMPLE.pythia outFile.lhe \n" << endl;
    return 1;
  }

  // Check that the provided input name corresponds to an existing file.
  ifstream is(argv[1]);  
  if (!is) {
    cerr << " Command-line file " << argv[1] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Interface for conversion from Pythia8::Event to HepMC one. 
  //  HepMC::Pythia8ToHepMC ToHepMC;
  // Switch off warnings for parton-level events.
  //  ToHepMC.set_print_inconsistency(false);
  //  ToHepMC.set_free_parton_warnings(false);

  // Specify file where HepMC events will be stored.                                                           
  //  HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);

  // Confirm that external files will be used for input and output.
  cout << " PYTHIA settings will be read from file " << argv[1] << endl;
 
  // Generator. 
  Pythia pythia;

  // Read in commands from external file.
  pythia.readFile(argv[1]);    

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

  // Initialize. Beam parameters set in .pythia file.
  pythia.init();

  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nShow    = int(nEvent/100);
  int nAbort   = 10; 
  /*
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
  */

  // Create an LHAup object that can access relevant information in pythia.                                            
  LHAupFromPYTHIA8 myLHA(&pythia.process, &pythia.info);

  // Open a file on which LHEF events should be stored, and write header.                                              
  myLHA.openLHEF(argv[2]);

  // Store initialization info in the LHAup object.                                                                    
  myLHA.setInit();

  // Write out this initialization info on the file.                                                                   
  myLHA.initLHEF();

  // Begin event loop.
  int nPace = max(1, nEvent / max(1, nShow) ); 
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (nShow > 0 && iEvent%nPace == 0) 
      cout << " Now begin event " << iEvent << endl;

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }

    // Store event info in the LHAup object.                                                                           
    myLHA.setEvent();

    // Write out this event info on the file.                                                                          
    // With optional argument (verbose =) false the file is smaller.                                                   
    myLHA.eventLHEF();

    // End of event loop.
  }
  
  // Give statistics. 
  pythia.stat();

  // Update the cross section info based on Monte Carlo integration during run.                                        
  myLHA.updateSigma();
  
  // Write endtag. Overwrite initialization info with new cross sections.                                              
  myLHA.closeLHEF(true);

  // Done.
  return 0;
}
