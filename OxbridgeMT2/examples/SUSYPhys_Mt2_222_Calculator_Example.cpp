// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/SUSYPhys_Mt2_222_Calculator.h"
#include <iostream>

int main(int argc, char * argv[]) {

  // For the example. we now need some momenta and masses from which
  // to calculate MT2.
  //
  // In "reality" we would get these momenta from an ntuple or
  // from a physics event.
  //
  // As this is just an example program, we will instead get
  // some "example" momenta from the class "ExampleEvent" 
  // defined in "ExampleEvent.h"as follows:
  
  ExampleEvent exampleEvent;
  
  Mt2::TwoVector pT_Vis_A = exampleEvent.  pT_Vis_A();
  Mt2::TwoVector pT_Vis_B = exampleEvent.  pT_Vis_B();
  Mt2::TwoVector  pT_Miss = exampleEvent.   pT_Miss();
  double       invis_mass = exampleEvent.invis_mass();
  
  std::cout << "Going to calculate MT2 with\n"
	    << "   pT_Vis_A  = " << pT_Vis_A  << "\n"
	    << "   pT_Vis_B  = " << pT_Vis_B  << "\n"
	    << "   pT_Miss    = " << pT_Miss    << "\n"
	    << "   invis_mass = " << invis_mass << std::endl;

  // Now that we have some visiable momenta and some missing transverse
  // momentum we can calculate MT2.
  
  // First we create the object that is going to do the calculation
  // of MT2 for us.
  //
  // For this example we will use the "332" aglorithm that was originally
  // defined in SUSYPhys: (note that this algorithm assumed that the
  // visiable particles are massless ... or equivalently assumed they were
  // of neglible mass.)
  Mt2::SUSYPhys_Mt2_222_Calculator mt2Calculator;
  
  // Could tell the MT2 calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  
  // mt2Calculator.setDebug(true);
  
  // Now we can actually calculate MT2:
  const double mt2
    = mt2Calculator.mt2_222( pT_Vis_A, pT_Vis_B, pT_Miss, invis_mass);
  
  // Now we print out the result:
  std::cout << "ANSWER: mt2 = " << mt2 
	    << " for " << mt2Calculator.algorithmName() << " algorithm"
	    << std::endl; 

  return 0;
}
