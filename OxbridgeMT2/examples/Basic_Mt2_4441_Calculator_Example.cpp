// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/Mt2Calculators.h"
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
  
  Mt2::LorentzVector p_Vis_A     = exampleEvent.p_Vis_A    ();
  Mt2::LorentzVector p_Vis_B     = exampleEvent.p_Vis_B    ();
  Mt2::LorentzVector p_Vis_Other = exampleEvent.p_Vis_Other();
  double           invis_mass  = exampleEvent.invis_mass();
  double           rootS       = exampleEvent.rootS();

  std::cout << "Going to calculate MT2 with\n"
	    << "   p_Vis_A = "     << p_Vis_A     << "\n"
	    << "   p_Vis_B = "     << p_Vis_B     << "\n"
	    << "   p_Vis_Other = " << p_Vis_Other << "\n"
	    << "   rootS = "       << rootS       << "\n"
	    << "   invis_mass = "  << invis_mass  << std::endl;

  // Now that we have some visiable momenta and some missing transverse
  // momentum we can calculate MT2.
  
  // First we create the object that is going to do the calculation
  // of MT2 for us.
  //
  // For this example we will use the most basic 4441 aglorithm
  // defined in the MT2 library:
  Mt2::Basic_Mt2_4441_Calculator mt2Calculator;
  
  // Could tell the MT2 calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  
  // mt2Calculator.setDebug(true);
  
  
  // Now we can actually calculate MT2:
  const double mt2
    = mt2Calculator.mt2_4441(p_Vis_A, p_Vis_B, p_Vis_Other, rootS, invis_mass);
  
  // Now we print out the result:
  std::cout << "ANSWER: mt2 = " << mt2 
	    << " for " << mt2Calculator.algorithmName() << " algorithm"
	    << std::endl; 

  return 0;
}
