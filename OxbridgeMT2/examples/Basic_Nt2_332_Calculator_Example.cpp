// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/Basic_Nt2_332_Calculator.h"
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
  
  Mt2::LorentzTransverseVector ltv_Vis_A = exampleEvent. ltv_Vis_A();
  Mt2::LorentzTransverseVector ltv_Vis_B = exampleEvent. ltv_Vis_B();
  Mt2::TwoVector                 pT_Miss = exampleEvent.   pT_Miss();
  double                    parent1_mass = exampleEvent.invis_mass(); // kludge
  double                    parent2_mass = exampleEvent.invis_mass(); // kludge
  
  std::cout << "Going to calculate MT2 with\n"
	    << "   ltv_Vis_A  = " << ltv_Vis_A  << "\n"
	    << "   ltv_Vis_B  = " << ltv_Vis_B  << "\n"
	    << "   pT_Miss    = " << pT_Miss    << "\n"
	    << "   parent1_mass = " << parent1_mass << "\n" 
	    << "   parent2_mass = " << parent2_mass << std::endl;

  // Now that we have some visiable momenta and some missing transverse
  // momentum we can calculate MT2.
  
  // First we create the object that is going to do the calculation
  // of MT2 for us.
  //
  // For this example we will use a modification of the "332" aglorithm
  // that was originally defined in SUSYPhys.  Our modification (called
  // Basic_Nt2_332_Calculator) is basically the same as
  // SUSYPhys_Mt2_222_Calculator except that we remove the assumption
  // that visible particles are massless.
  Mt2::Basic_Nt2_332_Calculator mt2Calculator;
  
  // Could tell the MT2 calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  
  // mt2Calculator.setDebug(true);
  
  // Now we can actually calculate MT2:
  const double mt2 
    = mt2Calculator.nt2_332( ltv_Vis_A, ltv_Vis_B, pT_Miss, parent1_mass, parent2_mass);
  
  // Now we print out the result:
  std::cout << "ANSWER: mt2 = " << mt2 
	    << " for " << mt2Calculator.algorithmName() << " algorithm"
	    << std::endl; 

  return 0;
}
