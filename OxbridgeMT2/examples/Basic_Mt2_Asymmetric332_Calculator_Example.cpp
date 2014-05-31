// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/Basic_Mt2_Asymmetric332_Calculator.h"
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
  const double                invis_mass = exampleEvent.invis_mass();
  const double theta=atan(1.);
 
  std::cout << "Going to calculate MT2 with\n"
	    << "   ltv_Vis_A  = " << ltv_Vis_A  << "\n"
	    << "   ltv_Vis_B  = " << ltv_Vis_B  << "\n"
	    << "   pT_Miss    = " << pT_Miss    << "\n"
	    << "   invis_mass = " << invis_mass << "\n" 
	    << "   theta      = " << theta << " which is " << theta-3.14159/4 << " from pi/4" << std::endl;

  // Now that we have some visiable momenta and some missing transverse
  // momentum we can calculate MT2.
  
  // First we create the object that is going to do the calculation
  // of MT2 for us.
  //
  // For this example we will use a modification of the "332" aglorithm
  // that was originally defined in SUSYPhys.  Our modification (called
  // Basic_Mt2_Asymmetric332_Calculator) is basically the same as
  // SUSYPhys_Mt2_222_Calculator except that we remove the assumption
  // that visible particles are massless.
  Mt2::Basic_Mt2_Asymmetric332_Calculator mt2Calculator;
  
  // Could tell the MT2 calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  
  // mt2Calculator.setDebug(true);
  
  // Now we can actually calculate MT2:

  const double mt2 
    = mt2Calculator.mt2_Asymmetric332(theta, ltv_Vis_A, ltv_Vis_B, pT_Miss, invis_mass);
  
  // Now we print out the result:
  std::cout << "ANSWER: mt2 = " << mt2 
	    << " for " << mt2Calculator.algorithmName() << " algorithm with theta= "<< theta
	    << std::endl; 

  return 0;
}
