// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/MC_330_Calculator.h"
#include <iostream>

int main(int argc, char * argv[]) {

  // For the example. we now need some momenta and masses from which
  // to calculate MC.
  //
  // In "reality" we would get these momenta from an ntuple or
  // from a physics event.
  //
  // As this is just an example program, we will instead get
  // some "example" momenta from the class "ExampleEvent" 
  // defined in "ExampleEvent.h" as follows:
  
  ExampleEvent exampleEvent;
  
  Mt2::LorentzVector lv_Vis_A = exampleEvent. p_Vis_A();
  Mt2::LorentzVector lv_Vis_B = exampleEvent. p_Vis_B();
  Mt2::TwoVector                 pT_Miss = exampleEvent.   pT_Miss();
  double                      invis_mass = exampleEvent.invis_mass();
  
  std::cout << "Going to calculate MC with\n"
	    << "   lv_Vis_A  = " << lv_Vis_A  << "\n"
	    << "   lv_Vis_B  = " << lv_Vis_B  << "\n" << std::endl;

  // Now we can actually calculate MCT:
  const double mc 
    = mc_330( lv_Vis_A, lv_Vis_B);
  
  // Now we print out the result:
  std::cout << "ANSWER: mc = " << mc
	    << " for mc algorithm"
	    << std::endl; 

  return 0;
}
