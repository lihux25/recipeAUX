// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/Advanced_Mt2_332_Calculator.h"
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
  double                      invis_mass = exampleEvent.invis_mass();
  
  std::cout << "Going to calculate MT2 with\n"
	    << "   ltv_Vis_A  = " << ltv_Vis_A  << "\n"
	    << "   ltv_Vis_B  = " << ltv_Vis_B  << "\n"
	    << "   pT_Miss    = " << pT_Miss    << "\n"
	    << "   invis_mass = " << invis_mass << std::endl;

  // Now that we have some visiable momenta and some missing transverse
  // momentum we can calculate MT2.
  
  // First we create the object that is going to do the calculation
  // of MT2 for us.
  //
  // For this example we will use the so-called Advanced_Mt2_332_Calculator,
  // which by default is the same as the Basic_Mt2_332_Calculator, except that
  // it can be configured to use a fancy assistor to provide analytic answers
  // to non-balanced solutions.  Don't trust these yet, though if you are
  // using massless visible particles.
  Mt2::Advanced_Mt2_332_Calculator mt2CalculatorDefault; // no fancy assistor
  Mt2::Advanced_Mt2_332_Calculator mt2CalculatorNoAssistor(false); // no fancy assistor
  Mt2::Advanced_Mt2_332_Calculator mt2CalculatorAssistor(true); // uses fancy assistor
  
  // Could tell the MT2 calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  
  // mt2Calculator.setDebug(true);
  
  // Now we can actually calculate MT2:
  const double mt2Default
    = mt2CalculatorDefault   .mt2_332( ltv_Vis_A, ltv_Vis_B, pT_Miss, invis_mass);
  const double mt2NoAssistor
    = mt2CalculatorNoAssistor.mt2_332( ltv_Vis_A, ltv_Vis_B, pT_Miss, invis_mass);
  const double mt2Assistor
    = mt2CalculatorAssistor  .mt2_332( ltv_Vis_A, ltv_Vis_B, pT_Miss, invis_mass);
  
  // Now we print out the result:
  std::cout << "ANSWER:"
	    << "\nmt2(default)\t\t= " << mt2Default 
	    << "\nmt2(noAssistor)\t\t= " << mt2NoAssistor 
	    << "\nmt2(fancyAssistor)\t= " << mt2Assistor 
	    << "\nfor " << mt2CalculatorDefault.algorithmName() << " algorithm. [NB: All answers should be basically the same if everything is working, however the first should be identical to all DPs.  The first two do not use the fancy assistor, while the second does."
	    << std::endl; 

  return 0;
}
