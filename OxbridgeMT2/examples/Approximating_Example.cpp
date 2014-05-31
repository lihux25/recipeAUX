// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/Basic_Mt2_332_Calculator.h"
#include "Mt2/Analytic_Mt2_330_Calculator.h"
#include "Mt2/Mt2ApproximatingAdapter_332_from_330.h"

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
  // For this example we will use a modification of the "332" aglorithm
  // that was originally defined in SUSYPhys.  Our modification (called
  // Basic_Mt2_332_Calculator) is basically the same as
  // SUSYPhys_Mt2_222_Calculator except that we remove the assumption
  // that visible particles are massless.
  Mt2::   Basic_Mt2_332_Calculator                            basic_mt2_332_calculator;
  Mt2::Analytic_Mt2_330_Calculator                         analytic_mt2_330_calculator;
  Mt2::Mt2ApproximatingAdapter_332_from_330  analyticPretendingToBe_mt2_332_calculator(analytic_mt2_330_calculator);

  // Could tell the MT2 calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  
  // mt2Calculator.setDebug(true);
  
  // Now we can actually calculate MT2:
  const double mt2_proper_332
    =                  basic_mt2_332_calculator.mt2_332( ltv_Vis_A, ltv_Vis_B, pT_Miss, invis_mass);
  const double mt2_pretend_332
    = analyticPretendingToBe_mt2_332_calculator.mt2_332( ltv_Vis_A, ltv_Vis_B, pT_Miss, invis_mass);
  
  // Now we print out the result:
  std::cout << "ANSWER: mt2_proper_332  = " << mt2_proper_332
	    << " for " << basic_mt2_332_calculator.algorithmName() << " algorithm"
	    << std::endl; 
  std::cout << "ANSWER: mt2_pretend_332 = " << mt2_pretend_332
	    << " for " << analyticPretendingToBe_mt2_332_calculator.algorithmName() << " algorithm"
	    << std::endl; 

  return 0;
}
