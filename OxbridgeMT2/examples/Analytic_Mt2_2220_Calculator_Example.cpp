// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/Analytic_Mt2_2220_Calculator.h"
#include "Mt2/Basic_Mt2_332_Calculator.h"
#include "Mt2/ChengHanBisect_Mt2_332_Calculator.h"
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

  exampleEvent.resetSpecial222();
  
 //while(true) 
 {
 //Mt2::TwoVector pT_Vis_A = exampleEvent.  pT_Vis_A();
 // Mt2::TwoVector pT_Vis_B = exampleEvent.  pT_Vis_B();
 // Mt2::TwoVector  pT_Miss = exampleEvent.   pT_Miss();
 // double       invis_mass = exampleEvent.invis_mass();
  //Mt2::TwoVector pT_Vis_A(1.2*cos(0.9),1.2*sin(0.9));
  //Mt2::TwoVector pT_Vis_B(1.7*cos(1.1),1.7*sin(1.1));
  Mt2::TwoVector pT_Vis_A(atan(3.14159/4.8)*cos(4.705),atan(3.14159/4.8)*sin(4.705));
  Mt2::TwoVector pT_Vis_B(1*cos(-4.705),1*sin(-4.705));
  Mt2::LorentzTransverseVector ltva(pT_Vis_A);
  Mt2::LorentzTransverseVector ltvb(pT_Vis_B);
  Mt2::TwoVector  pT_Miss(1,0); 
  double       invis_mass = 0;
  
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
  // defined in Analytic: (note that this algorithm assumed that the
  // visiable particles are massless ... or equivalently assumed they were
  // of neglible mass.)
  Mt2::Analytic_Mt2_2220_Calculator mt2Calculator;
  Mt2::Basic_Mt2_332_Calculator comparisonCalculator;
  Mt2::ChengHanBisect_Mt2_332_Calculator comparisonCalculator2;
  
  // Could tell the MT2 calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  
  // mt2Calculator.setDebug(true);
  
  // Now we can actually calculate MT2:
  double mt2;

  try {
    mt2 = mt2Calculator.mt2_2220( pT_Vis_A, pT_Vis_B, pT_Miss);
  } catch (...) {
    mt2 = -100;
  }
  const double comparisonCalcByOtherMethod
    = comparisonCalculator.mt2_332( ltva, ltvb, pT_Miss,invis_mass);
  const double comparisonCalcByOtherMethod2
    = comparisonCalculator2.mt2_332(ltva, ltvb, pT_Miss,invis_mass);
  
  // Now we print out the result:
  std::cout << "COMPARE: " << mt2 << " " << comparisonCalcByOtherMethod << " " << comparisonCalcByOtherMethod2 << " " << ((fabs(comparisonCalcByOtherMethod)<0.1) ? " SMALL" : "") << std::endl;

  std::cout << "ANSWER: mt2 = " << mt2 
	    << " for " << mt2Calculator.algorithmName() << " algorithm\n"
 	<< " which compares to " << comparisonCalcByOtherMethod << " for " << comparisonCalculator.algorithmName()
	    << std::endl; 

	std::cout << "-----------------------------" << std::endl;

	exampleEvent.resetRandom222();

}
  return 0;
}
