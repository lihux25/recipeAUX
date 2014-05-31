// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/Basic_MtGen_332_Calculator.h"
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
  
  Mt2::LorentzTransverseVector ltv_Vis_A = exampleEvent. ltv_Vis_A();
  Mt2::LorentzTransverseVector ltv_Vis_B = exampleEvent. ltv_Vis_B();
  Mt2::LorentzVector lv_Vis_A = exampleEvent. p_Vis_A();
  Mt2::LorentzVector lv_Vis_B = exampleEvent. p_Vis_B();
  Mt2::TwoVector                 pT_Miss = exampleEvent.   pT_Miss();
  double                      invis_mass = exampleEvent.invis_mass();
  
  std::cout << "Going to calculate MT2 with\n"
	    << "   ltv_Vis_A  = " << ltv_Vis_A  << "\n"
	    << "   ltv_Vis_B  = " << ltv_Vis_B  << "\n"
	    << "   lv_Vis_A  = " << lv_Vis_A  << "\n"
	    << "   lv_Vis_B  = " << lv_Vis_B  << "\n"
	    << "   pT_Miss    = " << pT_Miss    << "\n"
	    << "   invis_mass = " << invis_mass << std::endl;

  // Now that we have some visiable momenta and some missing transverse
  // momentum we can calculate MT2.
  
  // First we create the object that is going to do the calculation
  // of MT2 for us.
  //
  // For this example we will use a modification of the "332" aglorithm
  // that was originally defined in SUSYPhys.  Our modification (called
  // Analytic_Mt2_332_Calculator) is basically the same as
  // SUSYPhys_Mt2_222_Calculator except that we remove the assumption
  // that visible particles are massless.

  // Create an Mt2_332 calculator which we will ask the MtGen calculator to use later:
  Mt2::ChengHanBisect_Mt2_332_Calculator mt2_332_Calculator;

  // Ask MtGen to use the above Mt2_332_Calculator:
  Mt2::Basic_MtGen_332_Calculator mtGenCalculator(mt2_332_Calculator);
  
  // Could tell the MTGen calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  mtGenCalculator.setDebug(true);

  Mt2::TwoVector ptMiss;
  std::vector<Mt2::LorentzTransverseVector> interestingTParticles;
  std::vector<Mt2::LorentzVector> interestingParticles;
  for (int i=0; i<3; ++i) {
    interestingParticles.push_back(lv_Vis_A);
    interestingParticles.push_back(lv_Vis_B);
    interestingTParticles.push_back(ltv_Vis_A);
    interestingTParticles.push_back(ltv_Vis_B);
    ptMiss -= ltv_Vis_A.transverse();
    ptMiss -= ltv_Vis_B.transverse();
  }
  // A 16 jet event takes about half a second on pcgd 2007/07/12

  // Now we can actually calculate MT2:
  const double mtGenT 
    = mtGenCalculator.mtGen_332( interestingTParticles, ptMiss, invis_mass);
  const double mtGen
    = mtGenCalculator.mtGen_332( interestingParticles, ptMiss, invis_mass);
  
  // Now we print out the result:
  std::cout << "ANSWER: mtGenT = " << mtGenT << std::endl
	    << "ANSWER: mtGen = " << mtGen
	    << " for " << mtGenCalculator.algorithmName() << " algorithm"
	    << std::endl; 

  return 0;
}
