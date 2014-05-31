// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/Basic_MPairProd_Calculator.h"
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
  
  Mt2::LorentzVector p_Vis_A = exampleEvent. p_Vis_A();
  Mt2::LorentzVector p_Vis_B = exampleEvent. p_Vis_B();
  
  std::cout << "Going to calculate MPairProd with a number of copies of\n"
	    << "   p_Vis_A  = " << p_Vis_A  << "\n"
	    << "   p_Vis_B  = " << p_Vis_B  << std::endl;

  // Now that we have some visiable momenta and some missing transverse
  // momentum we can calculate MT2.
  
  // First we create the object that is going to do the calculation
  // of MT2 for us.
  //
  // For this example we will use a modification of the "330" aglorithm
  // that was originally defined in SUSYPhys.  Our modification (called
  // Analytic_Mt2_330_Calculator) is basically the same as
  // SUSYPhys_Mt2_222_Calculator except that we remove the assumption
  // that visible particles are massless.
  Mt2::Basic_MPairProd_Calculator mPairProdCalculator;
  
  // Could tell the MT2 calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  
  mPairProdCalculator.setDebug(true);
  std::vector<Mt2::LorentzVector> interestingParticles;
  for (int i=0; i<13; ++i) {
    interestingParticles.push_back(p_Vis_A);
    interestingParticles.push_back(p_Vis_B);
  }
  // A 16 jet event takes about half a second on pcgd 2007/07/12

  // Now we can actually calculate MT2:
  const double mPairProd  
    = mPairProdCalculator.mPairProd( interestingParticles);
  
  // Now we print out the result:
  std::cout << "ANSWER: mPairProd = " << mPairProd
	    << " for " << mPairProdCalculator.algorithmName() << " algorithm"
	    << std::endl; 

  return 0;
}
