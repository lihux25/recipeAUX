// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/Analytic_Mt2_330_Calculator.h"
#include "Mt2/Basic_Mt2_332_Calculator.h"
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
  

  unsigned int count=0;
  while(++count <= 100000) {
	exampleEvent.resetRandom332();
 
  const Mt2::LorentzTransverseVector & ltv_Vis_A = exampleEvent.ltv_Vis_A   ();
  const Mt2::LorentzTransverseVector & ltv_Vis_B = exampleEvent.ltv_Vis_B   ();
  const Mt2::TwoVector               & pT_ISR    = exampleEvent.pT_Vis_Other();  
  const Mt2::TwoVector               & pT_Miss   = exampleEvent.pT_Miss     ();
  const double                      invis_mass   = exampleEvent.invis_mass  ();
 
 
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
  // For this example we will use a modification of the "330" aglorithm
  // that was originally defined in SUSYPhys.  Our modification (called
  // Analytic_Mt2_330_Calculator) is basically the same as
  // SUSYPhys_Mt2_222_Calculator except that we remove the assumption
  // that visible particles are massless.
  Mt2::Analytic_Mt2_330_Calculator mt2Calculator_Analytic;
  
  // Could tell the MT2 calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  
  // mt2Calculator.setDebug(true);
  
  // Now we can actually calculate MT2:
  const double mt2_Analytic
    = mt2Calculator_Analytic.mt2_330( ltv_Vis_A, ltv_Vis_B, invis_mass);

  const Mt2::Analytic_Mt2_330_Calculator::Notes notes= mt2Calculator_Analytic.notes();
 
  Mt2::Basic_Mt2_332_Calculator mt2Calculator_Basic;
  const double mt2_Basic
    = mt2Calculator_Basic.mt2_332(ltv_Vis_A, ltv_Vis_B, pT_Miss, invis_mass);
  
  // Now we print out the result:
  std::cout << "ANSWER: mt2_Analytic\t" << mt2_Analytic
	    << " mt2_Basic\t" << mt2_Basic
	<< " diff\t" << mt2_Basic - mt2_Analytic
	<< " pT_ISR\t" << pT_ISR.pt()
	<< "\tsolnType " << notes.soln
	<< "\tcalcdRootS " << notes.rootS
	<< "\tgamma " << notes.gamma
 << " with\t"
            << " ltv_Vis_A\t" << ltv_Vis_A
            << " ltv_Vis_B\t" << ltv_Vis_B 
            << " pT_Miss\t" << pT_Miss    
            << " invis_mass\t" << invis_mass << std::endl
	    << std::endl;

  }
  return 0;
}
