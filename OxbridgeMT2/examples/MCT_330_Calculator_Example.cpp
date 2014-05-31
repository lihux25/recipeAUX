// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/MCT_330_Calculator.h"
#include "Mt2/MCTT_332_Calculator.h"
#include "Mt2/MCTll_332_Calculator.h"
#include <iostream>

int main(int argc, char * argv[]) {

  // For the example. we now need some momenta and masses from which
  // to calculate MCT.
  //
  // In "reality" we would get these momenta from an ntuple or
  // from a physics event.
  //
  // As this is just an example program, we will instead get
  // some "example" momenta from the class "ExampleEvent" 
  // defined in "ExampleEvent.h" as follows:
  
  ExampleEvent exampleEvent;

  exampleEvent.resetZhenA_332();

  for (int i=0; i<20; ++i) {

    if (i==0) {
      Mt2::LorentzTransverseVector ltv_Vis_A = exampleEvent.ltv_Vis_A();
      Mt2::LorentzTransverseVector ltv_Vis_B = exampleEvent.ltv_Vis_B();
      Mt2::TwoVector                     utm = exampleEvent.pT_Vis_Other();
      double                      invis_mass = exampleEvent.invis_mass();
  
      std::cout << "Going to calculate MCT and MCTT with\n"
		<< "   ltv_Vis_A  = " << ltv_Vis_A  << "\n"
		<< "   ltv_Vis_B  = " << ltv_Vis_B  << "\n"
		<< "   utm        = " << utm  << "\n" << std::endl;

      // Now we can actually calculate MCT:
      const double mct 
	= mct_330( ltv_Vis_A, ltv_Vis_B);

      const double mctT 
	= mctt_332 (ltv_Vis_A, ltv_Vis_B, utm);
  
      const double mctll
	= mctll_332(ltv_Vis_A, ltv_Vis_B, utm);
  
      // Now we print out the result:
      std::cout << "ANSWER: mct = " << mct
		<< " for mct algorithm"
		<< std::endl; 
      std::cout << "ANSWER: mctT = " << mctT
		<< " for mctT algorithm"
		<< std::endl; 
      std::cout << "ANSWER: mct|| = " << mctll
		<< " for mct|| algorithm"
		<< std::endl; 

    } else {

      if (i<=10) {
	static bool first=true;
	if (first) {
	  first = false;
	  std::cout << "Events of type 330" << std::endl;
	}
	exampleEvent.resetRandom330();
      } else {
	static bool first=true;
	if (first) {
	  first = false;
	  std::cout << "Events of type 332" << std::endl;
	}
	exampleEvent.resetRandom332();
      }
      Mt2::LorentzTransverseVector ltv_Vis_A = exampleEvent.ltv_Vis_A();
      Mt2::LorentzTransverseVector ltv_Vis_B = exampleEvent.ltv_Vis_B();
      Mt2::TwoVector                     utm = exampleEvent.pT_Vis_Other();
      double                      invis_mass = exampleEvent.invis_mass();
  

      // Now we can actually calculate MCT:
      const double mct 
	= mct_330( ltv_Vis_A, ltv_Vis_B);

      const double mctT 
	= mctt_332 (ltv_Vis_A, ltv_Vis_B, utm);
  
      const double mctll
	= mctll_332(ltv_Vis_A, ltv_Vis_B, utm);
  
      // Now we print out the result:
      std::cout << "LIST: mct = " << mct
		<< " mctT = " << mctT 
		<< " mct|| = " << mctll
		<< std::endl; 


    }
  }
  return 0;
}
