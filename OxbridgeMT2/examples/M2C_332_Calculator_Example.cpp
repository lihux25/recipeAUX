// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#include "ExampleEvent.h"
#include "Mt2/Basic_M2C_332_Calculator.h"

// Executes a single call to the M2C function.

// M2C part of the Code written by Mario A. Serna Jr
// Last updated on 08 Feb 2008
// The variable M2C was introduced in
//   G.Ross and M.Serna, "Mass Determination of New States at Hadron Colliders", 
//   arXiv:0712.0943 [hep-ph] 06 Dec 2007

int main(int argc, char *argv[]) 
  {
  // For the example. we now need some momenta and masses from which
  // to calculate M2C.

  // As this is just an example program, we will instead get
  // some "example" momenta 
  
  
  ExampleEvent exampleEvent;
  // If you want to use the example event uncomment these
//  Mt2::LorentzTransverseVector ltv_Vis_A=exampleEvent. ltv_Vis_A();
//  Mt2::LorentzTransverseVector ltv_Vis_B=exampleEvent. ltv_Vis_A();
//  Mt2::TwoVector               pT_Miss = exampleEvent. pT_Miss();

   // here is a manually entered event.
    Mt2::LorentzTransverseVector ltv_Vis_A(60.1104,3.71501,15.7146);
    Mt2::LorentzTransverseVector ltv_Vis_B(48.804,-43.3881,12.8803);
    Mt2::TwoVector               pT_Miss(39.6731,-28.5949);
    double DeltaM = 133-69; // GeV  for Demo 

  double m2C_result;
  
   // display the case being calculated
    std::cout << "Going to calculate M2C with\n"
	    << "   ltv_Vis_A  = " << ltv_Vis_A  << "\n"
	    << "   ltv_Vis_B  = " << ltv_Vis_B  << "\n"
	    << "   pT_Miss    = " << pT_Miss    << "\n"
        << "   Delta M    = " << DeltaM     << "\n"
	    << std::endl;
  
  // calculates and returns the result
  m2C_result=M2C_332Calculator(ltv_Vis_A, 
					          ltv_Vis_B,
					          pT_Miss, 
					           DeltaM);
  
  // displays the result
  std::cout << "M2C Result is " << m2C_result << std::endl;
  
  return 0;
 
}
