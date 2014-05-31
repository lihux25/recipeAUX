// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "Mt2/Basic_Mt2_332_Calculator.h"
#include <iostream>

int main(int argc, char * argv[]) {

  // First we create the object that is going to do the calculation
  // of MT2 for us.  You can do this once early on, and re-use it
  // multiple times.
  //
  // For this example we will use the "Basic_Mt2_332_Calculator" which is
  // the algorithm we recommend people use by default.
  Mt2::Basic_Mt2_332_Calculator mt2Calculator;
  
  // Could tell the MT2 calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  
  // mt2Calculator.setDebug(true);
  
  // Now we can actually calculate MT2.  Let's do it three times for the same "inputs".

  for (int i=0; i<3; ++i) {

    // The input parameters associated with the particle
    // (or collection of particles) associated with the 
    // first "side" of the event:
    const double massOfSystemA =  100; // GeV
    const double pxOfSystemA   =  410; // GeV
    const double pyOfSystemA   =   20; // GeV
    
    // The input parameters associated with the particle
    // (or collection of particles) associated with the 
    // second "side" of the event:
    const double massOfSystemB =  150; // GeV
    const double pxOfSystemB   = -210; // GeV
    const double pyOfSystemB   = -300; // GeV
    
    // The missing transverse momentum:
    const double pxMiss        = -200; // GeV
    const double pyMiss        =  280; // GeV

    // The mass of the "inivisible" particle presumed to have
    // been produced at the end of the decay chain in each
    // "half" of the event:
    const double invis_mass    = 100; // GeV
    
    // Now put the inputs together into the input structures that the library wants.
    
    /*
      Note: in the next two lines (constructing "vis_A" and "vis_B"),
      the ORDER of the arguments to the constructor of
      Mt2::LorentzTransverseVector is very important. 
      You need to be careful as, when the TwoVector comes
      first, the second arguments is taken to be a mass:
      
      LorentzTransverseVector(const TwoVector& momentum, double mass);
      
      but when the TwoVector comes second, the first arguemt is an ET=Sqrt(m^2+pt^2):
      
      LorentzTransverseVector(double Et, const TwoVector& momentum);
      
      You have been warned!
    */
    
    Mt2::LorentzTransverseVector  vis_A(Mt2::TwoVector(pxOfSystemA, pyOfSystemA), massOfSystemA);
    Mt2::LorentzTransverseVector  vis_B(Mt2::TwoVector(pxOfSystemB, pyOfSystemB), massOfSystemB);
    Mt2::TwoVector                pT_Miss(pxMiss, pyMiss);
    
    std::cout << "Going to calculate MT2 with\n"
	      << "   ltv_Vis_A  = " << vis_A  << "\n"
	      << "   ltv_Vis_B  = " << vis_B  << "\n"
	      << "   pT_Miss    = " << pT_Miss    << "\n"
	      << "   invis_mass = " << invis_mass << std::endl;
    
    // Now that we have some visiable momenta and some missing transverse
    // momentum we can calculate MT2.
    
    const double mt2 
      = mt2Calculator.mt2_332(vis_A, vis_B, pT_Miss, invis_mass);
    
    // Now we print out the result:
    std::cout << "ANSWER: mt2 = " << mt2 
	      << " for " << mt2Calculator.algorithmName() << " algorithm"
	      << std::endl; 
    
  }
  return 0;
}
