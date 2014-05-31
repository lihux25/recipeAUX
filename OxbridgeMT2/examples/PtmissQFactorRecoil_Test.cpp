// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "ExampleEvent.h"
#include "Mt2/Advanced_Mt2_332_Calculator.h"
#include <iostream>

inline double sq(const double x) { return x*x; }

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
 
  while(true){

    exampleEvent.resetRandom332();
    
    Mt2::LorentzTransverseVector ltv_Vis_A = exampleEvent. ltv_Vis_A();
    Mt2::LorentzTransverseVector ltv_Vis_B = exampleEvent. ltv_Vis_B();
    const double Q = 4.*(exampleEvent.flipshoot());
    Mt2::TwoVector                 visSum = (ltv_Vis_A.transverse() + ltv_Vis_B.transverse());
    Mt2::TwoVector                 pT_Miss = visSum * Q;
    double                      chi = exampleEvent.invis_mass();
    
    std::cout << "Going to calculate MT2 with\n"
	      << "   ltv_Vis_A  = " << ltv_Vis_A  << "\n"
	      << "   ltv_Vis_B  = " << ltv_Vis_B  << "\n"
	      << "   pT_Miss    = " << pT_Miss    << "\n"
	      << "   invis_mass = " << chi << std::endl;
    
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
    Mt2::Advanced_Mt2_332_Calculator mt2Calculator(true);
    
    // Could tell the MT2 calculating object to be verbose, and print out
    // debug messages while it is thinking ... but we won't:
    
    // mt2Calculator.setDebug(true);
    
    // Now we can actually calculate MT2:
    const double mt2 
      = mt2Calculator.mt2_332( ltv_Vis_A, ltv_Vis_B, pT_Miss, chi);

    
    Mt2::SolutionType type = mt2Calculator.lastSolutionType();

    const double ma = ltv_Vis_A.mass();
    const double mb = ltv_Vis_B.mass();
    const double AT = ltv_Vis_A.contralinearDot(ltv_Vis_B);

    /* The following is correct for Q=0 or Q=1 or ma=mb, and has "tuned" powers (the second (1.+Q) term and teh pow(fabs(1+Q),2) term) to make it do "well" elsewhere. */

    const double other =

      chi*chi 
      + (1.+Q)*(ma*ma + mb*mb)/2.  
      +   (1.+Q)*sq(mb*mb-ma*ma) / (2.*(2.*AT - ma*ma - mb*mb))
      - AT*Q
+
      sqrt(
	   (AT*AT - ma*ma *mb*mb)
	   *
	   (      

	    sq((mb*mb-ma*ma)/(2.*AT-ma*ma-mb*mb))*pow(fabs(1.+Q),2) 
	    +
	    (4.* chi*chi)/(2.* AT - ma*ma - mb*mb)
	    +
	    Q*Q
	    )
	   );

    // Now we print out the result:
    std::cout << "ANSWER: mt2 = " << mt2 
	      << " for " << mt2Calculator.algorithmName() << " algorithm"
	      << std::endl; 
    std::cout << "COMPARE " << type << " " << mt2 << " " << sqrt(other) << " " << Q << " " << 2.*AT - ma*ma - mb*mb <<" " << AT*AT - ma*ma *mb*mb << std::endl;
  }
  return 0;
}
