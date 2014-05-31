// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr



/* Use it like this:

#include "Mt2/AlphaT_Multijet_Calculator.h"

Mt2::AlphaT_Multijet_Calculator alphaTCalculator;

while (next event ...) {

        alphaTCalculator.clear();

        Mt2::LorentzTransverseVector jet1 = get(somehow);
        Mt2::LorentzTransverseVector jet2 = get(somehow);
        Mt2::LorentzTransverseVector jet3 = get(somehow);
        Mt2::LorentzTransverseVector jet4 = get(somehow);
        Mt2::LorentzTransverseVector jet5 = get(somehow);
        ...
        ...

        alphaTCalculator.push_back(jet1);
        alphaTCalculator.push_back(jet2);
        alphaTCalculator.push_back(jet3);
        alphaTCalculator.push_back(jet4);
        alphaTCalculator.push_back(jet5);
        ...
        ...

        const double alphaT = alphaTCalculator.alphaT();

}

*/

#include "ExampleEvent.h"
#include "Mt2/AlphaT_Multijet_Calculator.h"
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
  Mt2::AlphaT_Multijet_Calculator alphaTCalculator;

  alphaTCalculator.clear();
  
  // Could tell the MTGen calculating object to be verbose, and print out
  // debug messages while it is thinking ... but we won't:
  // mtGenCalculator.setDebug(true);

  Mt2::TwoVector ptMiss;
  //std::vector<Mt2::LorentzTransverseVector> interestingTParticles;
  //std::vector<Mt2::LorentzVector> interestingParticles;
  for (int i=0; i<10; ++i) {
    alphaTCalculator.push_back(ltv_Vis_A);
    alphaTCalculator.push_back(ltv_Vis_B);
  }

  // Now we can actually calculate MT2:
  const double alphaT
    = alphaTCalculator.alphaT();
  
  // Now we print out the result:
  std::cout << "ANSWER: alphaT = " << alphaT 
    //<< std::endl
    //	    << "ANSWER: mtGen = " << mtGen
		   //  << " for " << mtGenCalculator.algorithmName() << " algorithm"
	    << std::endl; 

  alphaTCalculator.clear();
  std::cout << (alphaTCalculator.alphaT()==-1 ? "Test passed!" : "Test failed!" )  << "  (You should get -1 if you attempt to calculate alphaT on no jets.)" << std::endl;

  alphaTCalculator.clear();
  alphaTCalculator.push_back(Mt2::LorentzTransverseVector(5,1,0));
  std::cout << (alphaTCalculator.alphaT()==-1 ? "Test passed!" : "Test failed!" )  << "  (You should get -1 if you attempt to calculate alphaT on a single jet.)" << std::endl;

  alphaTCalculator.clear();
  alphaTCalculator.push_back(Mt2::LorentzTransverseVector(1,+1,0));
  alphaTCalculator.push_back(Mt2::LorentzTransverseVector(1,-1,0));
  std::cout << (alphaTCalculator.alphaT()==0.5 ? "Test passed!" : "Test failed!" )  << "  (You should get 0.5 if you attempt to calculate alphaT on a perfectly balanced back-to-back QCD event.)" << std::endl;


  alphaTCalculator.clear();
  alphaTCalculator.push_back(Mt2::LorentzTransverseVector(1,+1,0));
  alphaTCalculator.push_back(Mt2::LorentzTransverseVector(0.2,-0.2,0));
  std::cout << (alphaTCalculator.alphaT()<=0.5 ? "Test passed!" : "Test failed!" )  << "  (You should get less than 0.5 if you attempt to calculate alphaT on an imbalanced back-to-back QCD event.  We got "<<alphaTCalculator.alphaT()<<")" << std::endl;

  alphaTCalculator.clear();
  alphaTCalculator.push_back(Mt2::LorentzTransverseVector(0.2,+0.2,0));
  alphaTCalculator.push_back(Mt2::LorentzTransverseVector(1,-1,0));
  std::cout << (alphaTCalculator.alphaT()<=0.5 ? "Test passed!" : "Test failed!" )  << "  (You should get less than 0.5 if you attempt to calculate alphaT on an imbalanced back-to-back QCD event.  We got "<<alphaTCalculator.alphaT()<<")" << std::endl;

  alphaTCalculator.clear();
 alphaTCalculator.push_back(Mt2::LorentzTransverseVector(5,3,4));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(13,5,12));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(25,7,24));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(17,8,15));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(41,9,40));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(61,11,60));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(37,12,35));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(85,13,84));
std::cout << "For example A massless jet challenge sent to Pascal Nef alphaT is " << alphaTCalculator.alphaT() << " and Pascal in a mail of 2nd Aug 2010 (talks) says he gets 4.56027 " << std::endl;

  alphaTCalculator.clear();
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(7,3,4));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(15,5,12));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(27,7,24));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(19,8,15));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(43,9,40));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(63,11,60));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(39,12,35));
alphaTCalculator.push_back(Mt2::LorentzTransverseVector(87,13,84));
std::cout << "For example B massive  jet challenge sent to Pascal Nef alphaT is " << alphaTCalculator.alphaT() << " and Pascal in a mail of 2nd Aug 2010 (talks) says he gets 4.56027 " << std::endl;
 
  return 0;
}
