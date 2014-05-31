// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

// This program needs standard input to be a file of the form
// such as provided in ../data/20100421-Sakurai.txt
// i.e. of the form of a white space separated list of
// values "ETA pxa pya   ETB pxb pyb    pxmiss pymiss".
// Comment lines (eg those introduced by a "#" character) must be
// stripped out before passing the file to this program (for simplicity).
// On the command line this program takes one or more MT2_332 
// calculator names, and then calculates MT2 using those names, 
// and prints them to std::out in the order they were provided
// on the command line.
// To separate this output from the minuit messages which cannot 
// be suppressed, the output is on lines which start with the 
// word "AlgorithmComparison".

#include "ExampleEvent.h"
#include "Mt2/ChengHanBisect_Mt2_332_Calculator.h"
#include "Mt2/Basic_Mt2_332_Calculator.h"
#include "Mt2/Basic_MCT2_332_Calculator.h"
#include "Mt2/Mt2Calculators.h"
#include <iostream>
#include <string>
#include <list>

void usage() {
  std::cout 
    << "\nUsage example:\n"
    << "       cat ../data/20100421-Sakurai.txt |  grep -v \"^#\" | ./AlgorithmComparison_Example Basic_Mt2_332 ChengHanBisect_Mt2_332 | grep ^AlgorithmComparison | less -S \n\n"
    << "See also the big comment at the top of AlgorithmComparison_Example.cpp\n" << std::endl;
}

int main(int argc, char * argv[]) {
  
  if (argc<=1) {
    usage();
    return 1;
  }
  
  std::list<std::string> calcNames;
  
  for(int i=1; i<argc; ++i) {
    calcNames.push_back(std::string(argv[i]));
  }
  
  double ETA, pxa ,pya ,  ETB, pxb, pyb ,   pxmiss, pymiss;
  
  while (std::cin >> ETA) {
    std::cin >> pxa  >> pya >> ETB >> pxb >> pyb >>  pxmiss >> pymiss;
       

    Mt2::LorentzTransverseVector ltv_Vis_A(/*ETA,*/pxa,pya);
    Mt2::LorentzTransverseVector ltv_Vis_B(/*ETB,*/pxb,pyb);
    Mt2::TwoVector               pT_Miss(pxmiss,pymiss);
    double                      invis_mass = 0;
  
    // Now that we have some visiable momenta and some missing transverse
    // momentum we can calculate MT2.
  
    {
      // get copyright printed here
      static Mt2::ChengHanBisect_Mt2_332_Calculator calculator;
    }

    std::cout << "AlgorithmComparison ";
    for(std::list<std::string>::const_iterator it = calcNames.begin();
	it!= calcNames.end();
	++it) {
      
      const std::string calcName = *it;
      
      if (calcName == "ChengHanBisect_Mt2_332") {
	static Mt2::ChengHanBisect_Mt2_332_Calculator calculator;
	const double mt2 = calculator.mt2_332( ltv_Vis_A, ltv_Vis_B, pT_Miss, invis_mass);
	std::cout << mt2 << " ";
      } else if (calcName == "Basic_Mt2_332") {
	static Mt2::Basic_Mt2_332_Calculator calculator;
	const double mt2 = calculator.mt2_332( ltv_Vis_A, ltv_Vis_B, pT_Miss, invis_mass);
	std::cout << mt2 << " ";
      } else if (calcName == "Basic_MCT2_332") {
	static Mt2::Basic_MCT2_332_Calculator calculator;
	const double mct2 = calculator.mct2_332( ltv_Vis_A, ltv_Vis_B, pT_Miss, invis_mass);
	std::cout << mct2 << " ";
      } else {
	std::cout << "XXX ";
      }
    }
    std::cout << std::endl;

  }

  return 0;
}
