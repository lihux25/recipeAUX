// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "recipeAUX/OxbridgeMT2/interface/MCTll_332_Calculator.h"

namespace Mt2 {



  double mctll_332  (const LorentzTransverseVector& visibleA,  // 3 d.o.f. 
		     const LorentzTransverseVector& visibleB,  // 3 d.o.f.
                     const TwoVector& utm){
    const double msq = mctll_332_Sq(visibleA, visibleB, utm);
    if (msq<=0) {
      return 0;
    } else {
      return sqrt(msq);
    }
  }

  double mctll_332_Sq(const LorentzTransverseVector& a,  // 3 d.o.f. 
		      const LorentzTransverseVector& b,  // 3 d.o.f.
                      const TwoVector& utm) {

    const ResolvedLTV aResolved(a,utm);
    const ResolvedLTV bResolved(b,utm);

    return a.masssq() + b.masssq() + 2.0*(aResolved.para()).contralinearDot(bResolved.para());
  } 

}
