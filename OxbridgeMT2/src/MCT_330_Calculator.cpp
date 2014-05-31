// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#include "recipeAUX/OxbridgeMT2/interface/MCT_330_Calculator.h"

namespace Mt2 {

  double mct_330   (const LorentzTransverseVector& a, 
		    const LorentzTransverseVector& b) {
    const double mctsq = mct_330_Sq(a, b);
    // if less than zero, this is only evidence of inexact maths:
    if (mctsq<=0) {
      return 0;
    } else {
      return sqrt(mctsq);
    }
  }
  
  double mct_330_Sq(const LorentzTransverseVector& a,  
		    const LorentzTransverseVector& b) {
    return a.masssq() + b.masssq() + 2.0*a.contralinearDot(b);
  }

   
}
