// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#include "recipeAUX/OxbridgeMT2/interface/MC_330_Calculator.h"

namespace Mt2 {

  double mc_330   (const LorentzVector& a, 
		    const LorentzVector& b) {
    const double mcsq = mc_330_Sq(a, b);
    // if less than zero, this is only evidence of inexact maths:
    if (mcsq<=0) {
      return 0;
    } else {
      return sqrt(mcsq);
    }
  }
  
  double mc_330_Sq(const LorentzVector& a,  
		    const LorentzVector& b) {
    return a.masssq() + b.masssq() + 2.0*a.contralinearDot(b);
  }

}
