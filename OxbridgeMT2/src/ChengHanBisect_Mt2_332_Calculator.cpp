// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "recipeAUX/OxbridgeMT2/interface/ChengHanBisect_Mt2_332_Calculator.h"
#include "recipeAUX/OxbridgeMT2/interface/Mt2Units.h"

#include "mt2_bisect.h"

namespace Mt2 {

  ChengHanBisect_Mt2_332_Calculator::ChengHanBisect_Mt2_332_Calculator() : Mt2_332_Calculator("ChengHanBisect_Mt2_332") {
  }

  double ChengHanBisect_Mt2_332_Calculator::mt2_332_Sq(const LorentzTransverseVector& visibleA, // 3 d.o.f. 
						 const LorentzTransverseVector& visibleB, // 3 d.o.f.
						 const TwoVector& ptmiss,                 // 2 d.o.f.
						 const double mInvisible){

    double mt2=ChengHanBisect_Mt2_332_Calculator::mt2_332(visibleA, visibleB, ptmiss, mInvisible);
    return mt2*mt2;
  }
  

  double ChengHanBisect_Mt2_332_Calculator::mt2_332(const LorentzTransverseVector& visA, 
							       const LorentzTransverseVector& visB,
					      const TwoVector& ptmiss, 
					      const double mEachInvisible) {


    double pa[3] = { visA.mass(), visA.px(), visA.py() };
    double pb[3] = { visB.mass(), visB.px(), visB.py() };
    double pmiss[3] = { 0, ptmiss.px(), ptmiss.py() };
    double mn = mEachInvisible;

    mt2_bisect::mt2 mt2_event;

    mt2_event.set_momenta(pa,pb,pmiss);
    mt2_event.set_mn(mn);

    return mt2_event.get_mt2();
  }

} // end of Mt2 Namespace
