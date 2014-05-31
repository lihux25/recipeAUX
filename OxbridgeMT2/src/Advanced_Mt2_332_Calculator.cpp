// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "recipeAUX/OxbridgeMT2/interface/Analytic_Mt2_332_Assistor.h"
#include "recipeAUX/OxbridgeMT2/interface/Advanced_Mt2_332_Calculator.h"
#include "recipeAUX/OxbridgeMT2/interface/Basic_Mt2_332_Calculator.h"

namespace Mt2 {

  double Advanced_Mt2_332_Calculator::mt2_332(const LorentzTransverseVector& visA, 
		    			      const LorentzTransverseVector& visB,
					      const TwoVector& ptmiss, 
					      const double mEachInvisible) {
    const double mSq = mt2_332_Sq(visA, visB, ptmiss, mEachInvisible);

    if (mSq>=0) {
      return sqrt(mSq);
    } else {
      // something went wrong -- maybe just rounding errors;
      return sqrt(fabs(mSq));
    }
  }


  double Advanced_Mt2_332_Calculator::mt2_332_Sq(const LorentzTransverseVector& visA,
					         const LorentzTransverseVector& visB,
					         const TwoVector& ptmiss,
					         const double mEachInvisible){

    m_lastSolutionType = SolutionType(NotSpecified);

    const bool useFancyAssistor = m_useFancyAssistor; // You can enable this if you like, but the fancy assistor is still a bit naughty about how it deals with the limit of massless visible particles, so the fancy assistor is disabled by default.  Lester (in light of Mario's and Alan's (correct) comments below).
 
// Old comment from Mario/Alan: It was thought by C.Lester that false and true are both valid.  true selects fast evaluation of unbalanced cases.  false brute forces the calc.
// M.Serna found cases (see FunnyEvent.cpp in bin) in which true is WRONG (leads to discontinutites in MT2).
// So A.Barr changes the above line to false 5 May 08.
// August 2008: Lester found the cause of the fancy assistor bug, and fixed it.  However he agrees that the fancy assistor should not be used in Basic_332, and so has moved it out to Advanced_332 (where it is still disabled by default!!) so that it will only be used if people really want it.

    if (useFancyAssistor) {
      static Analytic_Mt2_332_Assistor assistor;
      const Analytic_Mt2_332_Assistor::Ans 
	assistance = assistor.assist(visA,visB,ptmiss,mEachInvisible);
      
      if (assistance.solutionType.type == ASideGlobalMin ||
	  assistance.solutionType.type == BSideGlobalMin) {
	m_lastSolutionType = assistance.solutionType;
	return (assistance.mt2)*(assistance.mt2);
      }
    }
    m_lastSolutionType = SolutionType(NotSpecified);

    static Basic_Mt2_332_Calculator calculator;
    return calculator.mt2_332_Sq(visA, visB, ptmiss, mEachInvisible);
  }
  
} // end of Mt2 Namespace
