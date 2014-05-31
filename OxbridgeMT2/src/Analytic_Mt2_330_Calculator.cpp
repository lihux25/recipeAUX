// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "recipeAUX/OxbridgeMT2/interface/Analytic_Mt2_330_Calculator.h"
#include "recipeAUX/OxbridgeMT2/interface/Mt2Units.h"
#include "recipeAUX/OxbridgeMT2/interface/Analytic_Mt2_332_Assistor.h"

namespace Mt2 {

  /*
    
  */


  double Analytic_Mt2_330_Calculator::mt2_330(const LorentzTransverseVector& visA, 
					      const LorentzTransverseVector& visB,
					      const double mChi) {
     return sqrt(mt2_330_Sq(visA, visB, mChi));
  }

  double Analytic_Mt2_330_Calculator::mt2_330_Sq(const LorentzTransverseVector& visA, 
					         const LorentzTransverseVector& visB,
					         const double mChi) {

    m_notes.clear();

    static const LorentzTransverseVector lambda(1,0,0);

    const LorentzTransverseVector sigma = visA + visB;
    const LorentzTransverseVector delta = visA - visB;

    const double mASq = visA.masssq();
    const double mBSq = visB.masssq();
    const double mChiSq = mChi*mChi;
    //const double mA = sqrt(mASq);
    //const double mB = sqrt(mBSq);

    // first rule out or use unbalanced solutions:
    {
      if (debug()) {std::cout << "Looking for unbalanced solutions." << std::endl;}
      static Analytic_Mt2_332_Assistor assistor;
      const Analytic_Mt2_332_Assistor::Ans 
	assistance = assistor.assist(visA,visB,-(sigma.transverse()),mChi);
      
      if (assistance.solutionType.type == ASideGlobalMin ||
	  assistance.solutionType.type == BSideGlobalMin) {

	if(debug()) {std::cout << "Found unbalanced solution at " << (assistance.solutionType.type==ASideGlobalMin?"A":"B")<<" side unconstrained min." << std::endl;}

	m_notes.soln       = assistance.solutionType;
	m_lastSolutionType = assistance.solutionType;
	return (assistance.mt2)*(assistance.mt2);
      }
    }


    // now look for ballanced solutions
    m_notes.soln       = SolutionType(BalancedSolution);
    m_lastSolutionType = SolutionType(BalancedSolution);
    
    if(debug()) {std::cout << "Looking for BALANCED solution" << std::endl;}

    const double lambdaDotSigma = lambda.dot(sigma);
    const double lambdaDotDelta = lambda.dot(delta);
    const double sigmaDotDelta  = sigma.dot(delta);
    //const double lambdaDotSigmaSq = lambdaDotSigma * lambdaDotSigma;
    //const double sigmaSq = sigma.masssq();


    // want to minimise    A x - lambda sqrt(  (x-B)^2 - C^2 )
    const double ATERM = 0.5*lambdaDotSigma
      + 0.5*lambdaDotDelta*(sigmaDotDelta-lambdaDotSigma*lambdaDotDelta)/(sigma.ptsq());
      
    const double J = (sigma.ptsq()-lambdaDotDelta*lambdaDotDelta)/(4.*(sigma.ptsq()));
    
#ifndef __CINT__
//#warning Should check/think about negative J possibility!
#endif
    const double LAMBDATERM = sqrt(J)*(  
			       delta.px()*sigma.py() - delta.py()*sigma.px()
			       )/sigma.pt();

    const double BTERM = lambdaDotSigma;
#ifndef __CINT__
//#warning Should check/think about negative C possibility!
#endif
    const double CTERM = sqrt(sigma.ptsq() + mChiSq/J);

    //for (unsigned int i=0; i<=1; ++i) {


      // For the standard event, a value of 1572 ish is needed for rootS.
      
      const double rootS =  BTERM + CTERM / sqrt(1.0 - LAMBDATERM*LAMBDATERM/(ATERM*ATERM)) ;

      m_notes.rootS = rootS;
      /*
      std::cout << "BT\t" << BTERM << "\t"
		<< "CTERM\t" << CTERM << "\t"
		<< "LT\t" << LAMBDATERM << "\t"
		<< "AT\t" << ATERM << "\t"
		<< "rootSMin\t" << lambdaDotSigma + sqrt(4.*mChiSq + sigma.ptsq())
			<< std::endl;
      */

    // Technically we have:
    // const double rootSMin = lambdaDotSigma + sqrt(4.*mChiSq + lambdaDotSigmaSq - sigmaSq);
    // but this is the same as the simpler:
    //const double rootSMin = lambdaDotSigma + sqrt(4.*mChiSq + sigma.ptsq());

    //for (double extraBit = 0.1; extraBit<4000; extraBit += .1) { 

      double mt2Sq = 0;



      {

      //const double rootS = rootSMin + extraBit;
      //const double rootS = 1572;
      const LorentzTransverseVector bigB = lambda*rootS - sigma;
      const double sigmaDotBigB = sigma.dot(bigB);
      const double bigBSq = bigB.masssq();
      
      const LorentzTransverseVector w = LorentzTransverseVector(sigma.px()*bigB.py() - sigma.py()*bigB.px(),
			-sigma.py()*bigB.Et() + sigma.Et()*bigB.py(), // strange signs from index lowering
			-sigma.Et()*bigB.px() + sigma.px()*bigB.Et() // strange signs from index lowering
			);

      //std::cout << sigma << " " << bigB << " " << w << " " << w.dot(sigma) << " " << w.dot(bigB) << " " << w.masssq() << std::endl;
      
      const double wSq = w.masssq(); 
      const double minusWSq = -wSq; // should be positive and equal to (sigma.dot(bigB))^2 - sigmaSq * bigBSq
      //std::cout << minusWSq << "\t" << pow(sigma.dot(bigB),2) - sigmaSq * bigBSq << std::endl;
      
      const LorentzTransverseVector wHat = w/sqrt(minusWSq);
      const LorentzTransverseVector H = (sigma*bigBSq - bigB*sigmaDotBigB)*((-0.5*(delta.dot(bigB+sigma)))/wSq);
      const double HSq = H.masssq();

      const LorentzTransverseVector secondHalfGamma = wHat * sqrt(0.25*(bigBSq-4.*mChiSq) + HSq);
      
      const LorentzTransverseVector gamma_1 = H + secondHalfGamma;
      const LorentzTransverseVector gamma_2 = H - secondHalfGamma;

      const double firstPart = mChiSq + 0.5*(mASq+mBSq)+0.5*sigmaDotBigB;

      const double mt2Sq_1 = firstPart + delta.dot(gamma_1);
      const double mt2Sq_2 = firstPart + delta.dot(gamma_2);

      if (mt2Sq_1 > mt2Sq_2) {
	mt2Sq = mt2Sq_2;
        m_notes.gamma = gamma_2;
      } else {
	mt2Sq = mt2Sq_1;
	m_notes.gamma = gamma_1;
      }

      if(debug()) {std::cout << "rootS " << rootS << " MT2 " << sqrt(mt2Sq) << std::endl;}
      // I have checked that for the standard event which has MT2 = 412.628 the above algorith does indeed give 412.628 for a value of rootS of between 1570.51  and 1575.51

      // All that remains is to find the correct value of rootS to substitite in.

      /*
	std::cout << rootS << " =\t" << mt2Sq << "\t" 
	<< sqrt(mt2Sq_1) << "\t" << sqrt(mt2Sq_2) << "\t"
	<< mt2Sq_1 << "\t" << mt2Sq_2 << "\t"
	<< std::endl; 

      */
      

      }

      return mt2Sq;
  }
  

} // end of Mt2 Namespace
