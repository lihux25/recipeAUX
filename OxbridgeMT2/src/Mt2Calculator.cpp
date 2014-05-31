// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#include "recipeAUX/OxbridgeMT2/interface/Mt2Calculator.h"
#include <sstream>
#include <iomanip>
#include <iostream>

namespace Mt2 {

  /*

  /// The Mt2Calculator takes ownership of the minimiser whose address is passed to it on construction.  That minimiser will be deleted when Mt2Calucation::~Mt2Calculator() (the destructor) is called.
  Mt2Calculator::Mt2Calculator(Mt2Minimiser* minimiser) :
  rootSMin(0),
  rootSMax(root_s_max),
  m_minimiser(minimiser)
  {
  initialise();
  }

  Mt2Calculator::~Mt2Calculator(){
  delete m_minimiser;
  }
  */

  bool Mt2Calculator::m_initialised = false;

  Mt2Calculator::Mt2Calculator(const std::string & algName) : m_debug(false), m_algName(algName) {
    if (!m_initialised){
      m_initialised = true;
      std::cout << "----------------------------------------------------------------------\n";
      std::cout << " M_T2 : a variable for measuring masses when missing energy is expected.\n";
      std::cout << " If you use MT2 or this library, please cite:\n";
      std::cout << " \t(o) C.G.Lester, D.J.Summers \n"
		<< " \t    Phys.Lett.B.463:99-103 (1999) hep-ph/9906349\n";
      std::cout << " \t(o) A.J.Barr, C.G.Lester, P.Stephens \n"
		<< " \t    J.Phys.G 29:2343-2363  (2003) hep-ph/0304226\n";
      std::cout << " If you use MTGEN please also cite:\n"
	        << " \t(o) C.G.Lester, A.J.Barr (2007) arXiv:0708.1028 (hep-ph)\n"
		<< " \t    JHEP 12(2007)102\n";
      std::cout << " If you use the ChengHanBisection MT2 calculator, please also cite:\n"
                << " \t(o) H.Cheng, Z.Han,  arXiv:0810.5178\n";
      std::cout << "----------------------------------------------------------------------\n";
    }
    m_lastSolutionType = SolutionType(NotSpecified);
    /*
      m_calls=0;
      m_evaluations=0;
      m_unbalancedEvaluations=0;
  
      m_initialised=true;
    */
  }

  template<class T> 
  T max(T first, T second){
    return first > second ? first : second;
  }

  template<class T> 
  T min(T first, T second){
    return first < second ? first : second;
  }



#ifdef USE_SUSPECT_CALCLUATOR
  double Suspect_Mt2_332_Calculator::mt2_332(const LorentzTransverseVector& visA, 
					     const LorentzTransverseVector& visB,
					     const TwoVector& ptmiss, 
					     double mInvisibles){

    TwoVector other = - ptmiss - visA.transverse() - visB.transverse();
    double energy = visA.Et() + visB.Et() + other.pt();
    LorentzTransverseVector totalVisible(energy, -ptmiss);
  
    return mt2(visA, visB, totalVisible, mInvisibles);
  }

  double Mt2Calculator::mt2(const LorentzTransverseVector& visA, 
			    const LorentzTransverseVector& visB,
			    const LorentzTransverseVector& totalVisible, 
			    double mInvisibles){
    m_calls++;
  
    if (debug()){
      std::cout << "M_T2 calculation inputs : \n"
		<< "alpha: " << visA << "\n" 
		<< "beta:  " << visB << "\n"
		<< "total visible: " << totalVisible << "\n"
		<< "mass of invisible = " << mInvisibles << std::endl;
    }

    /// cache lots of useful information, so it dosent have to be calculated every time.
    alpha    = visA;
    beta     = visB;
    bigSigma = totalVisible;
  
    // derived quantities
    mlspsq  = square(mInvisibles);
    sigma   = alpha + beta;
    delta   = alpha - beta;
    lambda  = LorentzTransverseVector(1,0,0);
  
    // calculate dot products
    lambdaDotBigSigma = lambda.dot(bigSigma);
    lambdaDotSigma    = lambda.dot(sigma);
    lambdaDotDelta    = lambda.dot(delta);
    deltaDotBigSigma  = delta.dot(bigSigma);
    deltaDotSigma     = delta.dot(sigma);
    sigmaDotBigSigma  = sigma.dot(bigSigma);
  
    // invariants
    bigSigmaSq        = bigSigma.masssq();
    sigmaSq           = sigma.masssq();
    deltaSq           = delta.masssq();
    alphaSqPlusBetaSq = alpha.masssq() + beta.masssq();

    // minimum invariant mass = total visible + pseudo particle from 2*lsps
    rootSMin = lambdaDotBigSigma + sqrt(4.*mlspsq 
					+ square(lambdaDotBigSigma)
					- bigSigmaSq);
    if (debug()) printVariables();

    return getMinimiser()->minimise(*this);
  }

  /*
    calculate Mt2 constrained at a particular value of rootS
  */
  double Mt2Calculator::trialMt2(double contribution) const {
    m_evaluations++;

    double rootS = rootSMin + contribution;

    if (debug()) std::cout << "rootS = " << rootS << " ";

    // B is the lorentz 1+2 vector of the sum of the missing particles 
    double BSQ = square(rootS) + bigSigmaSq - 2.*rootS*lambdaDotBigSigma;
    double sigmaDotB = lambdaDotSigma*rootS - sigmaDotBigSigma;

    double H = square(sigmaDotB) - sigmaSq*BSQ;
    double deltaDotB = lambdaDotDelta*rootS - deltaDotBigSigma;
    
    double tt1 = square(deltaDotB) * sigmaSq
      + square(lambdaDotSigma) * BSQ
      - 2.*deltaDotB*deltaDotSigma*sigmaDotB;
  
    double tt2 = square(deltaDotB + deltaDotSigma);
  
    double inner1 = (BSQ-4.*mlspsq)*H - BSQ*tt2;
    double inner2 = H*(-deltaSq) - tt1;
  
    // case where cos theta^2 > 1 so no balanced Pt solution
    if (inner1<0.){
      m_unbalancedEvaluations++;
      double rootBSQ = sqrt(BSQ);
      double eA = 0.5 * (sigmaDotB + deltaDotB)/rootBSQ;
      double eB = 0.5 * (sigmaDotB - deltaDotB)/rootBSQ;
      double pA = sqrt(fabs(square(eA) - alpha.masssq()));
      double pB = sqrt(fabs(square(eB) - beta.masssq()));
      double eChi = 0.5 * rootBSQ;
      double pChi = sqrt(fabs(square(eChi)-mlspsq));
      double m1sq = square(eChi+eA)-square(pChi+pA);
      double m2sq = square(eChi+eB)-square(pChi+pB);
      double msq = m1sq > m2sq ? m1sq : m2sq;
      if (debug()) std::cout << "U : " << sqrt(msq) << std::endl;
    
      return sqrt(msq);
    }
    if (inner2<0.) inner2=0.;
  
    double mt2sq = mlspsq 
      + 0.5*(alphaSqPlusBetaSq) 
      + 0.5*sigmaDotB 
      - (0.5/H) * (
		   (deltaDotB + deltaDotSigma) * 
		   (deltaDotB*sigmaDotB - BSQ*deltaDotSigma)
		   + sqrt(inner1*inner2)
		   );
    if (mt2sq<0.){
      return rootSMin - mt2sq*exp(-mt2sq);
    }

    if (debug()) std::cout << "M : " << sqrt(mt2sq) << std::endl;

    return sqrt(mt2sq);
  }

  // minimise using the nag c-library function 
  // nag_opt_one_var_no_deriv (e04abc)
  // http://www.nag.co.uk/numeric/CL/nagdoc_cl08/xhtml/E04/e04abc.xml

  Mt2Minimiser* Mt2Calculator::getMinimiser(){
    if (!m_minimiser) throw Mt2Exception("No minimiser defined!");
    return m_minimiser;
  }

  void Mt2Calculator::setDebug(bool debug) { 
    m_debug = debug; 
  }

  bool Mt2Calculator::debug() const{ 
    return m_debug; 
  }

  unsigned Mt2Calculator::calls(){
    return m_calls;
  }

  unsigned Mt2Calculator::evaluations(){
    return m_evaluations;
  }

  unsigned Mt2Calculator::unbalancedEvaluations(){
    return m_unbalancedEvaluations;
  }

  void Mt2Calculator::printVariables() const{
    using namespace std;
    cout << "Vectors:"
	 << "\n alpha    : " << setprecision(10) << alpha
	 << "\n beta     : " << setprecision(10) << beta
	 << "\n bigSigma : " << setprecision(10) << bigSigma
	 << "\n sigma    : " << setprecision(10) << sigma
	 << "\n delta    : " << setprecision(10) << delta
	 << "\n lambda   : " << setprecision(10) << lambda
	 << "\n" << endl;
  
    cout.width(10);
    cout << "Variables:" 
	 << "\n lam.bigSig : " << setprecision(10) << lambdaDotBigSigma
	 << "\n lam.sig    : " << setprecision(10) << lambdaDotSigma
	 << "\n lam.del    : " << setprecision(10) << lambdaDotDelta
	 << "\n del.bigSig : " << setprecision(10) << deltaDotBigSigma
	 << "\n del.sig    : " << setprecision(10) << deltaDotSigma
	 << "\n sig.bigSig : " << setprecision(10) << sigmaDotBigSigma      
	 << "\n bigSigSq   : " << setprecision(10) << bigSigmaSq
	 << "\n sigSq      : " << setprecision(10) << sigmaSq
	 << "\n delSq      : " << setprecision(10) << deltaSq
	 << "\n a^2+b^2    : " << setprecision(10) << alphaSqPlusBetaSq
	 << "\n sqrt(smin) : " << setprecision(10) << rootSMin 
	 << endl;
  }

#endif

}
