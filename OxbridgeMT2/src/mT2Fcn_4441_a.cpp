// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#include "recipeAUX/OxbridgeMT2/interface/mT2Fcn_4441_a.h"

#include <math.h>
#include <string>

namespace Mt2 {

  template <class T> T sq(const T t) {
    return t*t;
  }

  Mt2::LorentzVector mT2Fcn_4441_a::boostFirstToRFOfSecond(const Mt2::LorentzVector & a,
							 const Mt2::LorentzVector & b) {
    //const Hep3Vector boostVec = -(b.boostVector()); // "-" because we want to take a to rest frame of b, not generate b from its rest frame.
    //return boostOf(a,boostVec); // a CLHEP function
	return a.boostBy(-b.px/b.e,
			 -b.py/b.e,
			 -b.pz/b.e);
  }

  double mT2Fcn_4441_a::trip(const double a,
			     const double b,
			     const double c) {
    return a*a + b*b + c*c - 2.0*(a*b + a*c + b*c);
    // also trip = (A+B+C)*(A-B+C)*(A+B-C)*(A-B-C)
    // where A = sqrt(a), B = sqrt(b) etc ...
  }

double mT2Fcn_4441_a::operator()(const std::vector<double>& x) const {
  
  if (m_debug) {
    std::cout << "pxpypzeAlpha  " << pxpypzeAlpha << std::endl;
    std::cout << "pxpypzeBeta   " << pxpypzeBeta  << std::endl;
    std::cout << "pxpypzeG      " << pxpypzeG     << std::endl;
  }

  const double maSq = pxpypzeAlpha.m2();
  const double mbSq = pxpypzeBeta .m2();

  static const Mt2::LorentzVector pxpypzeLambda(Mt2::LorentzVector::InitEPxPyPz(1,0,0,0));

  if (m_debug) {
    std::cout << "pxpypzeLambda " << pxpypzeLambda << std::endl;
  }

  // Calculate quantities related to F
  const Mt2::LorentzVector pxpypzeF = 
    pxpypzeLambda*rootS - pxpypzeAlpha - pxpypzeBeta -  pxpypzeG;
  
  if (m_debug) {
    std::cout << "pxpypzeF      " << pxpypzeF     << std::endl;
  }

  const double mFSq  = pxpypzeF.m2();
  const double pTFSq = pxpypzeF.pT2();

  //std::cout << "mFSq " << mFSq << " pTFSq " << pTFSq << std::endl;

  // Choose a value for mpq - using method of page (33) of lab book 6.
  const double mpqMinSq = sq(2.0 * mChi);
  const double mpqSq = x[0]*mFSq + (1.0-x[0])*mpqMinSq;
  if (m_debug) {
    std::cout << "mpqMinSq " << mpqMinSq << " mpqSq " << mpqSq << std::endl;
  }

  // Choose a value for mh - using method of page (33) of lab book 6.
  const double mhMaxSq = 2.0*pTFSq + mFSq + mpqSq
    - 2.0*sqrt(sq(pTFSq) + pTFSq*(mFSq+mpqSq) + mFSq*mpqSq);
  const double mhSq = x[1]*mhMaxSq + (1.0-x[1])*0.0;

  if (m_debug) {
    std::cout << "mhMaxSq " << mhMaxSq << " mhSq " << mhSq << std::endl;
    std::cout << "mF-(mh+mpq) " << sqrt(mFSq)-(sqrt(mhSq)+sqrt(mpqSq)) << std::endl;
  }

  // Calculate the vectors that are h and p+q
  const double fz = pxpypzeF.pz;
  const double ef = pxpypzeF.e;

  // Sign should be +1 or -1 in the following ....

  const double hzRadical =  trip(mFSq,mhSq,mpqSq)-4.0*mhSq*pTFSq;
  const double hz = ( fz*(mFSq+mhSq-mpqSq) + sign*ef*sqrt(
							  fabs(hzRadical)
							  ) )/
    (2.0*(  pTFSq + mFSq   ));
  
  if (m_debug) {
    std::cout << " tripargs " << mFSq << " " << mhSq << " " << mpqSq << "\n" 
	      << " trip " << trip(mFSq,mhSq,mpqSq) << "\n"
	      << " fz " << fz << " sign " << sign << " hz " << hz << hzRadical << std::endl;
  }

  const Mt2::LorentzVector & pxpypzeH = Mt2::LorentzVector(Mt2::LorentzVector::InitEPxPyPz(sqrt(sq(hz) + mhSq),0,0,hz));
  const Mt2::LorentzVector pxpypzePQ = pxpypzeF - pxpypzeH;

  if (m_debug) {
    // Check that the mass of (P+Q) is mPQ as it is supposed to be:
    std::cout<< "fractional error on mPQ is "
	     << sqrt(pxpypzePQ.m2()/mpqSq)-1 << std::endl;
    std::cout<< "fractional error on mh is              "
	     << sqrt(pxpypzeH.m2() /mhSq )-1 << std::endl;
  }

  // It is not yet certain that it will be possible to find a direction for p and q adding to p+q such that (alpha+p)^2 == (beta+q)^2   ..... that depends on the value of "cosTheta" (see pg (28) and (29) of labbook 6 one third of the way into the book).

  const Mt2::LorentzVector alphaRF 
    = boostFirstToRFOfSecond(pxpypzeAlpha,pxpypzePQ);
  const Mt2::LorentzVector betaRF
    = boostFirstToRFOfSecond(pxpypzeBeta ,pxpypzePQ);
  const Mt2::LorentzVector pqRF
    = boostFirstToRFOfSecond(pxpypzePQ   ,pxpypzePQ);
  
  if (m_debug) {
    std::cout << "alphaRF_PQ " << alphaRF << std::endl;
    std::cout << " betaRF_PQ " <<  betaRF << std::endl;
    std::cout << "   pqRF_PQ with mpq " << sqrt(mpqSq) << " " << pqRF << std::endl;
  }

  // Now using page (28) of labbook 6:
  const double epRF = sqrt(mpqSq)*0.5; // i.e. mpq/2
  const double ppRF = sqrt(fabs(sq(epRF)-sq(mChi)));  //fabs there only to control rounding errors
  const double eaRF = alphaRF.e;
  const double ebRF =  betaRF.e;
      
  const Mt2::LorentzVector sigmaRF = alphaRF+betaRF;

  const double cosThetaPSigma = (epRF*(eaRF-ebRF)+0.5*(maSq-mbSq))/
                                     (ppRF*sqrt(sigmaRF.p2()));

  if (m_debug) {
    std::cout << "cosThetaPSigma = " << cosThetaPSigma << std::endl;
  }

  if (cosThetaPSigma>+1.0 || cosThetaPSigma<-1.0) {
    if (m_debug) {
      std::cout << "No configuration with (p+a)^2 == (q+b)^2 !";
    }
//#warning NEXT LINE NEEDS IMPROVEMENT!
    return 100000;   // FIXME_PLEASE
  } else {
    if (m_debug) {
      std::cout << "There is a configuration with (p+a)^2 == (q+b)^2" << std::endl;
    }

    const double sinThetaPSigma = sqrt(1.0-sq(cosThetaPSigma));


    const double cosThetaAlphaSigma = alphaRF.cosineOfSpatialSeparationAngleFrom(sigmaRF);
      //cos(alphaRF.vect().angle(  sigmaRF.vect()   ));							 
    const double sinThetaAlphaSigma = sqrt(1.0-sq(cosThetaAlphaSigma));

    if (m_debug) {
      std::cout << "c1,s2,c2,s2 "
		<< cosThetaPSigma << " "
		<< sinThetaPSigma << " "
		<< cosThetaAlphaSigma << " "
		<< sinThetaAlphaSigma << std::endl;
    }

    const double cosDeltaTheta
      = cosThetaPSigma * cosThetaAlphaSigma
      + sinThetaPSigma * sinThetaAlphaSigma;

    // So in this case where (for the given
    // values of mh and mpq) we 
    // can find a form of p for which 
    // m=(a+p)^2==(b+q)^2 then that value of m
    // occurs at the value of m given 
    // on page (28A) of lab book 6:

    const double mslepSq = maSq + sq(mChi) + 2.0*(
		      eaRF*epRF
                      -
                      sqrt(alphaRF.p2())*ppRF*cosDeltaTheta
						  );
    const double f=sqrt(mslepSq);
    if (m_debug) {
      std::cout << "mslep(candidate),x(1),x(2) = "
		<< f << " " 
		<< x[0] << " "
		<< x[1] << std::endl;      
    }
    
    return f;
  }
      
}


}
