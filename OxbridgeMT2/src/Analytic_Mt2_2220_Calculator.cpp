// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "recipeAUX/OxbridgeMT2/interface/Analytic_Mt2_2220_Calculator.h"
#include "recipeAUX/OxbridgeMT2/interface/Mt2Units.h"
#include "Quartic.h"

#include <math.h>
#include <string>
#include <set>
#include <complex>

namespace Mt2 {

  void Analytic_Mt2_2220_Calculator::insertIfOK(const std::complex<double> & c,
						std::vector<double> & goodSinValues) {
    bool isOk = false;

    if (c.imag()==0 || fabs(c.imag())<1.0E-10) {
      isOk = true;
    } else if (c.real() != 0 && fabs(c.imag()/c.real())<1.E-6) {
      isOk = true;
    }

    if (isOk) {
      goodSinValues.push_back(c.real());
    }
     
  }


  /*
    Now the Analytic_Mt2_2220_Calculator produces answers which (for many millions of randomly generated events) are all as good as, or better than, those results from Basic_Mt2_332, when applied to massless visibles at chi=0.  Note that this does not exclude the possibility that Analytic_Mt2_2220_Calculator has bugs or instabilities under special case inputs (eg a parallel to b woudl probably be bad, as would many other alignments between simple sums of input momenta.
  */


  /* Documentation note:

  Much of the calcualtions that support this method are to be found in mathematica notebooks produced on/around 2010/06/16 - 2010/06/24.  In particular, 

  /usera/lester/proj/Ma/20100616a-interactivegraphics.nb

  and

  /usera/lester/proj/Ma/20100624a-mt2-masslessvisiblesandmasslesschi-lookingfornicetriganswers-dumpdotprods-3-correced.nb

  In addition, the main pen hand/working (which details how the values of Nss, Ncc, Ns, Nc, Nsc and N1 are determined) are to be found stapled approximately 30 sides into the lab-book on Bayesian inference.
  */

  double Analytic_Mt2_2220_Calculator::mt2_2220(const TwoVector& visA, 
						const TwoVector& visB,
						const TwoVector& ptmiss){
    
    // WARNING -- we really need a different prototype as we IGNORE the input value of mEachInvisible ... as it must be zero here ... it is a requirement of the method.

    return sqrt(fabs(Analytic_Mt2_2220_Calculator::mt2_2220_Sq(visA,visB,ptmiss)));
  }
  
  double Analytic_Mt2_2220_Calculator::mt2_2220_Sq(const TwoVector& visA, 
						   const TwoVector& visB,
						   const TwoVector& ptmiss){{
    
      // WARNING -- FOR TESTING PURPOSES I AM HIDING THE INPUTS!
    
      /* For testing, the following example matches 20100622a-mt2 series mathematica notebooks */

      /*    const TwoVector ptmiss(1.0,0.0);
	    const TwoVector visA(TwoVector(cos(1.0),sin(1.0))*1.2);
	    const TwoVector visB(TwoVector(cos(2.1),sin(2.1))*0.8);
      */
   

      m_info = Info();
 

      const double a = visA.pt();
      const double b = visB.pt();
      const double magPtmiss = ptmiss.pt();

      const TwoVector movedP(-(ptmiss.py()),+(ptmiss.px()));
      const TwoVector movedB(-(visB.py()),+(visB.px()));
 
 
      const double adotp = visA.dot(ptmiss);
      const double bdotp = visB.dot(ptmiss);
      const double ahdotbh = visA.dot(visB)/sqrt((visA.ptsq())*(visB.ptsq()));

      const double eap = visA.dot(movedP);
      const double ebp = visB.dot(movedP);
      const double eahbh = visA.dot(movedB)/sqrt((visA.ptsq())*(movedB.ptsq()));

      // For the special case of ma=0, mb=0, chi=0, MT2 is *zero* whenever the ptmiss vector is "between" visA and visB, since this always admits a solution in which p and q are parallel to a and b (respectively) making MT_a = MT_b = 0.  So find coeffs of ptmiss in the a,b basis, and see if they are both positive.  If they are, return zero:

      if (eap/eahbh >= 0 &&  -ebp/eahbh >=0) {
	return 0;
      }

      if (m_debugMode) {
	std::cout << "example=({" 
		  << "  a -> " << a
		  << ", b -> " << b
		  << ", p -> " << magPtmiss
		  << ", thetaap -> ArcTan[" << adotp << "," << eap << "]"
		  << ", thetabp -> ArcTan[" << bdotp << "," << ebp << "]"
		  << "})/.correct" << std::endl;
      
	std::cout << "adotp eap " << adotp << " " << eap << std::endl;
	std::cout << "bdotp ebp " << bdotp << " " << ebp << std::endl;
	std::cout << "ahdotbh eahbh " << ahdotbh << " " << eahbh << std::endl;
      }

      const double deltadotp = adotp - bdotp;
      const double sigmadotp = adotp + bdotp;
      const double esigmap = eap + ebp;
      const double edeltap = eap - ebp;

      const double Kss = - deltadotp * eahbh;
      const double Kcc =  - esigmap * ahdotbh;
      const double Ks = sigmadotp;
      const double Kc = esigmap;
      const double Kcs = -sigmadotp * ahdotbh - edeltap * eahbh;
      const double K1 = 0;

      m_info.Kss = Kss;
      m_info.Kcc = Kcc;
      m_info.Ks = Ks;
      m_info.Kc = Kc;
      m_info.Kcs = Kcs;
      m_info.K1 = K1;


      if (m_debugMode) {
	std::cout << "{Nss, Ncc, Ns, Nc, Ncs, N1} " << Kss << " " << Kcc << " " << Ks << " " << Kc << " " << Kcs << " " << K1 << std::endl;
      }

      // For a polynomial where Nss is the coeff of sin^2 and Nsc is the coefficent of the sin cos term, etc, then if you substitute sin=sqrt(1-cos^2) and then rearrange for a polynomial exclusively in cos, you get a polynomial with the following coeffs ... see 20100622a-mt2-masslessinvisiblesandmasslesschi-lookingfornicetriganswers-dumpdotprods-1.nb

      /* For Cos coeffs
	 const double coeffCos0 = (K1 - Ks + Kss)*(K1 + Ks + Kss);
	 const double coeffCos1 = 2.*(-Kcs* Ks + Kc*(K1 + Kss));
	 const double coeffCos2 = square(Kc) - square(Kcs) + square(Ks) + 2.*(Kcc - Kss)*(K1 + Kss);
	 const double coeffCos3 = 2.*(Kcs*Ks + Kc*(Kcc - Kss));
	 const double coeffCos4 = square(Kcs) + square(Kcc - Kss);
      */

      // For Sine Coeffs
      const double coeffSin0 = (K1 - Kc + Kcc)*(K1 + Kc + Kcc);
      const double coeffSin1 = 2.*(-Kcs* Kc + Ks*(K1 + Kcc));
      const double coeffSin2 = square(Ks) - square(Kcs) + square(Kc) + 2.*(Kss - Kcc)*(K1 + Kcc);
      const double coeffSin3 = 2.*(Kcs*Kc + Ks*(Kss - Kcc));
      const double coeffSin4 = square(Kcs) + square(Kss - Kcc);

      // NB we could get the same for coeffs of sin by interchanging s with c (and notint that Ksc should be re-written  as Kcs).

      // now find the roots of this quartic by the method of Farrari.
     
      const double AA=coeffSin4;
      const double BB=coeffSin3;
      const double CC=coeffSin2;
      const double DD=coeffSin1;
      const double EE=coeffSin0;
  
      m_info.A = AA;
      m_info.B = BB;
      m_info.C = CC;
      m_info.D = DD;
      m_info.E = EE;


      std::vector<double> goodSinValues;
     
      if (m_rootFindingMethod == FERRARI) {
  
	std::vector<double> sinValues(4);

	// FIND ROOTS OURSELVES
	if (m_debugMode) {
	  std::cout << "{AA,BB,CC,DD,EE} " 
		    << AA/magPtmiss/magPtmiss << " " 
		    << BB/magPtmiss/magPtmiss << " " 
		    << CC/magPtmiss/magPtmiss << " " 
		    << DD/magPtmiss/magPtmiss << " " 
		    << EE/magPtmiss/magPtmiss << " " 
		    << std::endl;
	}

	const double aa = -(3./8.)*square(BB/AA) + CC/AA;
	const double bb = +(1./8.)*cube(BB/AA) - BB*CC/(2.* square(AA)) + DD/AA;
	const double gg = 
	  - (3./256.)*square(square(BB/AA)) 
	  + (1./16.)*CC*square(BB)/cube(AA) 
	  - BB*DD/(4. *square(AA)) + EE/AA;

	if (m_debugMode) {
	  std::cout << "{aa,bb,gg} " << aa << " " << bb << " " << gg << std::endl;
	}

	const double PP = -square(aa)/12. - gg;
	const double QQ = -cube(aa)/108. + aa*gg/3. - square(bb)/8.;
	const double BITINROOT = square(QQ)/4. + cube(PP)/27.;
	// BITINROOT might be negative, and it is OK if it is, because things should be real again by the time we contstruct yy.  So temprarily drift into complex numbers:
	if (m_debugMode) {
	  std::cout << "PP QQ BITINROOT " << PP <<" " << QQ << " " << BITINROOT << std::endl;
	}

	/* THIS IS FOR REAL ONLY
	   const double RR = -QQ/2. + sqrt(BITINROOT);

	   if (m_debugMode) {
	   std::cout << "PP QQ BITINROOT RR " << PP << " " << QQ << " " << BITINROOT << " " << RR << std::endl;
	   }

	   const double UU = pow(RR,1./3.);

	   const double yy = -(5./6.)*aa + UU - PP/(3.*UU);

	*/

	const std::complex<double> complexBITINROOT(BITINROOT,0);
	const std::complex<double> complexRR = sqrt(complexBITINROOT)-QQ/2.;
	const std::complex<double> complexUU = pow(complexRR,1./3.);

	// Note the special line should have a separate case to deal with UU=0.

	std::complex<double> complexyy = + complexUU - PP/(complexUU*3.) -(5./6.)*aa ;

	if (m_debugMode) {
	  std::cout << "complexBITINROOT " << complexBITINROOT << std::endl;
	  std::cout << "complexRR " << complexRR << std::endl;
	  std::cout << "complexUU " << complexUU << std::endl;
	  std::cout << "complexyy " << complexyy << std::endl;
	}

	double yy = complexyy.real();

	double tmp1 = aa + 2.*yy;
	if (tmp1<0) {
	  if (m_debugMode) {
	    std::cout << "WARNING: Debug message from Analytic_Mt2_2220_Calculator ... send this infor to Lester.  About to throw exception." << std::endl;
	    std::cout << "example=({" 
		      << "  a -> " << a
		      << ", b -> " << b
		      << ", p -> " << magPtmiss
		      << ", thetaap -> ArcTan[" << adotp << "," << eap << "]"
		      << ", thetabp -> ArcTan[" << bdotp << "," << ebp << "]"
		      << "})/.correct" << std::endl;
	    std::cout << "WARNING: Oh dear ... a thing that should not be negative is negative and has val " << tmp1 << " !  Attempting to fix." << std::endl;
	    std::cout << "WARNING: A this point PP = " << PP << " complexUU = " << complexUU << " QQ = " << QQ << " aa = " << aa << " complexyy = " << complexyy << " tmp1 = " << tmp1 << std::endl;
	  }
	  // This is most probably a sign that UU was getting very close to zero, and the subtraction UU-PP/(3*UU) became unstable.  In the limit of small UU, PP/(3UU) goes to the cube root of QQ (and QQ should be positive) .. so let's try that again:
	  complexyy = + complexUU - pow(std::complex<double>(QQ,0),1./3.) -(5./6.)*aa ;
	  yy = complexyy.real();
	  tmp1 = aa + 2.*yy;
       
	  // is tmp1 STILL less than zero?
	  if (tmp1<0) {
	    if (m_debugMode) {
	      std::cout << "WARNING: Oh dear, even after attempt to fix things, tmp1 is still negative.  Now complexyy = " << complexyy << " tmp1 = " << tmp1 << std::endl;
	    }
	    // let us hope this is rounding, so that we may proceed without crashing.
	    tmp1=0;
	  } else {
	    std::cout << "WARNING: appears to have been fixed" << std::endl;
	  }

	}

	const double WW = sqrt(tmp1);  

	if (m_debugMode) {
	  std::cout << "yy WW " << " " << yy << " " << WW << std::endl;
	}

	{
	  unsigned int rootcount = 0;
	  unsigned int goodrootcount = 0;
	  for (int sign1 = -1; sign1<=+1; sign1 += 2) {
	    for (int sign2 = -1; sign2<=+1; sign2 += 2) {
	   
	      const double tmp2 = -(3.*aa + 2.* yy + sign1 * 2.*bb/WW);
	   
	      const double root = 
		- BB/(4.* AA)
		+ (sign1 * WW 
		   - 
		   sign2 *sqrt(tmp2)
		   )/2. ;
	   
	      sinValues[rootcount] = root;	 
	      ++rootcount;
	      if (m_debugMode) {
		std::cout << "A ROOT = " << root << std::endl;
	      }

	      if (tmp2<0) {
		if (m_debugMode) {
		  std::cout << "ROOT " << rootcount << " is complex." << std::endl;
		}
	      }else { 
		if (m_debugMode) { 
		  std::cout << "ROOT " << rootcount << " is real." << std::endl;
		}
		goodSinValues.push_back(root);
		++goodrootcount;
	      }
	   
	    }
	  }
	}
     
	// by this stage, goodSinValues should have a set of good sin vales, roots, 
     
        // END OF finding quartic roots by OURSELVES
      } else {
	// get someone else to calculate the quartic roots:

	std::complex<double> r1;
	std::complex<double> r2;
	std::complex<double> r3;
	std::complex<double> r4;
       
       
	Mt2::quartic::r8poly4_root(AA,BB,CC,DD,EE,r1,r2,r3,r4);

	
	m_info.r1 = r1;
	m_info.r2 = r2;
	m_info.r3 = r3;
	m_info.r4 = r4;

	// Now how do we separate the real roots from the "almost" real roots from the complex roots?

	insertIfOK(r1,goodSinValues);
	insertIfOK(r2,goodSinValues);
	insertIfOK(r3,goodSinValues);
	insertIfOK(r4,goodSinValues);

       
      }
       
      std::set<double> answers;
      for (int i=0; i<(int)goodSinValues.size(); ++i) {
	if (m_debugMode) { 
	  std::cout << "GOOD ROOT " << i+1 << " = " << goodSinValues[i] << std::endl;
	}

	// this valuie for sin(theta) could correspond EITHER to 
	//       cos(theta) = +sqrt(1-sin^2(theta))
	// OR to
	//       cos(theta) = -sqrt(1-sin^2(theta))
	// since when we squared our multinomial in {ss, cc, sc, c, s} to get rid of the sqrt() signs, we lost the sign in front of them,
       
	const double s = goodSinValues[i];
	const double modc  = sqrt(1.-s*s); /* modulus of cosine */

	// In the case that the positive square root is needed, then the following two quantities should be identical, otherwise they should differ in sign:
       
	const double LHS = K1 + Kcc + s *(Ks + (-Kcc + Kss)* s);
	const double RHS = -(Kc + Kcs* s) * modc;
	if (m_debugMode) { 
	  std::cout << "LHS RHS " << LHS << " " << RHS << std::endl;
	}

	const bool usePositiveSquareRoot = (fabs(LHS-RHS)<fabs(LHS+RHS));

	if (m_debugMode) { 
	  if (usePositiveSquareRoot) {
	    std::cout << "NEED positive SQUARE ROOT" << std::endl;
	  } else {
	    std::cout << "NEED negative SQUARE ROOT" << std::endl;
	  }
	}
	 
	const double sinTheta = s;
	const double cosTheta = (usePositiveSquareRoot?+modc:-modc);

	if (m_debugMode) { 
	  std::cout << "cosTheta = " << cosTheta << " sinTheta = " << sinTheta << std::endl;
	}

	// Now the balanced condition implies that 2 a_mu p^mu = 2 b_mu q^mu hence


	//const double p_over_q =
	//	 (b*(1 - cosTheta * ahdotbh - sinTheta * eahbh)) /
	//	 (a*(1 - cosTheta * ahdotbh + sinTheta * eahbh)) ;
	const double p_over_q = -
	  (( cosTheta * eap + sinTheta * adotp) * b )/
	  (( cosTheta * ebp + sinTheta * bdotp) * a );

	if (m_debugMode) {  
	  std::cout << "a b = " << a << " " << b << std::endl;
	  std::cout << "p over q = " << p_over_q << std::endl;
	}

	// The following quantity should be zero as it (times q) is the cpt of p+q in the direc perp to ptmiss:
	const double shouldBeZero = p_over_q*(cosTheta * ebp + sinTheta * bdotp)/b + (cosTheta*eap + sinTheta*adotp)/a;   
	if (m_debugMode) { 
	  std::cout << "should be zero : " << shouldBeZero << std::endl;
	}

	// The following quantity should NOT be zero as it (times q) is the square of the cpt of p+q in the direc of ptmiss:
	const double thing = p_over_q*(cosTheta * bdotp - sinTheta * ebp)/b + (cosTheta*adotp - sinTheta*eap)/a;   
	if (m_debugMode) { 
	  std::cout << "should NOT be zero : " << thing << std::endl;
	}

	// thing*q = magPtmiss^2  ... hence ... q = magPtmiss^2 /thing;
       
	const double q= square(magPtmiss) /thing;
	const double p = p_over_q * q;

	if (m_debugMode) { 
	  std::cout << "p and q are " << p << " " << q << std::endl;       
	}
	if (p>=0 && q>=0) {
	  const double MT2SQ1 = 2.*a*p*(1 - cosTheta * ahdotbh - sinTheta * eahbh);
	  const double MT2SQ2 = 2.*b*q*(1 - cosTheta * ahdotbh + sinTheta * eahbh);
       
	  if (m_debugMode) { 
	    std::cout << "MT2SQ 1 and 2 are " << MT2SQ1 << " " << MT2SQ2 << std::endl;
	  }
	  // MT2SQ1 and MT2SQ2 should be equal, but due to rounding errors they may not be ... so lets return the mean:
	  answers.insert( (MT2SQ1+MT2SQ2)*0.5);
	}
      }
     
      if (answers.empty()) {
       
	// If we got here, then there was no real root with p and q >0 ... so something went wrong:
	if (m_debugMode) { 
	  std::cout << "Warning ... no physical root found in Analytic_Mt2_2220_Calculator!" << std::endl; 
	}
	throw no_mt2_value_found; // no answer found
      } else {

	// answers should be sorted, and the first one should be the smallest, so
	return *(answers.begin());

      }

    }} // end of mt2_2220_Sq method 
 


} // end of Mt2 Namespace
