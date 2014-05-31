// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#include "recipeAUX/OxbridgeMT2/interface/Basic_M2C_332_Calculator.h"
#include "recipeAUX/OxbridgeMT2/interface/ChengHanBisect_Mt2_332_Calculator.h"

//#include <iostream>

// M2C part of the Code written by Mario A. Serna Jr
// Last updated on 02 May 2008
// Previous updated on 08 Feb 2008
// The variable M2C was introduced in
//   G.Ross and M.Serna, "Mass Determination of New States at Hadron Colliders", 
//   arXiv:0712.0943 [hep-ph] 06 Dec 2007


// M2C_332Calculator
//   Applies to new state decaying to light-stable state via a three body decay.
//   Needs the edge of the m12 distribution in the three-body decay to get a measure on m2-m1.
//   Now for each event, given the 2+1 momenta on branch A & B and the missing pT 
//   and the mass difference, this will find the smallest mass of the parent (second to lightest new state)
//   consistant with these states 
// Assumes momenta are all given in unist of GeV.

double M2C_332Calculator(Mt2::LorentzTransverseVector& visA, // 2+1 momenta from branch A  (eo, px, py)
					     Mt2::LorentzTransverseVector& visB, // 2+1 momenta from branch B  (eo, px, py)
					     Mt2::TwoVector& ptmiss, // the 2 momenta , the missing transverse momenta (px, py)
					     double deltam,  // The value of m2-m1 read off from the three-body decay m12 end point 
                                             Mt2::Mt2_332_Calculator * mt2CalculatorToUse
                    )
{

  // The calculation consists of finding the zero of the function
  //   g(chi)= chi+ deltam - mt2(chi)

 
  double      invis_mass_new, invis_mass, invis_mass_trialA, invis_mass_trialB;  // variables to be used in search
  double      invis_mass_lo, invis_mass_hi; // keeps track of the bracketing of the zero
//  double      g, gtemp_new, gtemp,gLo,gHi,gtempA,gtempB,dgdchi; // g is the function defined above
  double      gtemp_new, gtemp,gLo,gHi,dgdchi; // g is the function defined above
  double      m2c_result;
//  double      ghighest,gstart;
  double      gstart;
  double      Sqrts=1000;
  
//  double      invis_mass_highest=Sqrts;
//  double      invis_mass_start=0;

  
  double      precision_invis_mass=0.0001;  // the precision at which we will stop the search for m2c
  double      precision_g=0.00001; // the precision at which we will be happy to say we found the zero of g
  int         iteration_count=0;
  int         max_iteration=10000;
  

   
  // Now that we have some visiable momenta and some missing transverse
  // momentum we can calculate MT2.
  
  // First we create the object that is going to do the calculation
  // of MT2 for us.

  Mt2::ChengHanBisect_Mt2_332_Calculator defaultMt2Calculator;


  Mt2::Mt2_332_Calculator & mt2Calculator = ((mt2CalculatorToUse)?(*mt2CalculatorToUse):defaultMt2Calculator);

  int bracketed=0; //this variable will kep track of if Lo and Hi have the solution bracketed.  
  
  // Now we can actually calculate M2C:

 //  std::cout << "Going to calculate M2C with\n"
//	    << "   ltv_Vis_A  = " << visA  << "\n"
//	    << "   ltv_Vis_B  = " << visB  << "\n"
//	    << "   pT_Miss    = " << ptmiss    << "\n"
//      << "   Delta M    = " << deltam     << "\n"
//	    << "   invis_mass = " << invis_mass << std::endl;
  
  // these two check for bounds.  One can only find a solution if the function crosses zero.

  invis_mass_lo=0;    // this is the lower bound on the invisible mass
  invis_mass_hi=deltam/2; // at the LHC, this is the highest concievable mass of a particle that could be produced
                      // therefore, we start with this as our upper bound.
  gLo = deltam + invis_mass_lo  - mt2Calculator.mt2_332( visA, visB, ptmiss, invis_mass_lo);
  gHi = deltam + invis_mass_hi  - mt2Calculator.mt2_332( visA, visB, ptmiss, invis_mass_hi);

  //std::cout << "gLo = " << gLo << std::endl;
  //std::cout << "gHi = " << gHi << std::endl;

  gstart=gLo;

  // We are only interested in the first crossing.
  
  //I'll start by checking the two regions of deltam manually to see if it is bracketd
  //   Bracketd points will be give by gLo and gHi
  // We could check the derivative at the low point and see if one can forcast a crossing.
  
//  invis_mass_mid=deltam; // this is a first mid point where I will check for a crossing.  
//  gMid= deltam + invis_mass_mid - mt2Calculator.mt2_332( visA, visB, ptmiss, invis_mass_mid);  
  
  if (gLo*gHi < 0)
    {
    bracketed=1;
    }
  else
   {
   iteration_count=0;
   while ((!bracketed)&&(invis_mass_hi<Sqrts)&&(iteration_count<max_iteration))
      {
      iteration_count++;
      gLo=gHi;
      invis_mass_lo=invis_mass_hi;

      invis_mass_hi=invis_mass_lo+deltam/2.0;
      gHi = deltam + invis_mass_hi  - mt2Calculator.mt2_332( visA, visB, ptmiss, invis_mass_hi);
      
      if (gLo*gHi < 0)
        {
        bracketed=1;
        }      
      std::cout << "invis mass hi " <<  invis_mass_hi << std::endl;
      } // end while
      
   } // end of initial region does not bracket it

    
  if (bracketed)
    {

    // If there is a zero, we have it bracketed, and now we can do a newton method search for the zero.
    
    // we begin with a random trial point.  All values in units of GeV.
    invis_mass = 10; 
    // calculated g at a trial point inbetween the two extremes.
    gtemp = deltam + invis_mass - mt2Calculator.mt2_332( visA, visB, ptmiss, invis_mass); 

    // if trial point is lower than 
    if (gtemp<0) 
      {
      invis_mass_lo=invis_mass;
      gLo=gtemp;
      }
    else
      {
      invis_mass_hi=invis_mass;
      gHi=gtemp;
      }
    
    while ((fabs(gtemp) > precision_g) &&   (iteration_count<max_iteration))
      {
      // Just a backup to ensure that the code is not stuck in an infinite loop.
      iteration_count++;
      
    // calculate derivative of g at this trial point
    // recall g(chi) = chi + deltam - mt2(chi)
    invis_mass_trialA = invis_mass - precision_invis_mass; 
    invis_mass_trialB = invis_mass + precision_invis_mass;  
    dgdchi = 1 - 
       ( (mt2Calculator.mt2_332( visA, visB, ptmiss, invis_mass_trialA) - 
          mt2Calculator.mt2_332( visA, visB, ptmiss, invis_mass_trialB))/ 
          (invis_mass_trialA-invis_mass_trialB) );
   
    // use the derivative to find the value at the origin (newton's method)
    invis_mass_new = invis_mass - gtemp/dgdchi ;
    gtemp_new = deltam + invis_mass_new - mt2Calculator.mt2_332( visA, visB, ptmiss, invis_mass_new); 
   
    // check to ensure invis_mass_new is within the brackets
    if ((invis_mass_new>=invis_mass_lo)&&(invis_mass_new<=invis_mass_hi))
      {
      invis_mass=invis_mass_new;
      gtemp=gtemp_new;
      }
    else // if newton method did not return an invisible mass within the brackets
        // then one picks one via bisection method
      {
      invis_mass=0.5*(invis_mass_lo+invis_mass_hi);
      gtemp = deltam + invis_mass - mt2Calculator.mt2_332( visA, visB, ptmiss, invis_mass);       
      }

    // now to place this new gtemp as a new, tighter bracket
      if (gtemp<0) 
        {
        invis_mass_lo=invis_mass;
        gLo=gtemp;
        }
      else
        {
        invis_mass_hi=invis_mass;
        gHi=gtemp;
        }
    
      } // end while for newton method search

  /* // If we reach the maximum number of iteration I return a zero as an error signal.
    if (iteration_count==max_iteration)
      invis_mass=-deltam;
    */

    // If we reach the maximum number of iterations, I print a warning but carry on:
    std::cerr << "Warning, maximum number of iterations reached in " << __FILE__ << " at line " << __LINE__ << " so accuracy of MCT calculation is questionable.  Inputs were [visA=" << visA << ", visB=" << visB<< ", ptmiss="<< ptmiss << ", deltaM=" <<deltam<<"]" << std::endl; 

    // this sets m2c_result to the final answer
    m2c_result=deltam+invis_mass;
    
    } // end of  if bracketed actually bracket the zero.
    else
    {
    // error is given by 0 as a result.
    if (  gstart > 0 )
      m2c_result=deltam;
    else
      m2c_result=0;
    }

    // returns the final answer
	return m2c_result;
}


