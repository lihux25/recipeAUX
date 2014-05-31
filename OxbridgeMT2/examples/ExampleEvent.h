// Example from the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#ifndef EXAMPLEEVENT_H
#define EXAMPLEEVENT_H

/**
 *   Imagine an LHC SUSY event in which two sleptons (sl_A and sl_B)
 *   each decay into a visible leptons and an invisible
 *   massive neutralino, where the invisiable neutralino
 *   has mass "m_invis".  Suppose the visible lepton from A has mass
 *   m_a_vis and the visible lepton from B has mass m_b_vis; 
 *
 *  I.e.
 *         p+ p-  --->  Initial+final-state-radiation "G"
 *                           +
 *                      Hard-process
 *
 * where
 *          Hard-process -->   sl_A   sl_B       (2 sleptons)
 *
 * and
 *           sl_A  -->  vis_A + invis_A
 * while
 *           sl_B  -->  vis_B + invis_B.
 *
 * We could represent the momenta in such an example as follows:
 *
 */


#include <stdlib.h>
#include "Mt2/Mt2Vectors.h"
#include <vector>

struct ExampleEvent {

  ExampleEvent() {

    // Here is a fairly uninteresting and arbitrary event

    m_rootS = 14000;

    const double vis_A_mass=100;
    const double vis_B_mass=150;

    m_invis_mass=100;

    m_vis_A_4mom.setVectM( 410,  20,-20, vis_A_mass);
    m_vis_B_4mom.setVectM(-210,-300, 44, vis_B_mass);

    m_otherVisible4Mom = Mt2::LorentzVector(Mt2::LorentzVector::InitEPxPyPz(0,0,0,0)); // maybe we should make this more interesting -- after all, most events have some transverse momentum for the hard process as a result of initial state or final state radiation.
    
  }
  

public:
  static double shoot() { return static_cast<double>(rand())/(static_cast<double>(RAND_MAX)+1); }
  static double flipshoot() { return (shoot()-0.5)*2.0; }
public:
  static std::vector<double> randomUnitVector(const unsigned int dim) {
    std::vector<double> direc(dim);
    while (true) {
      double magsq=0;
      for (unsigned int i=0; i<dim; ++i) {
	direc[i] = flipshoot();
	magsq += direc[i]*direc[i];
      }
      if (magsq<1 && magsq>0) {
	const double mag=sqrt(magsq);
	for (unsigned int i=0; i<dim; ++i) {
	  direc[i] /= mag;
	}
	return direc;
      }
    }
  }
  static std::vector<double> randomBoostVector(const unsigned int dim) {
    std::vector<double> direc = randomUnitVector(dim);
    double mag=10;
    while (mag>=1 || mag<0) {
      mag = shoot();
    }
    for (unsigned int i=0; i<dim; ++i) {
      direc[i] *= mag;
    }
    return direc;
  }
public:

  void resetRandom330Massless() {
	m_rootS = 14000;
	const double vis_A_mass = 0;
	const double vis_B_mass = 0;
	m_invis_mass = 300*shoot();
	m_vis_A_4mom.setVectM(200*flipshoot(), 200*flipshoot(), 200*flipshoot(), vis_A_mass);
	m_vis_B_4mom.setVectM(200*flipshoot(), 200*flipshoot(), 200*flipshoot(), vis_B_mass);
	m_otherVisible4Mom.setVectM(0,0,0,0);
  }

  void resetRandom330() {
	m_rootS = 14000;
	const double vis_A_mass = 300*shoot();
	const double vis_B_mass = 300*shoot();
	m_invis_mass = 300*shoot();
	m_vis_A_4mom.setVectM(200*flipshoot(), 200*flipshoot(), 200*flipshoot(), vis_A_mass);
	m_vis_B_4mom.setVectM(200*flipshoot(), 200*flipshoot(), 200*flipshoot(), vis_B_mass);
	m_otherVisible4Mom.setVectM(0,0,0,0);
  }

  void resetZhenA_332() {

    m_rootS = 14000;

    const double vis_A_mass = 4.672030912e-06;
    const double vis_B_mass = 9.5367431641e-07;

    m_invis_mass = 0.1;
    m_vis_A_4mom      .setVectM(-58.117035689, 195.26159258, 0, vis_A_mass);
    m_vis_B_4mom      .setVectM(+12.434843065, -64.46291702, 0, vis_B_mass);
    
    const double pmissx = 111.27094511;
    const double pmissy = -70.456482466;

    m_otherVisible4Mom.setVectM(
	- pmissx - m_vis_A_4mom.px -  m_vis_B_4mom.px,
	- pmissy - m_vis_A_4mom.py -  m_vis_B_4mom.py,
	0,
	0 );

  }

  void resetCowden() {

    double ae,ax,ay,az;
    double be,bx,by,bz;
    
    /*ae = 652665.76752;
    ax = -473414.59051;
    ay=-445084.48633;
    az=551464.72166;
    be=1183824.976;
    bx=217398.70017;
    by=-142189.5412;
    bz=1154597.1998;*/

ae=728651.04085;
ax=-473414.59051;
ay=-445084.48633;
az=-323974.8991;
be=296301.75491;
bx=27658.018322;
by=-27198.208313;
bz=-293645.52692;


    m_rootS = 14000000; // thisi is a cunning "mev" test

    const double m_a_vis_sq = ae*ae-ax*ax-ay*ay-az*az;
    const double m_b_vis_sq = be*be-bx*bx-by*by-bz*bz;

    if (m_a_vis_sq<0) {
      std::cout << "Warning cowden masq negative! " << m_a_vis_sq << std::endl;
    }
    if (m_b_vis_sq<0) {
      std::cout << "Warning cowden mbsq negative! " << m_b_vis_sq << std::endl;
    }

    const double vis_A_mass = sqrt(fabs(m_a_vis_sq));
    const double vis_B_mass = sqrt(fabs(m_b_vis_sq));

    //m_invis_mass = 117*1000;
m_invis_mass=0;


    m_vis_A_4mom      .setVectM(ax,ay,az, vis_A_mass);
    m_vis_B_4mom      .setVectM(bx,by,bz, vis_B_mass);
    
    const double pmissx=445953.9375; //= -437244.65625;
    const double pmissy=466779.15625; //=  241964.375;


    m_otherVisible4Mom.setVectM(
	- pmissx - m_vis_A_4mom.px -  m_vis_B_4mom.px,
	- pmissy - m_vis_A_4mom.py -  m_vis_B_4mom.py,
	0,
	0 );

  }

  void resetRandom222() {
    resetRandom332(0);
  }
  void resetRandom332(const double visMassMax=300) {
	m_rootS = 14000;
	const double vis_A_mass = visMassMax*shoot();
	const double vis_B_mass = visMassMax*shoot();
	const double other_mass = 300*shoot();
	m_invis_mass = 300*shoot();
	m_vis_A_4mom      .setVectM(300*flipshoot(), 300*flipshoot(), 300*flipshoot(), vis_A_mass);
	m_vis_B_4mom      .setVectM(300*flipshoot(), 300*flipshoot(), 300*flipshoot(), vis_B_mass);
	m_otherVisible4Mom.setVectM(300*flipshoot(), 300*flipshoot(), 300*flipshoot(), other_mass);
  }

  // masses and energies
  double vis_A_mass() const { return m_vis_A_4mom.m(); }
  double vis_B_mass() const { return m_vis_B_4mom.m(); }
  double invis_mass() const { return m_invis_mass; }
  double      rootS() const { return m_rootS; }

  // transverse 2-vector quantities:
  Mt2::TwoVector pT_Vis_A    () const { return transverse(m_vis_A_4mom      ); }
  Mt2::TwoVector pT_Vis_B    () const { return transverse(m_vis_B_4mom      ); }
  // Upstream Transverse Momentum two-vector component: 
  Mt2::TwoVector pT_Vis_Other() const { return transverse(m_otherVisible4Mom); }
  Mt2::TwoVector pT_Vis      () const { return transverse(p_Vis()); }
  Mt2::TwoVector pT_Miss     () const { return -pT_Vis(); }
  
  // four-vector quantities:
  Mt2::LorentzVector p_Vis_A() const { return m_vis_A_4mom; }
  Mt2::LorentzVector p_Vis_B() const { return m_vis_B_4mom; }
  Mt2::LorentzVector p_Vis_Other() const { return m_otherVisible4Mom; }
  Mt2::LorentzVector p_Vis() const {
    return m_vis_A_4mom + m_vis_B_4mom + m_otherVisible4Mom;
  }
  
  void resetSpecial222() {
	m_rootS = 14000;
	const double vis_A_mass = 0;
	const double vis_B_mass = 0;
	const double other_mass = 300*shoot();
	m_invis_mass = 100;
	m_vis_A_4mom      .setVectM(1.2*cos(1.0), 1.2*sin(1.0),0, vis_A_mass);
	m_vis_B_4mom      .setVectM(0.8*cos(2.1), 0.8*sin(2.1),0, vis_B_mass);
	m_otherVisible4Mom.setVectM(-1.2*cos(1.0)-0.8*cos(2.1)-1,-1.2*sin(1.0)-0.8*sin(2.1)-0,0, other_mass);

  }

  // Lorentz-transverse quantities
  //  (i.e. quantities neeting to know the 3 components:
  //       (1) px
  //       (2) py
  //       (3) et ==def== sqrt(px^2+py^2+m^2)
  //  )

  Mt2::LorentzTransverseVector ltv_Vis_A() const { return ltv(m_vis_A_4mom); };
  Mt2::LorentzTransverseVector ltv_Vis_B() const { return ltv(m_vis_B_4mom); };
  
private:

  double m_rootS;

  double m_invis_mass;

  Mt2::LorentzVector m_vis_A_4mom;
  Mt2::LorentzVector m_vis_B_4mom;
  Mt2::LorentzVector m_otherVisible4Mom;
  
private:
  static Mt2::TwoVector transverse(const Mt2::LorentzVector & v) {
    return Mt2::TwoVector(v.px, v.py);
  }
  static Mt2::LorentzTransverseVector ltv(const Mt2::LorentzVector & v) {
    return v.getLorentzTransverseVector();
  }
};

#endif
