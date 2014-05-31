
/*

   This is an example ROOT macro that demonstrates
   how to evalutate MT2 in ROOT semi-interactively.
   It is important to realise that you MUST use a 
   proper c++ compiler (eg g++, or to some extent
   ACLIC) if you want to use ROOT with this library,
   as CINT is not good enough on its own.
   Follow the instructions
   below which have been tested on root v 5.15/08.

   STEP 1:
        run follow the instructions in the INSTALL
	file to build and install the oxbridgekinematics
	library.
       
	Note that if you did:

	    ./configure
	    make
	    make install

	your header files and libraries will have been installed in a
	tree like this: (version numbers may differ)

	    /usr/local/include/oxbridgekinetics-1.0/Mt2/Basic_M2C_332_Calculator.h
	    /usr/local/include/oxbridgekinetics-1.0/Mt2/Mt2Vectors.h
	    /usr/local/include/oxbridgekinetics-1.0/Mt2/Mt2LorentzTransverseVector.h
	    /usr/local/lib/liboxbridgekinetics-1.0.so.1
	    /usr/local/lib/liboxbridgekinetics-1.0.so
	    
       etc, whereas if you did

	    ./configure --prefix=/somewhere/else
	    make
	    make install

       they will have gone to

	    /somewhere/else/include/oxbridgekinetics-1.0/Mt2/Basic_M2C_332_Calculator.h
	    /somewhere/else/include/oxbridgekinetics-1.0/Mt2/Mt2Vectors.h
	    /somewhere/else/include/oxbridgekinetics-1.0/Mt2/Mt2LorentzTransverseVector.h
	    /somewhere/else/lib/liboxbridgekinetics-1.0.so.1
	    /somewhere/else/lib/liboxbridgekinetics-1.0.so
	    

       If you didn't do the "make install" stage at all (but you
       did do the "make") then your library is still sitting in

            ./.libs/liboxbridgekinetics-1.0.so

       
        The instructions in "STEP 2" below assume you installed the
	library to /somewhere/else, but you should alter them to have
	the appropriate path if your did one of the other options above.

   STEP 2:
	Copy me (demo_1_for_ROOT_usage.C) to the directory from which you intend to 
	run ROOT, and then cd into that directory.

   STEP 3: 

	Start ROOT ... eg

		root -l

	and wait for the root prompt to appear.

   STEP 4:
	Run the demo interactively by issuing the following commands
	at the ROOT prompt:

        // Get Mt2 library, and its pre-requisite minuit2 library:
        gSystem->Load("libMinuit2.so");
        gSystem->Load("/somewhere/else/lib/liboxbridgekinetics-1.0.so");

	// tell the interpreter where the include files are for the MT2 library:
	gInterpreter->AddIncludePath("/somewhere/else/include/oxbridgekinetics-1.0");

	// read in the example code and build a library from it
        .L demo_1_for_ROOT_usage.C+
	
	// run the example code:
        demo1();
  
*/

#include "Mt2/Basic_Mt2_332_Calculator.h"

void demo1() {

  Mt2::Basic_Mt2_332_Calculator calc;

  for (unsigned int event=0; event<10; ++event) {
    double m_vis_A_mass=100;
    double m_vis_B_mass=150;
    double m_invis_mass=100;
    double pxMiss=-100;
    dobule pyMiss=+23;   
 
    Mt2::LorentzTransverseVector vis_A(Mt2::TwoVector( 410,   20), m_vis_A_mass);
    Mt2::LorentzTransverseVector vis_B(Mt2::TwoVector(-210, -300), m_vis_B_mass);
    Mt2::TwoVector pT_Miss(pxMiss, pyMiss);

    const double mT2 = calc.mt2_332(vis_A,
				    vis_B,
				    pT_Miss,
				    m_invis_mass);

    std::cout << mT2 << std::endl;
  }

}
