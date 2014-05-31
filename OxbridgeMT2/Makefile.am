## OxbridgeKinetics:

## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS = subdir-objects

## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

## Define a libtool archive target "liboxbridgekinetics-@EXAMPLE_API_VERSION@.la", with
## @EXAMPLE_API_VERSION@ substituted into the generated Makefile at configure
## time.
## The libtool archive file (.la) will be installed into the directory named
## by the predefined variable $(bindir), along with the actual shared library
## file (.so).
lib_LTLIBRARIES = liboxbridgekinetics-@EXAMPLE_API_VERSION@.la

## Define the source file list for the "liboxbridgekinetics-@EXAMPLE_API_VERSION@.la"
## target.  Note that @EXAMPLE_API_VERSION@ is not interpreted by Automake and
## will therefore be treated as if it were literally part of the target name,
## and the variable name derived from that.
## The file extension .cc is recognized by Automake, and makes it produce
## rules which invoke the C++ compiler to produce a libtool object file (.lo)
## from each source file.  Note that it is not necessary to list header files
## which are already listed elsewhere in a _HEADERS variable assignment.
liboxbridgekinetics_@EXAMPLE_API_VERSION@_la_SOURCES = \
src/Analytic_Mt2_330_Calculator.cpp \
src/Mt2Vectors.cpp \
src/Mt2Calculator.cpp \
src/Quartic.h \
src/mt2_bisect.h \
src/mt2_bisect.cpp \
src/ChengHanBisect_Mt2_332_Calculator.cpp \
src/Analytic_Mt2_2220_Calculator.cpp \
src/Mt2LorentzTransverseVector.cpp \
src/Mt2LorentzVector.cpp \
src/Mt2Util.cpp \
src/Mt2TwoVector.cpp \
src/M2C_332Calculator.cpp \
src/MC_330_Calculator.cpp \
src/MCT_330_Calculator.cpp \
src/MCTll_332_Calculator.cpp \
src/MCTT_332_Calculator.cpp


if GOT_MINUIT2
liboxbridgekinetics_@EXAMPLE_API_VERSION@_la_SOURCES += \
src/Advanced_Mt2_332_Calculator.cpp \
src/Basic_Mt2_332_Calculator.cpp \
src/Basic_Mt2_4441_Calculator.cpp \
src/Basic_Mt2_Asymmetric332_Calculator.cpp \
src/Basic_Nt2_332_Calculator.cpp \
src/Basic_MCT2_332_Calculator.cpp \
src/mT2Fcn_4441_a.cpp \
src/SUSYPhys_Mt2_222_Calculator.cpp
endif

## Instruct libtool to include ABI version information in the generated shared
## library file (.so).  The library ABI version is defined in configure.ac, so
## that all version information is kept in one place.
liboxbridgekinetics_@EXAMPLE_API_VERSION@_la_LDFLAGS = -version-info $(EXAMPLE_SO_VERSION)

## Now add ROOT support:
if GOT_MINUIT2
liboxbridgekinetics_@EXAMPLE_API_VERSION@_la_LDFLAGS += -R @ROOTLIBDIR@
liboxbridgekinetics_@EXAMPLE_API_VERSION@_la_LIBADD  = -lCint -lCore -lMinuit2 @ROOTAUXLIBS@  -L@ROOTLIBDIR@
AM_CPPFLAGS = -I@ROOTINCDIR@
##AM_LDFLAGS  = -L@ROOTLIBDIR@
endif

## Sometiems the ROOT libraries are so horribly compiled, that TMinuit needs all of these! Gah!!
if GOT_MINUIT2
liboxbridgekinetics_@EXAMPLE_API_VERSION@_la_LIBADD  += -lGraf -lRIO -lHist -lMatrix -lCore -lCint -lMathCore
endif

LDADD=liboxbridgekinetics-@EXAMPLE_API_VERSION@.la

## Define the list of public header files and their install location.  The
## nobase_ prefix instructs Automake to not strip the directory part from each
## filename, in order to avoid the need to define separate file lists for each
## installation directory.  This only works if the directory hierarchy in the
## source tree matches the hierarchy at the install location, however.
oxbridgekinetics_includedir = $(includedir)/oxbridgekinetics-$(EXAMPLE_API_VERSION)
nobase_oxbridgekinetics_include_HEADERS =  \
Mt2/Advanced_Mt2_332_Calculator.h \
Mt2/AlphaT_Multijet_Calculator.h \
Mt2/Analytic_Mt2_2220_Calculator.h \
Mt2/Analytic_Mt2_330_Calculator.h \
Mt2/Analytic_Mt2_332_Assistor.h \
Mt2/Basic_M2C_332_Calculator.h \
Mt2/Basic_MPairProd_Calculator.h \
Mt2/Basic_MtGen_330_Calculator.h \
Mt2/Basic_MtGen_332_Calculator.h \
Mt2/ChengHanBisect_Mt2_332_Calculator.h \
Mt2/Frugal_MPairProd_Calculator.h \
Mt2/Frugal_MtGen_330_Calculator.h \
Mt2/Frugal_MtGen_332_Calculator.h \
Mt2/MC_330_Calculator.h \
Mt2/MCT2_332_Calculator.h \
Mt2/MCT_330_Calculator.h \
Mt2/MCTll_332_Calculator.h \
Mt2/MCTT_332_Calculator.h \
Mt2/MPairProd_Calculator.h \
Mt2/Mt2_300_Calculator.h \
Mt2/Mt2_302_Calculator.h \
Mt2/Mt2_330_Calculator.h \
Mt2/Mt2_332_Calculator.h \
Mt2/Mt2ApproximatingAdapter_332_from_330.h \
Mt2/Mt2_Asymmetric332_Calculator.h \
Mt2/Mt2Calculator.h \
Mt2/Mt2Calculators.h \
Mt2/Mt2LorentzTransverseVector.h \
Mt2/Mt2LorentzVector.h \
Mt2/Mt2MinimiserGsl.h \
Mt2/Mt2Minimiser.h \
Mt2/Mt2MinimiserNag.h \
Mt2/Mt2TwoVector.h \
Mt2/Mt2Units.h \
Mt2/Mt2Util.h \
Mt2/Mt2Vectors.h \
Mt2/MtGen_330_Calculator.h \
Mt2/MtGen_332_Calculator.h \
Mt2/Nt2_332_Calculator.h \
Mt2/SolutionType.h

if GOT_MINUIT2
nobase_oxbridgekinetics_include_HEADERS += \
Mt2/Basic_MCT2_332_Calculator.h \
Mt2/Basic_Mt2_332_Calculator.h \
Mt2/Basic_Mt2_Asymmetric332_Calculator.h \
Mt2/Basic_Nt2_332_Calculator.h \
Mt2/SUSYPhys_Mt2_222_Calculator.h \
Mt2/mT2Fcn_4441_a.h
endif

## The generated configuration header is installed in its own subdirectory of
## $(libdir).  The reason for this is that the configuration information put
## into this header file describes the target platform the installed library
## has been built for.  Thus the file must not be installed into a location
## intended for architecture-independent files, as defined by the Filesystem
## Hierarchy Standard (FHS).
## The nodist_ prefix instructs Automake to not generate rules for including
## the listed files in the distribution on 'make dist'.  Files that are listed
## in _HEADERS variables are normally included in the distribution, but the
## configuration header file is generated at configure time and should not be
## shipped with the source tarball.
oxbridgekinetics_libincludedir = $(libdir)/oxbridgekinetics-$(EXAMPLE_API_VERSION)/include
nodist_oxbridgekinetics_libinclude_HEADERS = oxbridgekineticsconfig.h

## Install the generated pkg-config file (.pc) into the expected location for
## architecture-dependent package configuration information.  Occasionally,
## pkg-config files are also used for architecture-independent data packages,
## in which case the correct install location would be $(datadir)/pkgconfig.
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = oxbridgekinetics-$(EXAMPLE_API_VERSION).pc

## Define an independent executable script for inclusion in the distribution
## archive.  However, it will not be installed on an end user's system due to
## the noinst_ prefix.
dist_noinst_SCRIPTS = autogen.sh

noinst_PROGRAMS =\
example_helloWorld_Mt2_Without_Minuit2_Example \
example_Analytic_Mt2_330_Calculator_Example \
example_ChengHanBisect_Mt2_332_Calculator_Example \
example_MC_330_Calculator_Example \
example_AlphaT_Multijet_Calculator_Example \
example_Basic_MtGen_330_Calculator_Example \
example_Frugal_MtGen_332_Calculator_Example \
example_Basic_MtGen_332_Calculator_Example  \
example_Frugal_MtGen_330_Calculator_Example \
example_Frugal_MPairProd_Calculator_Example \
example_M2C_332_Calculator_Example \
example_MCT_330_Calculator_Example

if GOT_MINUIT2
noinst_PROGRAMS += \
example_helloWorld_Mt2_Example \
example_Analytic_Mt2_2220_Calculator_Example \
example_Advanced_Mt2_332_Calculator_Example \
example_Analytic_Mt2_330_Stressor_Example \
example_AlgorithmComparison_Example \
example_Approximating_Example \
example_Basic_MCT2_332_Calculator_Example \
example_Basic_MPairProd_Calculator_Example \
example_Basic_Mt2_332_Calculator_Example \
example_Basic_Mt2_4441_Calculator_Example \
example_Basic_Mt2_Asymmetric332_Calculator_Example \
example_Basic_Nt2_332_Calculator_Example \
example_ComparisonTester_Example \
example_HowBadIs_No_ISR_Approx_Example \
example_PtmissQFactorRecoil_Test_Example \
example_M2C_332_Calculator_Using_Minuit2_Example \
example_SUSYPhys_Mt2_222_Calculator_Example
endif


#example_ChengHanBisect_Mt2_332_Calculator_Example_LDADD=liboxbridgekinetics-@EXAMPLE_API_VERSION@.la





example_helloWorld_Mt2_Without_Minuit2_Example_SOURCES=examples/helloWorld_Mt2_Without_Minuit2_Example.cpp
example_Analytic_Mt2_330_Calculator_Example_SOURCES=examples/Analytic_Mt2_330_Calculator_Example.cpp examples/ExampleEvent.h
example_ChengHanBisect_Mt2_332_Calculator_Example_SOURCES=examples/ChengHanBisect_Mt2_332_Calculator_Example.cpp
example_MC_330_Calculator_Example_SOURCES=examples/MC_330_Calculator_Example.cpp
example_MCT_330_Calculator_Example_SOURCES=examples/MCT_330_Calculator_Example.cpp
example_AlphaT_Multijet_Calculator_Example_SOURCES=examples/AlphaT_Multijet_Calculator_Example.cpp
example_Frugal_MtGen_330_Calculator_Example_SOURCES=examples/Frugal_MtGen_330_Calculator_Example.cpp
example_Basic_MtGen_330_Calculator_Example_SOURCES=examples/Basic_MtGen_330_Calculator_Example.cpp
example_Frugal_MPairProd_Calculator_Example_SOURCES=examples/Frugal_MPairProd_Calculator_Example.cpp
example_Basic_MtGen_332_Calculator_Example_SOURCES=examples/Basic_MtGen_332_Calculator_Example.cpp
example_Frugal_MtGen_332_Calculator_Example_SOURCES=examples/Frugal_MtGen_332_Calculator_Example.cpp
example_M2C_332_Calculator_Example_SOURCES=examples/M2C_332_Calculator_Example.cpp

if GOT_MINUIT2
example_helloWorld_Mt2_Example_SOURCES=examples/helloWorld_Mt2_Example.cpp
example_Analytic_Mt2_330_Stressor_Example_SOURCES=examples/Analytic_Mt2_330_Stressor.cpp
example_Advanced_Mt2_332_Calculator_Example_SOURCES=examples/Advanced_Mt2_332_Calculator_Example.cpp
example_AlgorithmComparison_Example_SOURCES=examples/AlgorithmComparison_Example.cpp
example_Analytic_Mt2_2220_Calculator_Example_SOURCES=examples/Analytic_Mt2_2220_Calculator_Example.cpp
example_Approximating_Example_SOURCES=examples/Approximating_Example.cpp
example_Basic_MCT2_332_Calculator_Example_SOURCES=examples/Basic_MCT2_332_Calculator_Example.cpp
example_Basic_MPairProd_Calculator_Example_SOURCES=examples/Basic_MPairProd_Calculator_Example.cpp
example_Basic_Mt2_332_Calculator_Example_SOURCES=examples/Basic_Mt2_332_Calculator_Example.cpp
example_Basic_Mt2_4441_Calculator_Example_SOURCES=examples/Basic_Mt2_4441_Calculator_Example.cpp
example_Basic_Mt2_Asymmetric332_Calculator_Example_SOURCES=examples/Basic_Mt2_Asymmetric332_Calculator_Example.cpp
example_Basic_Nt2_332_Calculator_Example_SOURCES=examples/Basic_Nt2_332_Calculator_Example.cpp
example_ComparisonTester_Example_SOURCES=examples/ComparisonTester.cpp
example_HowBadIs_No_ISR_Approx_Example_SOURCES=examples/HowBadIs_No_ISR_Approx.cpp
example_PtmissQFactorRecoil_Test_Example_SOURCES=examples/PtmissQFactorRecoil_Test.cpp
example_M2C_332_Calculator_Using_Minuit2_Example_SOURCES=examples/M2C_332_Calculator_Using_Minuit2_Example.cpp
example_SUSYPhys_Mt2_222_Calculator_Example_SOURCES=examples/SUSYPhys_Mt2_222_Calculator_Example.cpp
endif

dist_data_DATA = examples/demo_1_for_ROOT_usage.C
