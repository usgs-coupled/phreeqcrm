README.TXT

        PhreeqcRM--A reaction module for transport simulators 
               based on the geochemical model PHREEQC

This file describes how to compile the PhreeqcRM library and use
it in C++, C, and Fortran programs. The library may be compiled
after generating either a Visual Studio solution for Windows with 
CMake or a Makefile for Linux with configure. All source and 
compilation-related files are available at the following web site:

http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/index.html.

-----------------------------------------------------------------------------------
TABLE OF CONTENTS
-----------------------------------------------------------------------------------

	A. Distribution file
	B. Documentation
	C. Downloading and unzipping the distribution file
	D. Linux compiling, testing, and installing the PhreeqcRM library
	E. Windows compiling, testing, and installing the PhreeqcRM library
	F. Using the library
	G. Contacts

-----------------------------------------------------------------------------------
A. Distribution file
-----------------------------------------------------------------------------------

	The following distribution package (containing the software, 
	compilation-related files, and an advection example in
        C++, C, and Fortran) is available for Linux and Windows:

        phreeqcrm-@VERSION@-@REVISION@.tar.gz      


-----------------------------------------------------------------------------------
B. Documentation
-----------------------------------------------------------------------------------

	Descriptions and examples of PhreeqcRM methods for C++, C, and Fortran
	are available in HTML format in the doc directory of the distribution. 

	Parkhurst, D.L. and Wissmeier, Laurin, 2015, PhreeqcRM--A reaction module 
	   for transport simulators based on the geochemical model PHREEQC: Advances in
	   Water Resources, in press.

	Charlton, S.R., and Parkhurst, D.L., 2011, Modules based on the geochemical model 
	   PHREEQC for use in scripting and programming languages: Computers & Geosciences,
	   vol. 37, no. 10, p. 1653-1663.

	   *The following two user's guides are available in electronic format. Portable
	   Document Format (PDF) files are included in the doc subdirectory of all
	   PHREEQC program distributions.

	Parkhurst, D.L., and Appelo, C.A.J., 2013, Description of input and examples
	   for PHREEQC version 3--A computer program for speciation, batch-reaction,
	   one- dimensional transport, and inverse geochemical calculations: U.S.
	   Geological Survey Techniques and Methods, book 6, chap. A43, 497 p.
	   http://pubs.usgs.gov/tm/06/a43/.

	Parkhurst, D.L., and Appelo, C.A.J., 1999, User's guide to PHREEQC (Version
	   2)-- A computer program for speciation, batch-reaction, one-dimensional
	   transport, and inverse geochemical calculations: U.S. Geological Survey
	   Water-Resources Investigations Report 99-4259, 312 p.
	   http://pubs.er.usgs.gov/publication/wri994259.

	   *The following two reports document the theory and implementation of
	   isotopes in PHREEQC. Portable Document Format (PDF) of Thorstenson and
	   Parkhurst (2002) is included in the doc subdirectory of the PHREEQC program
	   distribution.

	Thorstenson, D.C., and Parkhurst, D.L., 2002, Calculation of individual isotope
	   equilibrium constants for implementation in geochemical models: U.S.
	   Geological Survey Water-Resources Investigations Report 02-4172, 129 p.
	   http://pubs.er.usgs.gov/publication/wri024172.

	Thorstenson, D.C., and Parkhurst, D.L., 2004, Calculation of individual isotope
	   equilibrium constants for geochemical reactions: Geochimica et Cosmochimica
	   Acta, v. 68, no. 11, p. 2449-2465.

	   *Brief description of the program PhreeqcI.

	Charlton, S.R., and Parkhurst, D.L., 2002, PhreeqcI--A graphical user interface
	   to the geochemical model PHREEQC: U.S. Geological Survey Fact Sheet
	   FS-031-02, 2 p.


-----------------------------------------------------------------------------------
C. Downloading and unzipping the distribution file
-----------------------------------------------------------------------------------

	Follow the steps using gunzip, 7zip, or other compression program that can unzip a tar.gz file:

	Steps in extracting files                   	Explanation
	----------------------------------------    	-----------------------------------
	gunzip phreeqcrm-@VERSION@-@REVISION@.tar.gz	Uncompress the tar.gz file.
	tar -xvpof phreeqcrm-@VERSION@-@REVISION@.tar	Extract files from the tar file.

	The directory phreeqcrm-@VERSION@-@REVISION@ is created when the files are
	extracted; if this directory already exists, you may want to delete or rename
	it before extracting the files. The following directory structure is created 
	(the contents of each directory are shown to the right):

	phreeqcrm-@VERSION@-@REVISION@			Source code and configure/CMake files
	     ./config        				More configure files
	     ./database???      				Database files for PHREEQC
	     ./doc???           				Documentation files 
	     ./IPhreeqcPhast 				More source code for PhreeqcRM
	     ./Tests         				Advection example in C++, C, and Fortran

	Note:  A compiled executable is not included in the distribution.

-----------------------------------------------------------------------------------
D. Linux compiling, testing, and installing the PhreeqcRM library
-----------------------------------------------------------------------------------

	In general, to compile the software, you will need:

		(a) a C++ compiler, and
		(b) a minimal level of knowledge of configure or CMake, Make, 
			the compiler, and the Linux operating system.

	A Makefile is generated by configure, and the Makefile is used to
	compile, test, and install the software.

  D.1.  Change directory to the directory that was extracted from the tar file.

        cd phreeqcrm-@VERSION@-@REVISION@
        
  D.2.  Make a directory, for example, Release_openmp.
	
	mkdir Release_openmp
	
  D.3.  Change directory to Release_openmp.
	
	cd Release_openmp

  D.4.  Run configure

	Many of the options for configure can be seen by typing:
    
    	../configure --help
    
	Most common options:
    	--prefix=dir--specifies the directory for installation of the
    		library. Default is /usr/local, $HOME would install
    		to your home directory.
    	--with-mpi=yes--specifies that MPI (Message Passing Interface)
    		will be used for parallelization of PhreeqcRM. Default 
    		is no, in which case, OpenMP is used for parallelization.
    	--disable-openmp--specifies that OpenMP is not used. This option 
    	 	is used if you want neither kind of parallelization.
    	--enable-fortran--specifies that you have a Fortran compiler
    		and want to include Fortran files during "make check", 
    		which compiles, links, and runs programs that use 
    		PhreeqcRM.
    		
	Run configure as follows:	
    	
    	../configure [options] 

  D.5.  Compile the PhreeqcRM library
	
	make [-j 4]
	
    	Optionally, use -j n--where n is the number of compilations make runs in paralllel.
	
  D.6.  Compile and run the test case. Make check may take several minutes. 
  	The files in the directory PhreeqcRM/Tests are compiled and linked to the 
	PhreeqcRM library. The Fortran files are optionally included (--enable-fortran
	configure option). The source files in the Tests directory demonstrate
	the use of almost all PhreeqcRM features for each programming language.
	
	make check
	
  D.7.  Install the PhreeqcRM library. By default the program is installed in /usr/local/bin
	and /usr/local/share/doc/phreeqcrm???. Alternatively, the directory defined
	by --prefix in the configure command is used. 
    
    	make install
    	
	The locations of various files are given for default installation (no
	--prefix definition). 
    
	Libraries:
		/usr/local/bin/libphreeqcrm.a           (Static library)
		/usr/local/bin/phreeqcrm.so             (Shared object library)
		    
	Include files:
		C++ and C headers for user's code:
		/usr/local/include/PhreeqcRM.h          (C++ header)
		/usr/local/include/RM_interface_C.h     (C header)
		
		Optionally, user's code may need one or more of the following:
		/usr/local/include/IPhreeqc.hpp         (C++ header)
		/usr/local/include/IPhreeqc.h           (C header)
		/usr/local/include/IPhreeqcPhast.h      (C++ header)
		/usr/local/include/Var.h                (IPhreeqc VAR structure)
		/usr/local/include/IrmResult.h          (PhreeqcRM error code enum)
		/usr/local/include/IPhreeqcCallbacks.h  (IPhreeqc callbacks)
	
	Documentation??? (PDFs, README, NOTICE, phreeqcrm.txt):
		/usr/local/share/doc/phreeqcrm
    
	Databases???: 
		/usr/local/share/doc/phreeqcrm/database
		
	Fortran interfaces???
    
-----------------------------------------------------------------------------------
E. Windows compiling, testing, and installing the PhreeqcRM library
-----------------------------------------------------------------------------------

	In general, to compile the PhreeqcRM library, you will need:

		(a) CMake (build process manager)
		(b) a C++ compiler,
		(c) familiarity with the compiler and the Windows operating system.
   
  E.1. Download and install CMake

	CMake can be downloaded from http://www.cmake.org/. If you download
	the executable, CMake is installed when you execute it.

  E.2  Generate a Visual Studio solution

	CMake has many options other than Visual Studio, but here we assume
	that you are using some version of Visual Studio. On the first use
	of CMake, it may be necessary to close all Visual Studio instances. 
	
	Open CMake. Fill in the top two lines.
		Where is the source code: 
			The PhreeqcRM/src??? directory. 
		Where to build the binaries: 
			A new directory. It is suggested that the name
			of the directory include the characteristics of 
			the library that is built, including the VS
			version, OpenMP/MPI, 32/64 bit, dll/lib.
			for example, vs2012_openmp_64_lib.
	
	Click Configure. From the pulldown menu select the 64-bit version 
	for your Visual Studio. Normally, the default radio button should
	be sufficient. Click Finish. 
	
	Options in the top panel:
	(1) BUILD_SHARED_LIBS--If checked, .dll files will be
		created. If not checked, a .lib file will be created.
		Note that .dll files have a .lib file associated with them,
		but it is not the same .lib generated when
		BUILD_SHARED_LIBS is not checked.
	(2) CMAKE_INSTALL_PREFIX--Sets the the install directory where the 
		library and header files will be saved. It is suggested that
		an INSTALL directory be defined inside the directory chosen 
		for the binaries.
	(3) PHREEQCRM_BUILD_MPI--Sets whether to build libraries that use
		MPI (Message Passing Interface) for parallelization. If not 
		clicked, the libraries will use OpenMP for parallelization.
		If you check the PHREEQCRM_BUILD_MPI checkbox, you will need
		to have a version of MPI installed on your computer.
	(4) PHREEQCRM_FORTRAN_TESTING--Sets whether the test calculations
		(Tests directory) will include Fortran files (FortranAdvect
		project). A Fortran compiler is required. 
	
	Click Configure until the screen is not red. If the screen remains 
	red, you will have to deal with the CMake error messages.
	
	Click Generate.
	
	You should now have a Visual Studio solution file (.sln) in the 
	directory chosen for the binaries (second line of CMake screen). 
	At any point, you can delete the contents of the directory and
	start over. You can make solutions for different Visual Studios,
	or different parallelization in other directories.
	
  E.3. Compile and install PhreeqcRM

        Open the Visual Studio solution in the build directory defined in CMake.
        
        Build ALL_BUILD. PhreeqcRM is compiled.
        
        Build RUN_TESTS. Programs that test PhreeqcRM are compiled and run, 
        including the TestRM project, and optionally (CMake 
        PHREEQCRM_FORTRAN_TESTING option) the FortransAdvect project. 
        These projects use most of the features of PhreeqcRM.
        
        Build INSTALL. The libraries, header files, and Fortran module source
        files are installed in the install directory selected in CMake.
        
        	Libraries that can be installed:
        		PhreeqcRM.lib
        		PhreeqcRMd.lib
        		PhreeqcRM.dll
        		PhreeqcRMd.dll
        		
        	Include files:
			C++ and C headers for user's code:	
				PhreeqcRM.h          (C++ header)
				RM_interface_C.h     (C header)

			Optionally, may need one or more of the following:
				IPhreeqc.hpp         (C++ header)
				IPhreeqc.h           (C header)
				IPhreeqcPhast.h      (C++ header)
				Var.h                (IPhreeqc VAR structure)
				IrmResult.h          (PhreeqcRM error code enum)
				IPhreeqcCallbacks.h  (IPhreeqc callbacks)
        		
        	Src files???:
        		To use PhreeqcRM in Fortran:
        			RM_interface.F90
        		
        		If IPhreeqc methods are called in the user's program:
        			IPhreeqc_interface.F90

-----------------------------------------------------------------------------------        			
F. Using the library
-----------------------------------------------------------------------------------

  F.1. C++
     
     In C++, you will include the PhreeqcRM.h header file and 
     create a PhreeqcRM instance (multithreaded in these examples):

     #include "PhreeqcRM.h"
     PhreeqcRM  my_phreeqcrm(nxyz, nthreads);

     or

     PhreeqcRM * myphreeqc_ = new PhreeqcRM(nxyz, nthreads);

     You will need to link to the PhreeqcRM.lib that you installed.

  F.2. C

     In C, you will include the C interface header, RM_interface_C.h, and
     create a PhreeqcRM instance.

     #include "RM_interface_C.h"
     int rm_id;
     rm_id = RM_create(nxyz, nthreads);

     You will need to link to the PhreeqcRM.lib that you installed.

  F.3. Fortran

     In Fortran, you will need to include the source file RM_interface.F90
     in your project files. This file defines the PhreeqcRM Fortran module.

     USE PhreeqcRM
     integer rm_id
     rm_id = RM_create(nxyz, nthreads)

     You will need to link to the PhreeqcRM.lib that you installed.

     If you use any IPhreeqc methods directly (optional, example in 
     Advection_f90.F90), you need to include IPhreeqc_interface.F90 in your 
     project files, which defines the IPhreeqc Fortran module.

-----------------------------------------------------------------------------------
G. CONTACTS
-----------------------------------------------------------------------------------

Inquiries about this software distribution should be directed to:

e-mail:  h2osoft@usgs.gov or dlpark@usgs.gov

