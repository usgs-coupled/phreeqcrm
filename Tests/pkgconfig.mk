PKGCONFIG= pkg-config
CPPFLAGS:= $(shell $(PKGCONFIG) --cflags phreeqcrm)
LDLIBS:= $(shell $(PKGCONFIG) --libs phreeqcrm)

objects =\
	AdvectBMI_cpp.o\
	AdvectBMI_cpp_test.o\
	Advect_c.o\
	Advect_cpp.o\
	Gas_c.o\
	Gas_cpp.o\
	main.o\
	SimpleAdvect_c.o\
	SimpleAdvect_cpp.o\
	Species_c.o\
	Species_cpp.o\
	TestAllMethods_cpp.o\
	WriteYAMLFile_cpp.o\
	WriteYAMLFile_cpp_test.o

cleanfiles =\
	Advect_c_utility.txt\
	Advect_c.chem.txt\
	Advect_c.dmp\
	Advect_c.log.txt\
	Advect_cpp_units_utility.txt\
	Advect_cpp_units_worker.chem.txt\
	Advect_cpp_units_worker.log.txt\
	Advect_cpp_units_worker.txt\
	Advect_cpp_utility.txt\
	Advect_cpp.chem.txt\
	Advect_cpp.dmp\
	Advect_cpp.log.txt\
	AdvectBMI_cpp_units_utility.txt\
	AdvectBMI_cpp_units_worker.txt\
	AdvectBMI_cpp_utility.txt\
	AdvectBMI_cpp.chem.txt\
	AdvectBMI_cpp.dmp\
	AdvectBMI_cpp.log.txt\
	AdvectBMI_cpp.yaml\
	Advectcpp_utility.txt\
	Advectcpp.dmp\
	Gas_c_utility.txt\
	Gas_c.chem.txt\
	Gas_c.dmp\
	Gas_c.log.txt\
	Gas_cpp_utility.txt\
	Gas_cpp.chem.txt\
	Gas_cpp.dmp\
	Gas_cpp.log.txt\
	SimpleAdvect_c.chem.txt\
	SimpleAdvect_c.log.txt\
	SimpleAdvect_cpp.chem.txt\
	SimpleAdvect_cpp.log.txt\
	Species_c_utility.txt\
	Species_c.chem.txt\
	Species_c.dmp\
	Species_c.log.txt\
	Species_cpp_utility.txt\
	Species_cpp.chem.txt\
	Species_cpp.dmp\
	Species_cpp.log.txt\
	Units_Worker.chem.txt\
	Units_Worker.log.txt\
	Utility_c.txt\
	Utility_cpp.out


pkgconfig_test : $(objects)
	$(CXX) -o pkgconfig_test $(objects) $(LDLIBS)

clean :
	rm -f $(objects) $(cleanfiles) pkgconfig_test
