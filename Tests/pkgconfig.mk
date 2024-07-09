PKGCONFIG= pkg-config
CPPFLAGS:= $(shell $(PKGCONFIG) --cflags phreeqcrm)
LDFLAGS:= $(shell $(PKGCONFIG) --libs phreeqcrm)

OBJS =\
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
	TestAllMethods_c.o\
	TestAllMethods_cpp.o\
	WriteYAMLFile_cpp.o\
	WriteYAMLFile_cpp_test.o

CLEANFILES =\
	AdvectBMI_cpp.chem.txt\
	AdvectBMI_cpp.log.txt\
	AdvectBMI_cpp.yaml\
	AdvectBMI_cpp_test.chem.txt\
	AdvectBMI_cpp_test.log.txt\
	AdvectBMI_cpp_test.yaml\
	Advect_c.chem.txt\
	Advect_c.dmp\
	Advect_c.log.txt\
	Advect_c_utility.txt\
	Advect_cpp.chem.txt\
	Advect_cpp.dmp\
	Advect_cpp.log.txt\
	Advect_cpp_units_utility.txt\
	Advect_cpp_units_worker.chem.txt\
	Advect_cpp_units_worker.log.txt\
	Advect_cpp_utility.txt\
	Gas_c.chem.txt\
	Gas_c.log.txt\
	Gas_cpp.chem.txt\
	Gas_cpp.log.txt\
	SimpleAdvect_c.chem.txt\
	SimpleAdvect_c.log.txt\
	SimpleAdvect_cpp.chem.txt\
	SimpleAdvect_cpp.log.txt\
	Species_c.chem.txt\
	Species_c.dmp\
	Species_c.log.txt\
	Species_c_utility.txt\
	Species_cpp.chem.txt\
	Species_cpp.dmp\
	Species_cpp.log.txt\
	Species_cpp_utility.txt\
	TestAllMethods_c.chem.txt\
	TestAllMethods_c.log.txt\
	TestAllMethods_c.yaml\
	TestAllMethods_cpp.chem.txt\
	TestAllMethods_cpp.dump\
	TestAllMethods_cpp.log.txt\
	TestAllMethods_cpp.yaml

# Executable
TARGET := pkgconfig_test

# Linking rule
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $(TARGET)

clean :
	rm -f $(OBJS) $(CLEANFILES) pkgconfig_test
