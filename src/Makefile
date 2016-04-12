#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#
#  Makefile for ALPS code on different machines
#       Johannes Helmes, 2016


host=$(shell hostname)

ifeq ($(host), l29)
	ALPS_ROOT = /scratch/alps
	include $(ALPS_ROOT)/share/alps/include.mk
	SRC = /home/helmes/Work/toric_code_3D_ALPS/src
	BIN = /home/helmes/Work/toric_code_3D_ALPS/bin
	CXX = /home/helmes/intel/bin/icpc
	CXXFLAGS = -O2 -DNDEBUG -march=corei7 -std=c++11 
endif
ifeq ($(host), johannes-thinkpad)
	ALPS_ROOT = /opt/alps
	include $(ALPS_ROOT)/share/alps/include.mk
	SRC = /home/johannes/Work/toric_code_3D_ALPS/src
	BIN = /home/johannes/Work/toric_code_3D_ALPS/bin
	CXX = /usr/bin/c++
	CXXFLAGS = -O2 -std=c++11
endif
ifeq ($(host), cheops0)
	ALPS_ROOT = /projects/ag-trebst/alps/alps-svn-intel-intelmpi
	include $(ALPS_ROOT)/share/alps/include.mk
	SRC = /home/jhelmes/Work/toric_code_3D_ALPS/src
	BIN = /home/jhelmes/Work/toric_code_3D_ALPS/bin
	CXXFLAGS = -O2 -DNDEBUG -march=native -std=c++0x 
endif



all: toriccode
	@echo $(host)
	@echo $(BIN)

#######################  Patterns ################################
$(BIN)/%.o: $(SRC)/%.cpp $(SRC)/%.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@ $(COMFLAG)

$(BIN)/main%.o: $(SRC)/main%.C
	$(CXX) $(CPPFLAGS)  $(CXXFLAGS) -c $< -o $@ 
###################### End Patterns ##############################

FTSPRandomSwitch: $(BIN)/toric3DFTSingleSpin.o $(BIN)/FTSPRandomSwitch.o $(BIN)/mainFTSPRandomSwitch.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN)/$@ $(BIN)/mainFTSPRandomSwitch.o $(BIN)/FTSPRandomSwitch.o $(BIN)/toric3DFTSingleSpin.o $(LDFLAGS) $(LIBS) $(LINKFLAG)

FTSPRandomIntegration: $(BIN)/toric3DFTSingleSpin.o $(BIN)/FTSPRandomIntegration.o $(BIN)/mainFTSPRandomIntegration.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN)/$@ $(BIN)/mainFTSPRandomIntegration.o $(BIN)/FTSPRandomIntegration.o $(BIN)/toric3DFTSingleSpin.o $(LDFLAGS) $(LIBS) $(LINKFLAG)

FTSPMembraneRandomIntegration: $(BIN)/site.o $(BIN)/measurer.o $(BIN)/updater.o $(BIN)/toric3DFTMembraneSingleSpin.o $(BIN)/FTSPMembraneRandomIntegration.o $(BIN)/mainFTSPMembraneRandomIntegration.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN)/$@ $(BIN)/site.o $(BIN)/measurer.o $(BIN)/updater.o $(BIN)/mainFTSPMembraneRandomIntegration.o $(BIN)/FTSPMembraneRandomIntegration.o $(BIN)/toric3DFTMembraneSingleSpin.o $(LDFLAGS) $(LIBS) $(LINKFLAG) $(COMFLAG)

toriccode: $(BIN)/site.o $(BIN)/measurer.o $(BIN)/updater.o $(BIN)/toriccode.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN)/$@ $(BIN)/site.o $(BIN)/measurer.o $(BIN)/updater.o $(BIN)/toriccode.o $(LDFLAGS) $(LIBS) $(LINKFLAG) $(COMFLAG)

clean:
	$(RM) $(RMFLAGS) $(BIN)/*.o FTRandomIntegration FTRandomSwitch toric3DFiniteTemp FTSPRandomIntegration FTSPMembraneRandomIntegration FTSPRandomSwitch toric3DFTSingleSpin toric3DFTMembraneSingleSpin
