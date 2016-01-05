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
	SRC = /home/helmes/Work/toric_code_3D_ALPS/src
	BIN = /home/helmes/Work/toric_code_3D_ALPS/bin

	CXXFLAGS = -O2 -DNDEBUG -march=corei7
endif
ifeq ($(host), cheops0)
	ALPS_ROOT = /projects/ag-trebst/alps/alps-svn-intel-intelmpi
	SRC = /home/jhelmes/Work/toric_code_3D_ALPS/src
	BIN = /home/jhelmes/Work/toric_code_3D_ALPS/bin

	CXXFLAGS = -O2 -DNDEBUG -march=native
endif


include $(ALPS_ROOT)/share/alps/include.mk


all: FTSPRandomIntegration FTSPRandomSwitch
	@echo $(host)
	@echo $(BIN)

#######################  Patterns ################################
$(BIN)/%.o: $(SRC)/%.cpp $(SRC)/%.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++0x -c $< -o $@ $(COMFLAG)

$(BIN)/main%.o: $(SRC)/main%.C
	$(CXX) $(CPPFLAGS)  $(CXXFLAGS) -std=c++0x -c $< -o $@ 
###################### End Patterns ##############################

FTSPRandomSwitch: $(BIN)/toric3DFTSingleSpin.o $(BIN)/FTSPRandomSwitch.o $(BIN)/mainFTSPRandomSwitch.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN)/$@ $(BIN)/mainFTSPRandomSwitch.o $(BIN)/FTSPRandomSwitch.o $(BIN)/toric3DFTSingleSpin.o $(LDFLAGS) $(LIBS) $(LINKFLAG)

FTSPRandomIntegration: $(BIN)/toric3DFTSingleSpin.o $(BIN)/FTSPRandomIntegration.o $(BIN)/mainFTSPRandomIntegration.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN)/$@ $(BIN)/mainFTSPRandomIntegration.o $(BIN)/FTSPRandomIntegration.o $(BIN)/toric3DFTSingleSpin.o $(LDFLAGS) $(LIBS) $(LINKFLAG)

clean:
	$(RM) $(RMFLAGS) $(BIN)/*.o FTRandomIntegration FTRandomSwitch toric3DFiniteTemp FTSPRandomIntegration FTSPMembraneRandomIntegration FTSPRandomSwitch toric3DFTSingleSpin toric3DFTMembraneSingleSpin
