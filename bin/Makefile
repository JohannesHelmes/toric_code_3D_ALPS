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

	CXXFLAGS = -O2 -DNDEBUG -march=corei7
endif
ifeq ($(host), cheops0)
	ALPS_ROOT = /projects/ag-trebst/alps/alps-svn-intel-intelmpi
	SRC = /home/jhelmes/Work/toric_code_3D_ALPS/src

	CXXFLAGS = -O2 -DNDEBUG -march=native
endif


include $(ALPS_ROOT)/share/alps/include.mk


all: FTRandomIntegration FTRandomSwitch
	@echo $(host)

#######################  Patterns ################################
%.o:: $(SRC)/%.cpp $(SRC)/%.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++0x -c $< -o $@ $(COMFLAG)

main%.o: $(SRC)/main%.C
	$(CXX) $(CPPFLAGS)  $(CXXFLAGS) -std=c++0x -c $<
###################### End Patterns ##############################

FTSPRandomSwitch: toric3DFTSingleSpin.o FTSPRandomSwitch.o mainFTSPRandomSwitch.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ mainFTSPRandomSwitch.o FTSPRandomSwitch.o toric3DFTSingleSpin.o $(LDFLAGS) $(LIBS) $(LINKFLAG)

FTSPRandomIntegration: toric3DFTSingleSpin.o FTSPRandomIntegration.o mainFTSPRandomIntegration.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ mainFTSPRandomIntegration.o FTSPRandomIntegration.o toric3DFTSingleSpin.o $(LDFLAGS) $(LIBS) $(LINKFLAG)

clean:
	$(RM) $(RMFLAGS) *.o FTRandomIntegration FTRandomSwitch toric3DFiniteTemp FTSPRandomIntegration FTSPMembraneRandomIntegration FTSPRandomSwitch toric3DFTSingleSpin toric3DFTMembraneSingleSpin
