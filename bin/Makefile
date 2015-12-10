#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#
#  Makefile for ALPS code on local thp-machine
#       Johannes Helmes


ALPS_ROOT = /scratch/alps

include $(ALPS_ROOT)/share/alps/include.mk

CXXFLAGS = -O2 -DNDEBUG -march=corei7

all: FTRandomIntegration FTRandomSwitch

toric3DFTSingleSpin.o: ../src/toric3DFTSingleSpin.cpp ../src/toric3DFTSingleSpin.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++0x -c ../src/toric3DFTSingleSpin.cpp -o toric3DFTSingleSpin.o $(COMFLAG)

FTSPRandomSwitch.o: ../src/FTSPRandomSwitch.cpp ../src/FTSPRandomSwitch.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++0x -c ../src/FTSPRandomSwitch.cpp -o FTSPRandomSwitch.o $(COMFLAG)
FTSPRandomSwitch: toric3DFTSingleSpin.o FTSPRandomSwitch.o mainFTSPRandomSwitch.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ mainFTSPRandomSwitch.o FTSPRandomSwitch.o toric3DFTSingleSpin.o $(LDFLAGS) $(LIBS) $(LINKFLAG)

FTSPRandomIntegration.o: ../src/FTSPRandomIntegration.cpp ../src/FTSPRandomIntegration.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++0x -c ../src/FTSPRandomIntegration.cpp -o FTSPRandomIntegration.o $(COMFLAG)
FTSPRandomIntegration: toric3DFTSingleSpin.o FTSPRandomIntegration.o mainFTSPRandomIntegration.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ mainFTSPRandomIntegration.o FTSPRandomIntegration.o toric3DFTSingleSpin.o $(LDFLAGS) $(LIBS) $(LINKFLAG)

mainFTSPRandomSwitch.o: ../src/mainFTSPRandomSwitch.C
	$(CXX) $(CPPFLAGS)  $(CXXFLAGS) -std=c++0x -c ../src/mainFTSPRandomSwitch.C
mainFTSPRandomIntegration.o: ../src/mainFTSPRandomIntegration.C
	$(CXX) $(CPPFLAGS)  $(CXXFLAGS) -std=c++0x -c ../src/mainFTSPRandomIntegration.C


clean:
	$(RM) $(RMFLAGS) *.o FTRandomIntegration FTRandomSwitch toric3DFiniteTemp FTSPRandomIntegration FTSPMembraneRandomIntegration FTSPRandomSwitch toric3DFTSingleSpin toric3DFTMembraneSingleSpin
