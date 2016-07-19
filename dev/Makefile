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
	SRC = /home/helmes/Work/toric_code_3D_ALPS/dev
	BIN = /home/helmes/Work/toric_code_3D_ALPS/bin
	CXX = /home/helmes/intel2016/bin/icpc
	CXXFLAGS = -O2 -DNDEBUG -march=corei7 -std=c++11
endif
ifeq ($(host), johannes-thinkpad)
	ALPS_ROOT = /opt/alps
	include $(ALPS_ROOT)/share/alps/include.mk
	SRC = /home/johannes/Work/toric_code_3D_ALPS/dev
	BIN = /home/johannes/Work/toric_code_3D_ALPS/bin
	CXX = /usr/bin/c++
	CXXFLAGS = -O2 -std=c++11
endif
ifeq ($(host), cheops0)
	ALPS_ROOT = /projects/ag-trebst/alps/alps-svn-intel-intelmpi
	include $(ALPS_ROOT)/share/alps/include.mk
	SRC = /home/jhelmes/Work/toric_code_3D_ALPS/dev
	BIN = /home/jhelmes/Work/toric_code_3D_ALPS/bin
	CXXFLAGS = -O2 -DNDEBUG -march=native -std=c++11 # -gcc-name=/opt/rrzk/compiler/gnu/4.8.2/bin/gcc
endif



all: toriccode_dev xcube_dev haahcode_dev
	@echo $(host)
	@echo $(BIN)
	@echo $(LDFLAGS)

#######################  Patterns ################################
$(BIN)/%_dev.o: $(SRC)/%.cpp $(SRC)/%.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@ $(COMFLAG)

###################### End Patterns ##############################


toriccode_dev: $(BIN)/site_dev.o $(BIN)/measurer_dev.o $(BIN)/updater_dev.o $(BIN)/toriccode_dev.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN)/$@ $(BIN)/site_dev.o $(BIN)/measurer_dev.o $(BIN)/updater_dev.o $(BIN)/toriccode_dev.o $(LDFLAGS) $(LIBS) $(LINKFLAG) $(COMFLAG)

xcube_dev: $(BIN)/site_dev.o $(BIN)/measurer_dev.o $(BIN)/updater_dev.o $(BIN)/xcube_dev.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN)/$@ $(BIN)/site_dev.o $(BIN)/measurer_dev.o $(BIN)/updater_dev.o $(BIN)/xcube_dev.o $(LDFLAGS) $(LIBS) $(LINKFLAG) $(COMFLAG)

haahcode_dev: $(BIN)/site_dev.o $(BIN)/measurer_dev.o $(BIN)/updater_dev.o $(BIN)/haahcode_dev.o 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(BIN)/$@ $(BIN)/site_dev.o $(BIN)/measurer_dev.o $(BIN)/updater_dev.o $(BIN)/haahcode_dev.o $(LDFLAGS) $(LIBS) $(LINKFLAG) $(COMFLAG)

clean:
	$(RM) $(RMFLAGS) $(BIN)/*_dev.o toriccode_dev xcube_dev haah_dev
