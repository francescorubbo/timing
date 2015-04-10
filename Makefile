# --------------------------------------------- #
# Makefile for Recluster code                        #
# Pascal Nef, March 6th 2014                    #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS =   -O2 -Wall -std=c++0x -g

.PHONY: clean debug all

all: setup Timing

setup:
	mkdir -p lib

Timing:  lib/Timing.so lib/TimingAnalysis.so
	$(CXX) lib/Timing.so lib/TimingAnalysis.so -o $@ \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs` -lEG -lEGPythia8 \
	-I./include -L./lib \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lSubjetJVF  -lVertexJets \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options 

lib/Timing.so: src/Timing.C lib/TimingAnalysis.so   
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I./include -L./lib \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

lib/TimingAnalysis.so : src/TimingAnalysis.cc include/TimingAnalysis.h 
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets \
	-I./include \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags --libs` 

clean:
	rm -rf Timing
	rm -rf lib
	rm -f *~

install:
	install Timing -t ${HOME}/local/bin
	install setupTiming.sh -t ${HOME}/local/bin
	install scripts/Timing.sh -t ${HOME}/local/bin