CXX = g++
#O3 for max optimization (go to 0 for debug)
CXXFLAGS = -Wall -Werror -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

INCLUDE=-I $(PWD)
ROOT=`root-config --cflags --glibs`

FASTJETPATH=/usatlas/u/cfmcginn/Packages/FastJet/fastjet-install
FASTJETCS=/usatlas/u/cfmcginn/Packages/FastJet/fjcontrib-1.042/ConstituentSubtractor
FASTJET=`$(FASTJETPATH)/bin/fastjet-config --cxxflags  --libs --plugins --rpath`
#FASTJET=-I/home/cfmcginn/Packages/FastJet/fastjet-install/include -Wl,-rpath,/home/cfmcginn/Packages/FastJet/fastjet-install/lib -lm -L/home/cfmcginn/Packages/FastJet/fastjet-install/lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone
FJCONTRIB=-lConstituentSubtractor  

PYTHIA8=-I$(PYTHIA8PATH)/include -O2 -pedantic -W -Wall -Wshadow -fPIC -L$(PYTHIA8PATH)/lib -Wl,-rpath,$(PYTHIA8PATH)/lib -lpythia8 -ldl

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_OUTPUT=mkdir -p $(PWD)/output
MKDIR_PDF=mkdir -p $(PWD)/pdfDir

all: mkdirBin mkdirPdf mkdirOutput bin/clusterToCS.exe bin/makeClusterHist.exe bin/plotClusterHist.exe bin/deriveSampleWeights.exe bin/deriveCentWeights.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

bin/constituentTest.exe: src/constituentTest.C
	$(CXX) $(CXXFLAGS) src/constituentTest.C $(ROOT) $(PYTHIA8) -ltbb $(FASTJET) $(FJCONTRIB) $(INCLUDE) -o bin/constituentTest.exe

bin/clusterToCS.exe: src/clusterToCS.C $(FASTJETCS)/ConstituentSubtractor.cc
	$(CXX) $(CXXFLAGS) src/clusterToCS.C $(FASTJETCS)/ConstituentSubtractor.cc $(ROOT) $(FASTJET) $(FJCONTRIB) $(INCLUDE) -fopenmp -o bin/clusterToCS.exe

bin/makeClusterHist.exe: src/makeClusterHist.C
	$(CXX) $(CXXFLAGS) src/makeClusterHist.C $(ROOT) $(INCLUDE) -o bin/makeClusterHist.exe

bin/plotClusterHist.exe: src/plotClusterHist.C
	$(CXX) $(CXXFLAGS) src/plotClusterHist.C $(ROOT) $(INCLUDE) -o bin/plotClusterHist.exe

bin/deriveSampleWeights.exe: src/deriveSampleWeights.C
	$(CXX) $(CXXFLAGS) src/deriveSampleWeights.C $(INCLUDE) $(ROOT) -o bin/deriveSampleWeights.exe

bin/deriveCentWeights.exe: src/deriveCentWeights.C
	$(CXX) $(CXXFLAGS) src/deriveCentWeights.C $(INCLUDE) $(ROOT) -o bin/deriveCentWeights.exe

clean:
	rm -f ./*~
	rm -f ./#*#
	rm -f bash/*~
	rm -f bash/#*#
	rm -f bin/*.exe
	rm -rf bin
	rm -f configs/*~
	rm -f configs/#*#
	rm -f include/*~
	rm -f include/#*#
	rm -f input/*~
	rm -f input/#*#
	rm -f src/*~
	rm -f src/#*#
