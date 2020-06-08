cCXX = g++
#O3 for max optimization (go to 0 for debug)
CXXFLAGS = -Wall -Werror -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

define QTDIRERR
 QTDIR is not set at all. Please set this environment variable to point to your build - this should be either
export QTDIR=$(PWD)
or
source setEnv.sh
if you have made appropriate changes.
For more, see README for full setup recommendations
endef

define FJCONTRIBERR
 FJCONTRIB__HOME is not set at all. Please set this environment variable to point to your contrib package. If local, try your fastjet-install, i.e
export FJCONTRIB__HOME=`fastjet-config --cxxflags`
If lxplus or acf/rcf, try
lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt fjcontrib" 
and see README for full setup recommendations
endef

ifndef QTDIR
$(error "$(QTDIRERR)")	
endif

ifndef FJCONTRIB__HOME
$(error "$(FJCONTRIBERR)")	
endif

INCLUDE=-I$(QTDIR)
LIB=-L$(QTDIR)/lib
ROOT=`root-config --cflags --glibs`

FASTJET=`fastjet-config --cxxflags  --libs --plugins --runpath`
FJCONTRIB=-I$(FJCONTRIB__HOME)/include -L$(FJCONTRIB__HOME)/lib -lConstituentSubtractor  

PYTHIA8=-I$(PYTHIA8PATH)/include -O2 -pedantic -W -Wall -Wshadow -fPIC -L$(PYTHIA8PATH)/lib -Wl,-rpath,$(PYTHIA8PATH)/lib -lpythia8 -ldl

MKDIR_BIN=mkdir -p $(QTDIR)/bin
MKDIR_LIB=mkdir -p $(QTDIR)/lib
MKDIR_OBJ=mkdir -p $(QTDIR)/obj
MKDIR_OUTPUT=mkdir -p $(QTDIR)/output
MKDIR_PDF=mkdir -p $(QTDIR)/pdfDir

all: mkdirBin mkdirLib mkdirObj mkdirOutput mkdirPdf obj/checkMakeDir.o obj/constituentBuilder.o obj/globalDebugHandler.o obj/rhoBuilder.o obj/sampleHandler.o obj/configParser.o obj/centralityFromInput.o obj/towerWeightTwol.o lib/libCSATLAS.so bin/makeClusterTree.exe bin/makeClusterHist.exe bin/plotClusterHist.exe bin/deriveSampleWeights.exe bin/deriveCentWeights.exe bin/validateRho.exe bin/validateRhoHist.exe bin/validateRhoPlot.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirLib:
	$(MKDIR_LIB)

mkdirObj:
	$(MKDIR_OBJ)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

bin/constituentTest.exe: src/constituentTest.C
	$(CXX) $(CXXFLAGS) src/constituentTest.C $(ROOT) $(PYTHIA8) $(FASTJET) $(FJCONTRIB) $(INCLUDE) -o bin/constituentTest.exe

obj/checkMakeDir.o: src/checkMakeDir.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/checkMakeDir.C -o obj/checkMakeDir.o $(INCLUDE)

obj/globalDebugHandler.o: src/globalDebugHandler.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/globalDebugHandler.C -o obj/globalDebugHandler.o $(ROOT) $(INCLUDE)

obj/constituentBuilder.o: src/constituentBuilder.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/constituentBuilder.C -o obj/constituentBuilder.o $(FASTJET) $(ROOT) $(INCLUDE)

obj/rhoBuilder.o: src/rhoBuilder.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/rhoBuilder.C -o obj/rhoBuilder.o $(INCLUDE)

obj/sampleHandler.o: src/sampleHandler.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/sampleHandler.C -o obj/sampleHandler.o $(INCLUDE)


obj/configParser.o: src/configParser.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/configParser.C -o obj/configParser.o $(INCLUDE) $(ROOT)

obj/centralityFromInput.o: src/centralityFromInput.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/centralityFromInput.C -o obj/centralityFromInput.o $(INCLUDE) $(ROOT)

obj/towerWeightTwol.o: src/towerWeightTwol.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/towerWeightTwol.C -o obj/towerWeightTwol.o $(INCLUDE) $(ROOT)

lib/libCSATLAS.so:
	$(CXX) $(CXXFLAGS) -fPIC -shared -o lib/libCSATLAS.so obj/checkMakeDir.o obj/globalDebugHandler.o obj/constituentBuilder.o obj/rhoBuilder.o obj/configParser.o obj/centralityFromInput.o obj/sampleHandler.o $(FASTJET) $(ROOT) $(INCLUDE)

bin/makeClusterTree.exe: src/makeClusterTree.C
	$(CXX) $(CXXFLAGS) src/makeClusterTree.C -o bin/makeClusterTree.exe $(FJCONTRIB) $(FASTJET) $(ROOT) $(INCLUDE) $(LIB) -lCSATLAS

bin/clusterToCS.exe: src/clusterToCS.C
	$(CXX) $(CXXFLAGS) src/clusterToCS.C $(ROOT) $(FJCONTRIB) $(FASTJET) $(INCLUDE) -fopenmp -o bin/clusterToCS.exe

bin/makeClusterHist.exe: src/makeClusterHist.C
	$(CXX) $(CXXFLAGS) src/makeClusterHist.C $(ROOT) $(INCLUDE) $(LIB) -lCSATLAS -o bin/makeClusterHist.exe

bin/plotClusterHist.exe: src/plotClusterHist.C
	$(CXX) $(CXXFLAGS) src/plotClusterHist.C $(ROOT) $(INCLUDE) $(LIB) -lCSATLAS -o bin/plotClusterHist.exe

bin/deriveSampleWeights.exe: src/deriveSampleWeights.C
	$(CXX) $(CXXFLAGS) src/deriveSampleWeights.C $(INCLUDE) $(ROOT) $(LIB) -lCSATLAS -o bin/deriveSampleWeights.exe

bin/deriveCentWeights.exe: src/deriveCentWeights.C
	$(CXX) $(CXXFLAGS) src/deriveCentWeights.C $(INCLUDE) $(ROOT) $(LIB) -lCSATLAS -o bin/deriveCentWeights.exe

bin/validateRho.exe: src/validateRho.C
	$(CXX) $(CXXFLAGS) src/validateRho.C $(INCLUDE) $(ROOT) $(FASTJET) $(LIB) -lCSATLAS -o bin/validateRho.exe

bin/validateRhoHist.exe: src/validateRhoHist.C
	$(CXX) $(CXXFLAGS) src/validateRhoHist.C $(INCLUDE) $(ROOT) $(FASTJET) $(LIB) -lCSATLAS -o bin/validateRhoHist.exe

bin/validateRhoPlot.exe: src/validateRhoPlot.C
	$(CXX) $(CXXFLAGS) src/validateRhoPlot.C $(INCLUDE) $(ROOT) $(FASTJET) $(LIB) -lCSATLAS -o bin/validateRhoPlot.exe

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
	rm -f lib/*.so
	rm -rf lib
	rm -f obj/*.o
	rm -rf obj
	rm -f src/*~
	rm -f src/#*#
