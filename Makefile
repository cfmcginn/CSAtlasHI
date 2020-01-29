CXX = g++
#O3 for max optimization (go to 0 for debug)
CXXFLAGS = -Wall -Werror -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

ifndef FJCONTRIB__HOME
$(error FJCONTRIB__HOME is not set at all. Please set this env. variable to point to your contrib package)
endif

INCLUDE=-I$(PWD)
LIB=-L$(PWD)/lib
ROOT=`root-config --cflags --glibs`

FASTJET=`fastjet-config --cxxflags  --libs --plugins --runpath`
FJCONTRIB=-I$(FJCONTRIB__HOME)/include -L$(FJCONTRIB__HOME)/lib -lConstituentSubtractor  

PYTHIA8=-I$(PYTHIA8PATH)/include -O2 -pedantic -W -Wall -Wshadow -fPIC -L$(PYTHIA8PATH)/lib -Wl,-rpath,$(PYTHIA8PATH)/lib -lpythia8 -ldl

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_LIB=mkdir -p $(PWD)/lib
MKDIR_OBJ=mkdir -p $(PWD)/obj
MKDIR_OUTPUT=mkdir -p $(PWD)/output
MKDIR_PDF=mkdir -p $(PWD)/pdfDir

all: mkdirBin mkdirLib mkdirObj mkdirOutput mkdirPdf obj/globalDebugHandler.o obj/checkMakeDir.o obj/constituentBuilder.o obj/rhoBuilder.o obj/configParser.o obj/centralityFromInput.o lib/libCSATLAS.so bin/makeClusterTree.exe bin/makeClusterHist.exe bin/plotClusterHist.exe bin/deriveSampleWeights.exe bin/deriveCentWeights.exe

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

obj/configParser.o: src/configParser.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/configParser.C -o obj/configParser.o $(INCLUDE) $(ROOT)

obj/centralityFromInput.o: src/centralityFromInput.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/centralityFromInput.C -o obj/centralityFromInput.o $(INCLUDE) $(ROOT)

lib/libCSATLAS.so:
	$(CXX) $(CXXFLAGS) -fPIC -shared -o lib/libCSATLAS.so obj/checkMakeDir.o obj/globalDebugHandler.o obj/constituentBuilder.o obj/rhoBuilder.o obj/configParser.o obj/centralityFromInput.o $(FASTJET) $(ROOT) $(INCLUDE)

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
	rm -f lib/*.o
	rm -rf lib
	rm -f obj/*.o
	rm -rf obj
	rm -f src/*~
	rm -f src/#*#
