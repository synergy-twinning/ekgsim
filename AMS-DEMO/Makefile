COMPILEFLAGS = -g -c
LINKFLAGS = -O3 -s -L../include -L../../lib -L../copyOfLibs -lrt -ldl -lIni -lArguments -lRandom
CPPFLAGS = -O3 -s -I. -I.. -I../include -I../../include -I../copyOfLibs -std=c++0x
MPICXX = mpic++.mpich
MPICXX = mpic++
CXX = c++
GENERAL_FILES = ParallelFramework.h ParallelNumericOptimizer.h GeneralOptimizationAlgorithm.h GeneralOptimizationAlgorithm.cpp pDemo.h pDe.h DeBase.h VectorArithmetics.h internalEval.h TypeWrapper.h Initializer.h
H_FILES = Array.h ExternalEvaluation.h Initializer.h NumericOptimizer.h pDemo.h TypeWrapper.h DeBase.h GeneralOptimizationAlgorithm.h internalEval.h ParallelFramework.h utilities.h De.h HyperVolume.h jobDistributer.h ParallelNumericOptimizer.h TestFunctions.h VectorArithmetics.h Demo.h Individual.h mpiWrapper.h pDe.h timer.h

all: directories AMS-DEMO DEMO
	cp AMS-DEMO bin/
	cp DEMO bin/
	
directories: bin

bin:
	mkdir -p bin
	
echo:
	\# $(MPICXX) ${COMPILEFLAGS} ${CPPFLAGS} main.cpp


AMS-DEMO: main.o utilities.o GeneralOptimizationAlgorithm.o
	$(MPICXX) -o AMS-DEMO GeneralOptimizationAlgorithm.o utilities.o main.o -I/usr/mpich/include $(LINKFLAGS)

GeneralOptimizationAlgorithm.o: $(GENERAL_FILES) 
	$(MPICXX) ${COMPILEFLAGS} ${CPPFLAGS} GeneralOptimizationAlgorithm.cpp

utilities.o: utilities.h timer.h utilities.cpp
	g++ ${COMPILEFLAGS} ${CPPFLAGS} utilities.cpp

main.o: main.cpp ${PCH_FILES}
	$(MPICXX) ${COMPILEFLAGS} ${CPPFLAGS} main.cpp
	

DEMO: main.nompi.o utilities.o GeneralOptimizationAlgorithm.nompi.o
	$(CXX) -o DEMO GeneralOptimizationAlgorithm.nompi.o utilities.o main.nompi.o -I/usr/mpich/include $(LINKFLAGS) 

GeneralOptimizationAlgorithm.nompi.o: $(GENERAL_FILES)
	$(CXX) ${COMPILEFLAGS} ${CPPFLAGS} -DNO_MPI GeneralOptimizationAlgorithm.cpp -o GeneralOptimizationAlgorithm.nompi.o

main.nompi.o: main.cpp ${PCH_FILES}
	$(CXX) ${COMPILEFLAGS} ${CPPFLAGS} -DNO_MPI main.cpp -o main.nompi.o 

clean:
	rm -rf *.o AMS-DEMO DEMO


