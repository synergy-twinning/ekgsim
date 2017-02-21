LinkerName             :=/usr/bin/mpic++
SharedObjectLinkerName :=/usr/bin/mpiAR       := /usr/bin/ar rcu
IntermediateDirectory  :=./Release
Preprocessors          :=-DNDEBUG 
LinkOptions            :=
LibPath                := -L.
MakeDirCommand         :=mkdir -p
ObjectsFileList        :="copyOfLibs.txt"
IncludePath            := -I. 

AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++
CC       := /usr/bin/gcc
CXXFLAGS :=  -O2 -std=c++11 -Wall $(Preprocessors)
CFLAGS   :=  -O2 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/asc++ -shared -fPIC



##
## Main Build Targets 
##
.PHONY: all clean MakeIntermediateDirs
all: libIni.a libRandom.a libArguments.a

## libIni
libIni.a: $(IntermediateDirectory)/Ini.o
	@$(MakeDirCommand) $(@D)
	$(AR) libIni.a $(IntermediateDirectory)/Ini.o
	
$(IntermediateDirectory)/Ini.o: Ini.cpp Ini.h
	@$(MakeDirCommand) $(@D)
	@$(CXX) $(CXXFLAGS) $(IncludePath) -c $< -o $(IntermediateDirectory)/Ini.o
	
## libRandom
libRandom.a: $(IntermediateDirectory)/Random.o
	@$(MakeDirCommand) $(@D)
	$(AR) libRandom.a $(IntermediateDirectory)/Random.o
	
$(IntermediateDirectory)/Random.o: Random.cpp Random.h
	@$(MakeDirCommand) $(@D)
	@$(CXX) $(CXXFLAGS) $(IncludePath) -c $< -o $(IntermediateDirectory)/Random.o

## libArguments
libArguments.a: $(IntermediateDirectory)/Arguments.o
	@$(MakeDirCommand) $(@D)
	$(AR) libArguments.a $(IntermediateDirectory)/Arguments.o
	
$(IntermediateDirectory)/Arguments.o: Arguments.cpp Arguments.h
	@$(MakeDirCommand) $(@D)
	@$(CXX) $(CXXFLAGS) $(IncludePath) -c $< -o $(IntermediateDirectory)/Arguments.o

MakeIntermediateDirs:
	@test -d ./Release || $(MakeDirCommand) ./Release
	
	
clean:
	rm *.a ./Release/*.o
