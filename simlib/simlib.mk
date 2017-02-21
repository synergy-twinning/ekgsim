##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Release
ProjectName            :=simlib
ConfigurationName      :=Release
WorkspacePath          :=.
ProjectPath            :=.
IntermediateDirectory  :=./Release
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Matjaž
Date                   :=21/02/17
CodeLitePath           :=/home/matjaz/.codelite
LinkerName             :=/usr/bin/g++
SharedObjectLinkerName :=/usr/bin/g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/lib$(ProjectName).a
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="simlib.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch)../copyOfLibs $(IncludeSwitch)AMD-DEMO 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)copyOfLibs 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++
CC       := /usr/bin/gcc
CXXFLAGS :=  -std=c++11 $(Preprocessors)
CFLAGS   :=   $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/simulator.cpp$(ObjectSuffix) $(IntermediateDirectory)/lib_main.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(IntermediateDirectory) $(OutputFile)

$(OutputFile): $(Objects)
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(AR) $(ArchiveOutputSwitch)$(OutputFile) @$(ObjectsFileList) $(ArLibs)
	@$(MakeDirCommand) "/home/matjaz/Todo/BogdanFilipic/EkgSim/.build-release"
	@echo rebuilt > "/home/matjaz/Todo/BogdanFilipic/EkgSim/.build-release/simlib"

MakeIntermediateDirs:
	@test -d ./Release || $(MakeDirCommand) ./Release


./Release:
	@test -d ./Release || $(MakeDirCommand) ./Release

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/simulator.cpp$(ObjectSuffix): simulator.cpp $(IntermediateDirectory)/simulator.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "$(ProjectPath)/simulator.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/simulator.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/simulator.cpp$(DependSuffix): simulator.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/simulator.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/simulator.cpp$(DependSuffix) -MM simulator.cpp

$(IntermediateDirectory)/simulator.cpp$(PreprocessSuffix): simulator.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/simulator.cpp$(PreprocessSuffix)simulator.cpp

$(IntermediateDirectory)/lib_main.cpp$(ObjectSuffix): lib_main.cpp $(IntermediateDirectory)/lib_main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "$(ProjectPath)/lib_main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/lib_main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/lib_main.cpp$(DependSuffix): lib_main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/lib_main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/lib_main.cpp$(DependSuffix) -MM lib_main.cpp

$(IntermediateDirectory)/lib_main.cpp$(PreprocessSuffix): lib_main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/lib_main.cpp$(PreprocessSuffix)lib_main.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Release/


