##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=AMS-DEMO
ConfigurationName      :=Debug
WorkspacePath          :=/home/matjaz/Todo/BogdanFilipic/EkgSim
ProjectPath            :=/home/matjaz/Todo/BogdanFilipic/EkgSim/AMS-DEMO
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Matjaž
Date                   :=10/02/17
CodeLitePath           :=/home/matjaz/.codelite
LinkerName             :=/usr/bin/mpic++
SharedObjectLinkerName :=/usr/bin/mpic++ -shared -fPIC
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
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=$(PreprocessorSwitch)NO_MPI 
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="AMS-DEMO.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch)../copyOfLibs 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)rt $(LibrarySwitch)dl $(LibrarySwitch)Ini $(LibrarySwitch)Arguments $(LibrarySwitch)Random 
ArLibs                 :=  "rt" "dl" "Ini" "Arguments" "Random" 
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)../copyOfLibs 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/mpic++
CC       := /usr/bin/mpicc
CXXFLAGS :=  -g -O0 -std=c++11 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/GeneralOptimizationAlgorithm.cpp$(ObjectSuffix) $(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IntermediateDirectory)/utilities.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/GeneralOptimizationAlgorithm.cpp$(ObjectSuffix): GeneralOptimizationAlgorithm.cpp $(IntermediateDirectory)/GeneralOptimizationAlgorithm.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/matjaz/Todo/BogdanFilipic/EkgSim/AMS-DEMO/GeneralOptimizationAlgorithm.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/GeneralOptimizationAlgorithm.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/GeneralOptimizationAlgorithm.cpp$(DependSuffix): GeneralOptimizationAlgorithm.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/GeneralOptimizationAlgorithm.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/GeneralOptimizationAlgorithm.cpp$(DependSuffix) -MM GeneralOptimizationAlgorithm.cpp

$(IntermediateDirectory)/GeneralOptimizationAlgorithm.cpp$(PreprocessSuffix): GeneralOptimizationAlgorithm.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/GeneralOptimizationAlgorithm.cpp$(PreprocessSuffix)GeneralOptimizationAlgorithm.cpp

$(IntermediateDirectory)/main.cpp$(ObjectSuffix): main.cpp $(IntermediateDirectory)/main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/matjaz/Todo/BogdanFilipic/EkgSim/AMS-DEMO/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.cpp$(DependSuffix): main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/main.cpp$(DependSuffix) -MM main.cpp

$(IntermediateDirectory)/main.cpp$(PreprocessSuffix): main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.cpp$(PreprocessSuffix)main.cpp

$(IntermediateDirectory)/utilities.cpp$(ObjectSuffix): utilities.cpp $(IntermediateDirectory)/utilities.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/matjaz/Todo/BogdanFilipic/EkgSim/AMS-DEMO/utilities.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/utilities.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/utilities.cpp$(DependSuffix): utilities.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/utilities.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/utilities.cpp$(DependSuffix) -MM utilities.cpp

$(IntermediateDirectory)/utilities.cpp$(PreprocessSuffix): utilities.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/utilities.cpp$(PreprocessSuffix)utilities.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


