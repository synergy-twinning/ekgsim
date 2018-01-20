##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Release
ProjectName            :=ekgSim
ConfigurationName      :=Release
WorkspacePath          :=.
ProjectPath            :=.
IntermediateDirectory  :=./Release
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Matjaž
Date                   :=19/01/18
CodeLitePath           :=/home/matjaz/.codelite
LinkerName             :=mpic++
SharedObjectLinkerName :=mpic++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.o.i
DebugSwitch            :=-gstab
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=$(PreprocessorSwitch)NDEBUG 
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E 
ObjectsFileList        :="ekgSim.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch)copyOfLibs $(IncludeSwitch)simlib 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)Ini $(LibrarySwitch)Arguments $(LibrarySwitch)Random $(LibrarySwitch)dl $(LibrarySwitch)simlib 
ArLibs                 :=  "Ini" "Arguments" "Random" "dl" "simlib" 
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)copyOfLibs $(LibraryPathSwitch)simlib/Release 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := ar rcus
CXX      := mpic++
CC       := mpicc
CXXFLAGS :=  -O2 -std=c++11 -Wall -Wno-unused-but-set-variable -Wno-deprecated-declarations -Wno-unused-function $(Preprocessors)
CFLAGS   :=  -O2 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/copyOfLibs_Ini.cpp$(ObjectSuffix) $(IntermediateDirectory)/copyOfLibs_Random.cpp$(ObjectSuffix) $(IntermediateDirectory)/copyOfLibs_Arguments.cpp$(ObjectSuffix) $(IntermediateDirectory)/sim.cpp$(ObjectSuffix) $(IntermediateDirectory)/main.cpp$(ObjectSuffix) 



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
	@test -d ./Release || $(MakeDirCommand) ./Release


$(IntermediateDirectory)/.d:
	@test -d ./Release || $(MakeDirCommand) ./Release

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/copyOfLibs_Ini.cpp$(ObjectSuffix): copyOfLibs/Ini.cpp $(IntermediateDirectory)/copyOfLibs_Ini.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "copyOfLibs/Ini.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/copyOfLibs_Ini.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/copyOfLibs_Ini.cpp$(DependSuffix): copyOfLibs/Ini.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/copyOfLibs_Ini.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/copyOfLibs_Ini.cpp$(DependSuffix) -MM copyOfLibs/Ini.cpp

$(IntermediateDirectory)/copyOfLibs_Ini.cpp$(PreprocessSuffix): copyOfLibs/Ini.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/copyOfLibs_Ini.cpp$(PreprocessSuffix) copyOfLibs/Ini.cpp

$(IntermediateDirectory)/copyOfLibs_Random.cpp$(ObjectSuffix): copyOfLibs/Random.cpp $(IntermediateDirectory)/copyOfLibs_Random.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "copyOfLibs/Random.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/copyOfLibs_Random.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/copyOfLibs_Random.cpp$(DependSuffix): copyOfLibs/Random.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/copyOfLibs_Random.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/copyOfLibs_Random.cpp$(DependSuffix) -MM copyOfLibs/Random.cpp

$(IntermediateDirectory)/copyOfLibs_Random.cpp$(PreprocessSuffix): copyOfLibs/Random.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/copyOfLibs_Random.cpp$(PreprocessSuffix) copyOfLibs/Random.cpp

$(IntermediateDirectory)/copyOfLibs_Arguments.cpp$(ObjectSuffix): copyOfLibs/Arguments.cpp $(IntermediateDirectory)/copyOfLibs_Arguments.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "copyOfLibs/Arguments.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/copyOfLibs_Arguments.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/copyOfLibs_Arguments.cpp$(DependSuffix): copyOfLibs/Arguments.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/copyOfLibs_Arguments.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/copyOfLibs_Arguments.cpp$(DependSuffix) -MM copyOfLibs/Arguments.cpp

$(IntermediateDirectory)/copyOfLibs_Arguments.cpp$(PreprocessSuffix): copyOfLibs/Arguments.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/copyOfLibs_Arguments.cpp$(PreprocessSuffix) copyOfLibs/Arguments.cpp

$(IntermediateDirectory)/sim.cpp$(ObjectSuffix): sim.cpp $(IntermediateDirectory)/sim.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "sim.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/sim.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/sim.cpp$(DependSuffix): sim.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/sim.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/sim.cpp$(DependSuffix) -MM sim.cpp

$(IntermediateDirectory)/sim.cpp$(PreprocessSuffix): sim.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/sim.cpp$(PreprocessSuffix) sim.cpp

$(IntermediateDirectory)/main.cpp$(ObjectSuffix): main.cpp $(IntermediateDirectory)/main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.cpp$(DependSuffix): main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/main.cpp$(DependSuffix) -MM main.cpp

$(IntermediateDirectory)/main.cpp$(PreprocessSuffix): main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.cpp$(PreprocessSuffix) main.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Release/


