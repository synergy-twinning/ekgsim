##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=ekgSim
ConfigurationName      :=Debug
WorkspacePath          :=/home/matjaz/Todo/BogdanFilipic/EkgSim
ProjectPath            :=/home/matjaz/Todo/BogdanFilipic/EkgSim
IntermediateDirectory  :=Debug
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
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="ekgSim.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch)copyOfLibs $(IncludeSwitch)/home/matjaz/include $(IncludeSwitch)simlib 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)Ini $(LibrarySwitch)Arguments $(LibrarySwitch)Random $(LibrarySwitch)dl $(LibrarySwitch)simlib 
ArLibs                 :=  "Ini" "Arguments" "Random" "dl" "simlib" 
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)copyOfLibs $(LibraryPathSwitch)simlib/Debug 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/mpic++
CC       := /usr/bin/mpicc
CXXFLAGS :=  -g -O0 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IntermediateDirectory)/sim.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d ".build-debug/simlib" $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

".build-debug/simlib":
	@$(MakeDirCommand) ".build-debug"
	@echo stam > ".build-debug/simlib"




MakeIntermediateDirs:
	@test -d Debug || $(MakeDirCommand) Debug


$(IntermediateDirectory)/.d:
	@test -d Debug || $(MakeDirCommand) Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/main.cpp$(ObjectSuffix): main.cpp $(IntermediateDirectory)/main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/matjaz/Todo/BogdanFilipic/EkgSim/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.cpp$(DependSuffix): main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/main.cpp$(DependSuffix) -MM main.cpp

$(IntermediateDirectory)/main.cpp$(PreprocessSuffix): main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.cpp$(PreprocessSuffix)main.cpp

$(IntermediateDirectory)/sim.cpp$(ObjectSuffix): sim.cpp $(IntermediateDirectory)/sim.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/matjaz/Todo/BogdanFilipic/EkgSim/sim.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/sim.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/sim.cpp$(DependSuffix): sim.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/sim.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/sim.cpp$(DependSuffix) -MM sim.cpp

$(IntermediateDirectory)/sim.cpp$(PreprocessSuffix): sim.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/sim.cpp$(PreprocessSuffix)sim.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r Debug/


