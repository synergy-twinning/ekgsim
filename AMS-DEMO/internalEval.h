#ifndef INTERNALEVAL_H_INCLUDED
#define INTERNALEVAL_H_INCLUDED


#include <string>
#include <vector>
#include "timer.h"


#if defined WIN32 || defined WIN64
    #include <windows.h>
    #include <process.h>
    typedef void (WINAPI * PPROC) (double *, double *, int, int);
    #define LIBHANDLE HANDLE
    #define GetProcedure(dl, name) GetProcAddress(reinterpret_cast<HINSTANCE>(dl), name)
    #define CloseDynalink(dl) FreeLibrary(reinterpret_cast<HMODULE>(dl))
    #define OpenDynalink(name) LoadLibrary(name);
    #define ErrorDynalink GetLastError
    const char* cec2007LibName = "fsuite.dll";
    #define usleep(usec)
#else
    // linux
    #include <dlfcn.h>
    #include <pthread.h>
    #include <unistd.h>
    typedef void (*PPROC) (double *, double *, int, int);
    #define LIBHANDLE void*
    #define GetProcedure dlsym
    #define CloseDynalink dlclose
    #define OpenDynalink(name) dlopen(name, RTLD_NOW);
    #define ErrorDynalink dlerror
    const char* cec2007LibName = "./fsuite.so";
#endif


struct TestFunctionCEC2007 {
    PPROC function;
    LIBHANDLE libraryHandle;
    int selectedFunction;
    bool loaded;
    int sleepTime;

    TestFunctionCEC2007() : loaded(false), sleepTime(0) {init();}
    ~TestFunctionCEC2007() {if (loaded) CloseDynalink(libraryHandle);}

    void init() {
        libraryHandle = OpenDynalink(cec2007LibName);
        loaded = (libraryHandle);
        if (!loaded)
            throw ErrorDynalink();
    }

    void setSleepTime(int t) {
        sleepTime = t;
    }

    void selectFunction(const std::string& name) {
        selectedFunction = -1;
        static const char* funcNames[] = {"OKA2", "SYMPART", "S_ZDT1", "S_ZDT2", "S_ZDT4", "R_ZDT4", "S_ZDT6", "S_DTLZ2", "R_DTLZ2", "S_DTLZ3", "WFG1", "WFG8", "WFG9"};
        for (size_t i = 0; i < sizeof(funcNames)/sizeof(char*); ++i) {
            if (name == funcNames[i]) {
                selectedFunction = i;
                break;
            }
        }

        if (selectedFunction >= 0)
            function = (PPROC) GetProcedure(libraryHandle, funcNames[selectedFunction]);
        else
            throw "TestFunctionCEC2007 - invalid function selected";
    }

    template<class Vector>
    inline void operator() (const Vector& input, Vector& output, double& violation, Vector& properties) const {
        PrecisionTimer pt;
        pt.start();
        function(const_cast<double*>(&(input[0])), &(output[0]), input.size(), output.size());
        violation = 0;
        /*
        // on new ubuntu distro, usleep has about 75us resolution but only 2000us on FC2
        if (sleepTime > 0)
            usleep(sleepTime);
        */
        pt.waitTotal(sleepTime);
    }
};


#endif // INTERNALEVAL_H_INCLUDED
