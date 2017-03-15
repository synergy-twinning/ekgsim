/*  
    Copyright (c) 2017 Institute Jožef Stefan, Jamova cesta 39, SI-1000, Ljubljana, Slovenija

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    Please cite the following works (bibtex source below):
        - DepolliAvbeljTrobec2008 for the simulator and simulation-based optimization
        - DepolliTrobecFilipic2013 for the AMS-DEMO optimizer
        - TrobecDepolliAvbelj2009 for the simulator

    @article{DepolliAvbeljTrobec2008,
        author = {Depolli, Matjaž and Avbelj, Viktor and Trobec, Roman},
        title = {Computer-Simulated Alternative Modes of {U}-Wave Genesis},
        journal = {Journal of Cardiovascular Electrophysiology},
        volume = {19},
        number = {1},
        publisher = {Blackwell Publishing Inc},
        issn = {1540-8167},
        url = {http://dx.doi.org/10.1111/j.1540-8167.2007.00978.x},
        doi = {10.1111/j.1540-8167.2007.00978.x},
        pages = {84--89},
        keywords = {U wave, ECG, action potential, repolarization, myocardium, computer simulation},
        year = {2008}
    }

    @article{DepolliTrobecFilipic2013,
        author = {Depolli, Matjaž and Trobec, Roman and Filipič, Bogdan},
        title = {Asynchronous master-slave parallelization of differential evolution for multiobjective optimization},
        journal = {Evolutionary Computation},
        volume = {21},
        number = {2},
        pages = {261-291},
        doi = {10.1162/EVCO_a_00076},
        issn = {1063-6560},
        url = {http://www.mitpressjournals.org/doi/abs/10.1162/EVCO_a_00076},
        year = {2013}
    }

    @inproceedings{TrobecDepolliAvbelj2009,
        title           = {Simulation of {ECG} repolarization phase with improved model of cell action potentials},
        author          = {Trobec, Roman and Depolli, Matja{\v{z}} and Avbelj, Viktor},
        booktitle       = {International Joint Conference on Biomedical Engineering Systems and Technologies},
        pages           = {325--332},
        year            = {2009},
        organization    = {Springer}
    }
*/

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
