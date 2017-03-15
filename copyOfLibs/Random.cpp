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

#include "Random.h"
#include <ctime>


namespace Random {
	
	/// class SingleRandom
	void CRand::setSeed(int seed) {
		srand(seed);
		rand();
	}
	
	
	unsigned int CRand::randomizeSeed() {
		unsigned int timeSeed = time(NULL);
		srand(timeSeed);
		// throw away first generated number
		rand();
		return timeSeed;
	}

    //#define rotl(r,n) (((r)<<(n)) | ((r)>>((8*sizeof(r))-(n))))
    template<class R, class N>
    R rotl(R r, N n) {
        return (((r)<<(n)) | ((r)>>((8*sizeof(r))-(n))));
    }
    
    void RersResrResdra::setSeed(uint32_t seed) {
        xx = 914489ULL; 
        yy = 8675416ULL; 
        zz = 439754684ULL;
        for (unsigned int n=((seed>>22)&0x3ff)+20; n>0; n--) xx = rotl(xx,52) - rotl(xx, 9);
        for (unsigned int n=((seed>>11)&0x7ff)+20; n>0; n--) yy = rotl(yy,24) - rotl(yy,45);
        for (unsigned int n=((seed    )&0x7ff)+20; n>0; n--) zz -= rotl(zz,38);
    }
    
    uint64_t RersResrResdra::newValue() {  // Combined period = 2^116.23
       xx = rotl(xx,8) - rotl(xx,29);                 //RERS,   period = 4758085248529 (prime)
       yy = rotl(yy,21) - yy;  yy = rotl(yy,20);      //RESR,   period = 3841428396121 (prime)
       zz = rotl(zz,42) - zz;  zz = zz + rotl(zz,14); //RESDRA, period = 5345004409 (prime)
       return xx ^ yy ^ zz;
    }
    
    void ThreeResr::setSeed(uint32_t seed) {
        xx = 590009ULL;
        yy = 8675416ULL;
        zz = 46017471ULL;
        for (unsigned int n=((seed>>22)&0x3ff)+20; n>0; n--) { xx = rotl(xx,43) - xx; xx = rotl(xx,27); }
        for (unsigned int n=((seed>>11)&0x7ff)+20; n>0; n--) { yy = rotl(yy,21) - yy; yy = rotl(yy,20); }
        for (unsigned int n=((seed    )&0x7ff)+20; n>0; n--) { zz = rotl(zz,51) - zz; zz = rotl(zz,26); }
    }
    
    uint64_t ThreeResr::newValue() {  // Combined period = 2^123.32
        xx = rotl(xx,43) - xx; xx = rotl(xx,27); // RESR, period =  9925159703554 = 2*53*93633582109
        yy = rotl(yy,21) - yy; yy = rotl(yy,20); // RESR, period =  3841428396121 (prime)
        zz = rotl(zz,51) - zz; zz = rotl(zz,26); // RESR, period =  348142888313 = 11*11*2877213953
        return xx ^ yy ^ zz;
    }

}
