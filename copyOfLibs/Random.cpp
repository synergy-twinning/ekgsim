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
