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

#ifndef RANDOM_H_INCLUDED
#define RANDOM_H_INCLUDED


#include <cstdlib>
#include <cstdint>


namespace Random {

    /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /// wrapper RandBaseInt
    ///
    /// implementation should provide:
    ///  - static const Int rndMax;
    ///  - implementation ctor should set this->rnd with random generated int between 0 and
    ///     rndMax (both inclusive)
    ///  - int nextRand() const, which should return another random number
    ///
    /// wrapping classes through curiously recurring template pattern (CRTP)
    ///
    template <class RandIntImplementation, class Int = int>
    class RandBaseInt {
    protected:
        Int rnd;

    public:
		operator Int() const 		{
		    return rnd;
        }

		inline double interval(double upperBound) const {
            return rnd * upperBound / double(RandIntImplementation::rndMax);
        }

        inline double interval(double lowerBound, double upperBound) const {
            return lowerBound + interval(upperBound - lowerBound);
        }

        inline double interval(double bounds[2]) const {
            return interval(bounds[0], bounds[1]);
        }

        // interval whose upper bound is exclusive (cannot be generated)
        inline double exclusiveInterval(double upperBound) const {
            return rnd * upperBound / (double(RandIntImplementation::rndMax)+1.0);
        }

        inline double exclusiveInterval(double lowerBound, double upperBound) const {
            return lowerBound + exclusiveInterval(upperBound - lowerBound);
        }

        inline double exclusiveInterval(double bounds[2]) const {
            return exclusiveInterval(bounds[0], bounds[1]);
        }

        inline size_t unsignedInt(size_t numOfPossibilities) const {
            Int maxr = RandIntImplementation::rndMax - (RandIntImplementation::rndMax % numOfPossibilities);
            Int r = rnd;
            while (r > maxr)
                r = static_cast<const RandIntImplementation*>(this)->nextRand();
            return r % numOfPossibilities;
        }

//        inline double gaussian() {
//            return sqrt(-2 * log(newDouble())) * cos(PIx2 * newDouble());
//        }
    };


    /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /// wrapper RandBaseFloat
    ///
    /// implementation should provide:
    ///  - implementation ctor should set this->rnd with random generated int between 0 and 1
    ///     (both inclusive)
    ///  - float nextRand() const, which should generate another random number
    ///
    /// wrapping classes through curiously recurring template pattern (CRTP)
    ///
    template <class RandFloatImplementation, class Float = float>
    class RandBaseFloat {
    protected:
        Float rnd;

    public:
        operator Float() const 		{
            return rnd;
        }

        inline static Float one() {
            static Float _one(1);
            return _one;
        }

		inline double interval(double upperBound) const {
            return rnd * upperBound;
        }

        inline double interval(double lowerBound, double upperBound) const {
            return lowerBound + interval(upperBound - lowerBound);
        }

        inline double interval(double bounds[2]) const {
            return interval(bounds[0], bounds[1]);
        }

        // interval whose upper bound is exclusive (cannot be generated)
        inline double exclusiveInterval(double upperBound) const {
            return interval(upperBound*Float(0.99999999));
        }

        inline double exclusiveInterval(double lowerBound, double upperBound) const {
            return lowerBound + exclusiveInterval(upperBound - lowerBound);
        }

        inline double exclusiveInterval(double bounds[2]) const {
            return exclusiveInterval(bounds[0], bounds[1]);
        }

        inline size_t exclusiveInterval(size_t numOfPossibilities) const {
            return (size_t)floor((numOfPossibilities-1) * rnd + Float(0.5));
        }
    };


	/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /// class CRand
    ///
    /// implementation of random based on generation of integer randoms provided by c std library
    ///
    class CRand : public RandBaseInt<CRand> {
    protected:
	    static const int rndMax = RAND_MAX;
	    friend class Random::RandBaseInt<CRand>;

	public:
		CRand() 					{this->rnd = rand();}
		CRand(int forceRndInt) 		{this->rnd = forceRndInt;}

		int nextRand() const        {return rand();}

		static void setSeed(int seed);
		static unsigned int randomizeSeed();
	};
	
	
	class RersResrResdra : public RandBaseInt<RersResrResdra> {
        /// source: http://drdobbs.com/go-parallel/article/229625477?pgno=2
        /// — Mark Overton is a lead Software Engineer at Northrop Grumman and has multiple patents in image processing. He        can be contacted at MarkDOverton@cox.net. 
        /// Combined period = 2^116.23
        uint64_t xx, yy, zz;
 
    public:
        RersResrResdra() {}
        void setSeed(uint32_t seed);
    protected:
        uint64_t newValue();
    };
 
    class ThreeResr : public RandBaseInt<ThreeResr> {
        /// source: http://drdobbs.com/go-parallel/article/229625477?pgno=2
        /// — Mark Overton is a lead Software Engineer at Northrop Grumman and has multiple patents in image processing. He        can be contacted at MarkDOverton@cox.net. 
        /// Combined period = 2^123.32
        uint64_t xx, yy, zz;
        uint64_t seed;
 
    public:
        ThreeResr() {}
        void setSeed(uint32_t seed);
 
    protected:
        uint64_t newValue();
    };


}


#endif // RANDOM_H_INCLUDED
