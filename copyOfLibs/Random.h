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
