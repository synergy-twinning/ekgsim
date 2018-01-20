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

#ifndef HYPERMATRIX_H
#define HYPERMATRIX_H


#include <assert.h>


namespace Pattern {


	// ***********************************************************************************************
	// definitions
	//
	template<class T, size_t N>
	class Vector;


	template<class T, size_t N, class Container>
	class Hypermatrix;


	// ---------------------------------------------------------------------------------------------
	// private namespace
	//
	namespace {

		// ***********************************************************************************************
		// definitions
		//
		template<class T, size_t N, class Container>
		class Hyperplane {
			Container&					data;
			const size_t* const			dimensions;
			size_t							offset;
			
			typedef Hyperplane<T, N-1, Container> Subplane;
			
		public:
			Hyperplane(Container& cont, const size_t* const dimensionVector, size_t offs = 0);
			
			inline Subplane operator[] (size_t nthDimension);
		};
		
		
		template<class T, class Container>
		class Hyperplane<T, 1, Container> {
			Container&					data;
			const size_t* const			dimensions;
			size_t							offset;
			
		public:
			Hyperplane(Container& cont, const size_t* const dimensionVector, size_t offs = 0);
			
			inline T& operator[] (size_t nthDimension);
			//const T& operator[] (size_t nthDimension) const;
		};
		
		
		template<class T, class Container>
		class Hyperplane<T, 1, const Container> {
			const Container&			data;
			const size_t* const			dimensions;
			size_t							offset;
			
		public:
			Hyperplane(const Container& cont, const size_t* const dimensionVector, size_t offs = 0);
			
			const T& operator[] (size_t nthDimension) const;
		};
		
	}
	//
	// end of private namespace
	// ---------------------------------------------------------------------------------------------
	
	
	// ***********************************************************************************************
	// template class ConstInit
	///
	/// used to initialize vectors with constant value (same value for all elements)
	///	
	template<class T>
    struct CInit {
        T val;
        
        CInit(const T t = 1) : val(t) {}
        operator T() const {return val;}
    };
    
    template<class T> CInit<T> cInit(T val) {
        return CInit<T>(val);
    }
    
    
    // ***********************************************************************************************
	// template class ZeroInit
	//	
	///
	/// specialization of ConstInit
	/// used to initialize vectors with 0 (all elements)
	///	
	template<class T>
    struct ZInit : CInit<T> {
        ZInit() : CInit<T>(0) {}
    };
    

	// ***********************************************************************************************
	// template class Vector
	//
	template<class T, size_t N>
	class Vector {
		T									data[N];
		
	public:
		enum {
			dimension = N
		};
		
		typedef CInit<T> ConstInit;
		typedef ZInit<T> ZeroInit;
	
	public:
		// ctors ++++++++++++++++++++++++++++++++++++++++++++
		Vector() {}
		
        explicit Vector(T val) {
            for (size_t i = 0; i < N; ++i)
				data[i] = val;
        }
		
		template<class T1>
		Vector(const CInit<T1> c) {
			for (size_t i = 0; i < N; ++i)
				data[i] = c;
		}
		
		
		template<class T2> 
		inline Vector(const Vector<T2, N>& other) {
			this->operator= (other);
		}
		
		// operators ++++++++++++++++++++++++++++++++++++++++
		inline Vector& operator= (const Vector& other);
		
		inline const T operator[] (size_t position) const;
		inline T& operator[] (size_t position);
		
		inline T* const c_array();
		inline const T* const c_array() const;
		
		template<class T2>
		inline Vector& operator= (const Vector<T2, N>& v2) {
			for (size_t i = 0; i < N; ++i)
				data[i] = (T2)v2[i];
			return *this;
		}
		
		template<class T2>
		inline Vector& operator+= (const Vector<T2, N>& v2) {
			for (size_t i = 0; i < N; ++i)
				data[i] += (T2)v2[i];
			return *this;
		}
		
		template<class T2>
		inline Vector& operator-= (const Vector<T2, N>& v2) {
			for (size_t i = 0; i < N; ++i)
				data[i] -= (T2)v2[i];
			return *this;
		}
		
		template<class T2>
		inline Vector& operator*= (T2 val) {
			for (size_t i = 0; i < N; ++i)
				data[i] *= val;
			return *this;
		}
			
		// other ++++++++++++++++++++++++++++++++++++++++++++
		inline static size_t size() {
			return N;
		}
		
		inline friend T sqrLength(const Vector& v) {
			T 								sum = v.data[0] * v.data[0];
			for (size_t i = 1; i < N; ++i) {
				sum += v.data[i] * v.data[i];
			}
			return sum;
		}
		
		inline friend T pow3Length(const Vector& v) {
			T 								l = sqrLength(v);
			return l * sqrt(l);
		}
		
		/// normalize the input vector (inplace)
		inline friend void normalize(Vector& v) {
			T l = sqrLength(v);
			if (l != 0)
				v *= 1.0/l;
			else 
				v *= 0.0;
		}
		
		/// cross product of two vectors
		inline friend Vector cross(const Vector &a, const Vector &b) {
			// this function should be enabled only for vectors of 3 elements ...
			if (N != 3)
				return Vector((T)(1.0/0.0));
				 
			Vector ab;
			ab[0] = a[1]*b[2] - a[2]*b[1];
			ab[1] = a[0]*b[2] - a[2]*b[0];
			ab[2] = a[0]*b[1] - a[1]*b[0];
			return ab;
		}
		
		inline friend Vector operator* (T scalar, const Vector& v) {
			Vector result = v;
			result*=scalar;
			return result;
		}
	};

	template<class T, size_t N>
	inline Vector<T, N> operator- (const Vector<T, N>& v1, const Vector<T, N>& v2);
	template<class T, size_t N>
	inline Vector<T, N> operator+ (const Vector<T, N>& v1, const Vector<T, N>& v2);
	template<class T, size_t N>	
	inline T operator* (const Vector<T, N>& v1, const Vector<T, N>& v2);
	template<class T, size_t N>	
	inline Vector<T, N> operator* (const Vector<T, N>& v1, T scalar);
	template<class T, class S, size_t N>	
	inline Vector<S, N> operator* (const Vector<T, N>& v1, S scalar);
	// dominatates / is-dominated-by operators (non strict)
	template<class T, size_t N>	
	inline bool operator>= (const Vector<T, N>& v1, const Vector<T, N>& v2);
	template<class T, size_t N>	
	inline bool operator<= (const Vector<T, N>& v1, const Vector<T, N>& v2);
	template<class T, size_t N>	
	inline Vector<T, N> abs(const Vector<T, N>& other);


	// ***********************************************************************************************
	// template class Region
	//
	/// Defines a box in space with boundaries parallel to space axes. This effectivelly means
	/// every coordinate in space gets its minimum and maximum value.
	///
	template<class Vec>
	class Region {
		Vec									minPoint, maxPoint;
		
	public:
		Region() {}
		Region(const Vec& vmin, const Vec& vmax) : minPoint(vmin), maxPoint(vmax) {}
		
		void set(const Vec& vmin, const Vec& vmax) {
			minPoint = (vmin);
			maxPoint = (vmax);
		}
		
		bool inside(const Vec& v) {
			for (size_t i = 0; i < Vec::size(); ++i) {
				if ((v[i] < minPoint[i]) || (v[i] >= maxPoint[i]))
					return false;
			}
			return true;
		}
		
		bool incIndex(Vec& index) {
			for (size_t dim = Vec::size(); dim > 0; --dim) {
				if (++index[dim-1] >= maxPoint[dim-1]) {
					index[dim-1] = minPoint[dim-1];
				} else {
					return true;
				}
			}
			return false;
		}
	};


	// ***********************************************************************************************
	// template class Hypermatrix
	//
	template<class T, size_t N, class Container>
	class Hypermatrix {
		typedef Hyperplane<T, N-1, Container> HyperplaneType;
		typedef Hyperplane<T, N-1, const Container> ConstHyperplaneType;
		
	public:
		typedef Vector<size_t, N> SizeVector;
		typedef T DataType;
		enum {Dimensionality = N};
				
		struct ConstInit {
			T val;
			
			explicit ConstInit(const T t = 1) : val(t) {}
			operator T() const {return val;}
		};
		
		struct ZeroInit : ConstInit {
			ZeroInit() : ConstInit(0) {}
		};
		
	protected:
		Container						data;
		SizeVector						dimensions;
		
	public:
		Hypermatrix() {}
		Hypermatrix& operator= (const Hypermatrix& other);
		
		bool resize(const SizeVector& newDimensions);
		bool resize(const SizeVector& newDimensions, const ConstInit& c);
		const SizeVector& size() const;
		
		inline HyperplaneType operator[] (size_t nthDimension);
		inline ConstHyperplaneType operator[] (size_t nthDimension) const;
		inline T& operator[] (const SizeVector& position);
		inline const T& operator[] (const SizeVector& position) const;
		
		inline size_t positionToIndex(const SizeVector& position) const;
		T& elementAtIndex(size_t i);
		const T& elementAtIndex(size_t i) const;
		
		size_t dimensionality() const {return N;}
	};
	
	
	template<class MatSrc, class MatDest>
	bool inject(const MatSrc& src, MatDest& dest, const typename MatDest::SizeVector& offs = MatDest::SizeVector(typename MatDest::SizeVector::ZeroInit())) {
		if ((src.dimensionality() == dest.dimensionality()) && (offs + src.size() <= dest.size())) {
			typedef typename MatSrc::SizeVector Vec;
			typedef typename Vec::ZeroInit VecZeroInit;
			VecZeroInit init0;
			Vec index(init0);
			Region<Vec> srcRegion(index, src.size());
			bool indexOk = true;
			while(indexOk) {
				dest[index + offs] = src[index];
				indexOk = srcRegion.incIndex(index);
			};
		}
		return true;
	}

}


// ------------------------------------------------------------------------------------------------
// implementation
// ------------------------------------------------------------------------------------------------


namespace Pattern {

	namespace {
		
		template<class T, size_t N, class Container>
		Hyperplane<T, N, Container>::Hyperplane(Container& cont, const size_t* const dimensionVector, size_t offs) 
			: data(cont),
			dimensions(dimensionVector + 1),
			offset(offs)
		{
		}
		
		
		template<class T, size_t N, class Container>
		typename Hyperplane<T, N, Container>::Subplane Hyperplane<T, N, Container>::operator[] (size_t nthDimension) {
			assert(nthDimension < dimensions[0]);
			return Subplane(data, dimensions, nthDimension + dimensions[0] * offset);
		}
		
		
		template<class T, class Container>
		Hyperplane<T, 1, Container>::Hyperplane(Container& cont, const size_t* const dimensionVector, size_t offs) 
			: data(cont),
			dimensions(dimensionVector + 1),
			offset(offs)
		{
		}
		
		
		template<class T, class Container>
		T& Hyperplane<T, 1, Container>::operator[] (size_t nthDimension) {
			assert(nthDimension < dimensions[0]);
			return data[nthDimension + dimensions[0] * offset];
		}
		
		
		template<class T, class Container>
		Hyperplane<T, 1, const Container>::Hyperplane(const Container& cont, const size_t* const dimensionVector, size_t offs) 
			: data(cont),
			dimensions(dimensionVector + 1),
			offset(offs)
		{
		}
		
		
		template<class T, class Container>
		const T& Hyperplane<T, 1, const Container>::operator[] (size_t nthDimension) const {
			assert(nthDimension < dimensions[0]);
			return data[nthDimension + dimensions[0] * offset];
		}
	}

	template<class T> 
	T abs(const T var) {
		return (var > 0 ? var : -var);
	}
	
	/// ***********************************************************************************************
	/// template class Vector
	///
	template<class T, size_t N>
	Vector<T, N>& Vector<T, N>::operator= (const Vector<T, N>& other) {
		for (size_t i = 0; i < N; ++i)
			data[i] = other.data[i];
		return *this;
	}
	
	
	template<class T, size_t N>
	const T Vector<T, N>::operator[] (size_t position) const {
		return data[position];
	}


	template<class T, size_t N>
	T& Vector<T, N>::operator[] (size_t position) {
		return data[position];
	}


	template<class T, size_t N>
	T* const Vector<T, N>::c_array() {
		return data;
	}


	template<class T, size_t N>
	const T* const Vector<T, N>::c_array() const {
		return data;
	}
	
	
	template<class T, size_t N>
	Vector<T, N> operator- (const Vector<T, N>& v1, const Vector<T, N>& v2) {
		Vector<T, N>					result;
		for (size_t i = 0; i < N; ++i)
			result[i] = v1[i] - v2[i];
		return result;
	}
	
	
	template<class T, size_t N>
	Vector<T, N> operator+ (const Vector<T, N>& v1, const Vector<T, N>& v2) {
		Vector<T, N>					result;
		for (size_t i = 0; i < N; ++i)
			result[i] = v1[i] + v2[i];
		return result;
	}
	
	
	template<class T, size_t N>	
	T operator* (const Vector<T, N>& v1, const Vector<T, N>& v2) {
		T									result(0);
		for (size_t i = 0; i < N; ++i)
			result += v1[i] * v2[i];
		return result;
	}
	
	
	template<class T, size_t N>	
	Vector<T, N> operator* (const Vector<T, N>& v1, T scalar) {
		Vector<T, N>					result;
		for (size_t i = 0; i < N; ++i)
			result[i] = v1[i] * scalar;
		return result;
	}

	
	template<class T, class S, size_t N>	
	Vector<S, N> operator* (const Vector<T, N>& v1, S scalar) {
		Vector<S, N>					result;
		for (size_t i = 0; i < N; ++i)
			result[i] = v1[i] * scalar;
		return result;
	}
	
	
	template<class T, size_t N>	
	inline bool operator>= (const Vector<T, N>& v1, const Vector<T, N>& v2) {
		bool dominates = true;
		for (size_t i = 0; i < N; ++i) {
			dominates &= (v1[i] >= v2[i]);
		}
		return dominates;
	}
	
	
	template<class T, size_t N>	
	inline bool operator<= (const Vector<T, N>& v1, const Vector<T, N>& v2) {
		return v2 >= v1;
	}
	
	
	template<class T, size_t N>
	Vector<T, N> abs(const Vector<T, N>& other) {
		Vector<T, N>					result;
		for (size_t i = 0; i < N; ++i)
			result[i] = abs(other[i]);
		return result;
	}
	

	/// ***********************************************************************************************
	/// template class Hypermatrix
	///
	template<class T, size_t N, class Container>
	bool Hypermatrix<T, N, Container>::resize(const Vector<size_t, N>& newDimensions) {
		dimensions = newDimensions;
		size_t								newSize = 1;
		for (size_t i = 0; i < N; ++i) 
			newSize *= newDimensions[i];
		data.resize(newSize);
		return true;
	}
	
	
	template<class T, size_t N, class Container>
	bool Hypermatrix<T, N, Container>::resize(const Vector<size_t, N>& newDimensions, const ConstInit& c) {
		dimensions = newDimensions;
		size_t								newSize = 1;
		for (size_t i = 0; i < N; ++i) 
			newSize *= newDimensions[i];
		data.resize(newSize, (T)c);
		return true;
	}
	
	
	template<class T, size_t N, class Container>
	Hypermatrix<T, N, Container>& Hypermatrix<T, N, Container>::operator= (const Hypermatrix<T, N, Container>& other) {
		dimensions = other.dimensions;
		data = other.data;
		return *this;
	}
	
	
	template<class T, size_t N, class Container>
	const typename Hypermatrix<T, N, Container>::SizeVector& Hypermatrix<T, N, Container>::size() const {
		return dimensions;
	}


	template<class T, size_t N, class Container>
	typename Hypermatrix<T, N, Container>::HyperplaneType Hypermatrix<T, N, Container>::operator[] (size_t nthDimension) {
		assert(nthDimension < dimensions[0]);
		return Hyperplane<T, N-1, Container>(data, dimensions.c_array(), nthDimension);
	}
	
	
	template<class T, size_t N, class Container>
	typename Hypermatrix<T, N, Container>::ConstHyperplaneType Hypermatrix<T, N, Container>::operator[] (size_t nthDimension) const {
		assert(nthDimension < dimensions[0]);
		return ConstHyperplaneType(data, dimensions.c_array(), nthDimension);
	}
	
	
	template<class T, size_t N, class Container>
	const T& Hypermatrix<T, N, Container>::operator[] (const typename Hypermatrix<T, N, Container>::SizeVector& position) const {
		return data[positionToIndex(position)];
	}
	
	
	template<class T, size_t N, class Container>
	T& Hypermatrix<T, N, Container>::operator[] (const typename Hypermatrix<T, N, Container>::SizeVector& position) {
		return data[positionToIndex(position)];
	}
	
	
	template<class T, size_t N, class Container>
	inline size_t Hypermatrix<T, N, Container>::positionToIndex(const typename Hypermatrix<T, N, Container>::SizeVector& position) const {
		size_t offset = position[0];
		assert(position[0] < dimensions[0]);
		for (size_t i = 1; i < N; ++i) {
			assert(position[i] < dimensions[i]);
			offset = offset * dimensions[i] + position[i];
		}
		return offset;
	}
	
	
	template<class T, size_t N, class Container>
	T& Hypermatrix<T, N, Container>::elementAtIndex(size_t i) {
		return data[i];
	}
	
	
	template<class T, size_t N, class Container>
	const T& Hypermatrix<T, N, Container>::elementAtIndex(size_t i) const {
		return data[i];
	}
	
}


#endif //HYPERMATRIX_H
