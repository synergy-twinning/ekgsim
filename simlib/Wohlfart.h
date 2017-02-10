#ifndef WOHLFART_H_INCLUDED
#define WOHLFART_H_INCLUDED

/**
    AP model by Wohlfart
**/

#include <algorithm>
#include <cmath>
#include <vector>


/**
	Original Wohlfart formula + 8th parameter that represents an offset in y
**/
class Wohlfart {
	double k[8];
	
public:
	Wohlfart() {}
	Wohlfart(double newK[8]) {
		setK(newK);
	}
	
	void setK(double newK[8]) {
		using namespace std;
		copy(newK, newK+8, k);
	}
	
	void setK(const std::vector<double>& newK) {
		using namespace std;
		copy(newK.begin(), newK.end(), k);
	}
	
	const double* getK() const {
		return k;
	}
	
	double* getK() {
		return k;
	}
	
	inline double operator[](double t) const {
		using namespace std;
		return  
			(k[1] * ((1.0 - k[2]) * exp(-k[3] * t) + k[2])) * 
			(exp(-k[4] * t) / (1.0 + exp(k[5] * (t - k[6])))) /
			(1.0 + exp(-k[0] * t)) +
			k[7];
	}
	
	double apd90() const {
		const size_t d_size = 1000;
		double d[d_size];
		
		for (size_t di = 0; di < d_size; ++di) {
			d[di] = (*this)[di];
		}
		
		double target = *std::max_element(d, d+d_size) * 0.1 + k[7] * 0.9;
		double time = -1.0;
		for (size_t di = 1; di < d_size; ++di) {
			if ((d[di-1] > target) && (d[di] < target)) {
				time = interpolate(double(di-1), (d[di-1] - target), 
					double(di), (target - d[di]));
			}
		}
		return time;
	}
	
private:
	template<class T>
	static double interpolate(T a, T wa, T b, T wb) {
		return (a * wa + b * wb) / (wa + wb);
	}
};

/*
void saveWohlfart(const Wohlfart& w, const char* fname = "AP.txt") {
	ofstream file(fname);
	file << "# saveWohlfart";
	for (int i = 0; i < 8; ++i) 
		file << " " << w.getK()[i];
	file << "\n";
	for (double d = 0; d < 700; d += 0.10)
		file << d << "\t" << w[d] << "\n";
}*/


/**
	Modified Wohlfart formula:
		- original 7th parameter is now 8th parameter (1-based count)
		- original 6th parameter is broken into two new paramters: 6th and 7th
			6th and 7th parameter now separately control phase 3 slope in its
			upper and lower part, respectively
		- added 9th parameter that represents an offset
	
	matlab formula: 
		apFunc = k(9) + ((1 + exp(-k(1) * t)) .^ -1) .* (k(2) * ((1 - k(3)) * exp(-k(4)*t) + k(3))) .* (exp(-k(5) * t) .* (1 - (1 + exp(-k(7) * (t - k(8)) + log(2^(k(7)/k(6)) - 1))).^-(k(6)/k(7))));
		
    to get value in certain time, call it as v = wohlfartPlus[t];
**/
class WohlfartPlus {
public:
	static const size_t numParams = 9;

protected:
	double k[numParams];
	
public:
	WohlfartPlus() {}
	WohlfartPlus(const double newK[numParams]) {
		setK(newK);
	}
	
	void setK(const double newK[numParams]) {
		using namespace std;
		copy(newK, newK+numParams, k);
	}
	
	void setK(const std::vector<double>& newK) {
		using namespace std;
		copy(newK.begin(), newK.end(), k);
	}
	
	inline const double* getK() const {return k;}
	
	inline double* getK() {return k;}
	
	/// returns value of the AP in time t (t >= 0)
	inline double operator[](double t) const {
		using namespace std;
		return  
			(1.0 / (1.0 + exp(-k[1] * t))) 
			* (k[2] * ((1.0 - k[3]) * exp(-k[4] * t) + k[3])) 
			* (exp(-k[5] * t) * (1 - pow((1 + exp(-k[7] * (t - k[8]) 
				+ log(pow(2, (k[7] / k[6])) - 1))), -(k[6] / k[7]))))
			+ k[0];
	}
	
	/// APD90 is the time in which AP repolarizes back to 90% of the normal value
	double apd90() const {
		const size_t d_size = 1000;
		double d[d_size];
		
		for (size_t di = 0; di < d_size; ++di) {
			d[di] = (*this)[di];
		}
		
		double target = k[0] + k[2] * 0.1;
		double time = -1.0;
		for (size_t di = 1; di < d_size; ++di) {
			if ((d[di-1] > target) && (d[di] <= target)) {
				time = interpolate(double(di-1), (d[di-1] - target), 
					double(di), (target - d[di]));
			}
		}
		return time;
	}
	
private:
	template<class T>
	static inline double interpolate(T a, T wa, T b, T wb) {
		return (a * wa + b * wb) / (wa + wb);
	}
};

#endif // WOHLFART_H_INCLUDED
