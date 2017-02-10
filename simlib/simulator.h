#ifndef SIMULATOR_H_INCLUDED
#define SIMULATOR_H_INCLUDED

/**
    routines for ECG simulation
**/


#include "matrix.h"
#include <stdexcept>


namespace SimLib {
	
	// ***********************************************************************************************
	// struct Settings
	//
	/// All the simulator settings that are read from ini file.
	///
	struct Settings {
		// [input]
		string                       	 	inputShapeFilename;
		string                        		inputExcitationSequenceFilename;
		string								inputPointsFilename;
		string								inputTransferFilename;
		// [APs]
		string								inputActionPotentials;
		string								inputActionPotentialsFilter;
		double								inputActionPotentialsFilterParam;
		double								inputActionPotentialsTimeStep;
		string								inputActionPotentialsFunction;
		double								inputActionPotentialsScale;
		// [output]
		string								outputExcitationSequence;
		string								outputApFilename;
		string								outputFilename;
		// [simulation]
		int 								simulationLength;
		int									simulationStart;
		string								neighbourhoodType;
		double								simulationTimeStep;

		void loadFromIni(const Ini::File& ini) {
			int								modelSec = ini.getSectionNumber("model");
			ini.loadVar(inputShapeFilename, "shape", modelSec);
			ini.loadVar(inputExcitationSequenceFilename, "excitation sequence", modelSec);
			ini.loadVar(inputPointsFilename, "points", modelSec);
			ini.loadVar(inputTransferFilename, "transfer", modelSec);

			int								apSec = ini.getSectionNumber("action potentials");
			ini.loadVar(inputActionPotentials, "input file", apSec);
			ini.loadVar(inputActionPotentialsFilter, "filter", apSec);
			ini.loadVar(inputActionPotentialsTimeStep, "resample time step", apSec);
			ini.loadVar(inputActionPotentialsFilterParam, "filter parameter", apSec);
			ini.loadVar(inputActionPotentialsFunction, "combination function", apSec);
			ini.loadVar(inputActionPotentialsScale, "expected length", apSec);

			int								outputSec = ini.getSectionNumber("output files");
			ini.loadVar(outputExcitationSequence, "excitation sequence", outputSec);
			ini.loadVar(outputApFilename, "action potentials filename", outputSec);
			outputFilename = "result.column";
			ini.loadVar(outputFilename, "results filename", outputSec);

			int								simulationSec = ini.getSectionNumber("simulation");
			ini.loadVar(simulationLength, "length", simulationSec);
			ini.loadVar(simulationStart, "start", simulationSec);
			ini.loadVar(neighbourhoodType, "neighbourhood type", simulationSec);
			simulationTimeStep = 1.0 / 3.0;
			ini.loadVar(simulationTimeStep, "time step", simulationSec);
		}

		Settings() :
            // resonable defaults for some of the settings
			inputExcitationSequenceFilename(""),
			inputActionPotentialsTimeStep(1),
			inputActionPotentialsScale(0),
			simulationLength(800),
			simulationStart(0),
			simulationTimeStep(10)
		{
		}
	};


    /// function for saving vector data as a file (*.txt with tab separated values or *.column with newline separated values, readable by matlab)
    /// timeStep (makes a difference only in *.column files) tells function weather timestamps should be included in the file (first element will have timestamp 0, every next element will have timestamp +timeStep) or not (timeStep < 0)
	template<class Vec>
	void exportVector(const Vec& vec, const string& filename, double timeStep = -1.0f) {
		string								ext = filenameExtension(filename);
		
		if (ext == ".txt") {
			// tab separated values (single line, ignoring time step)
			ofstream							outFile(filename.c_str());			
			
			if (vec.size() > 0)
				outFile << vec[0];
			for (size_t i = 0; i < vec.size(); ++i) 
				outFile << "\t" << vec[i];
		} else if (ext == ".column") {
			// new-line separated values + one comment line + time values
			ofstream							outFile(filename.c_str());
			
			outFile << "# Comment: -none-\n";
			if (timeStep > 0)
				outFile << "# Time[ms]\t";
			outFile << "# Value\n";
			
			double							time = 0;
			for (size_t i = 0; i < vec.size(); ++i) {
				if (i > 0)
					outFile << "\n";
				if (timeStep > 0) {
					outFile << std::fixed << std::setprecision(5) << time << '\t';
					time += timeStep;
				}
				outFile << std::scientific << std::setprecision(6) << vec[i];
			}
		} else {
			throw std::runtime_error("saving vectors to .txt & .column files only for the time being");
		}
	}

	
	/// export multiple vectors (overloaded parameters)
	template<class Vec>
	void exportVectors(const Vec& vec, const string& filename, double startTime, const std::string& comment = "") {
		exportVectors(vec, filename, startTime, -1.0f, comment);
	}
	
	/// export multiple vectors to *.column files. 
	/// vec should be a vector (indexable by [] and having method size_t size()) of elements, with each element having members name, comment, and method s size() and [] for indexing elements 
	template<class Vec>
	void exportVectors(const Vec& vec, const string& filename, double startTime, double timeStep = -1.0f, const std::string& comment = "") {
		string							ext = filenameExtension(filename);
		
		if (ext == ".column") {
			// new-line separated values + one comment line + time values
			ofstream						outFile(filename.c_str());
			
			outFile << "#Comment: " << comment << '\n';
			for (size_t i = 0; i < vec.size(); ++i) {
				if (vec[i].comment != "")
					outFile << "# " << vec[i].name << " : " << vec[i].comment << '\n';
			}
			outFile << "# \n";
			
			if (timeStep > 0)
				outFile << "#Time[ms]\t";
			else 
				outFile << "# ";
			for (size_t i = 0; i < vec.size(); ++i) {
				if (i > 0)
					outFile << '\t';
				outFile << vec[i].name;			
			}
			outFile << '\n';
			
			double							time = startTime;
			for (size_t i = 0; i < vec[0].size(); ++i) {
				if (i > 0)
					outFile << "\n";
				if (timeStep > 0) {
					outFile << std::fixed << std::setprecision(5) << std::setw(10) << time << '\t';
					time += timeStep;
				}
				for (size_t j = 0; j < vec.size(); ++j) {
					if (j > 0)
						outFile << '\t';
					outFile << std::scientific << std::setw(10) << vec[j][i];
					//outFile << std::fixed << std::setprecision(5) << std::setw(10) << vec[j][i];
				}
			}
		} else {
			throw std::runtime_error("saving vectors to .column files only for the time being");
		}
	}


    /// load vectors from *.txt file (whitespace separated values)
	template<class Vec>
	void importVector(const Vec& vec, const string& filename) {
		string							ext = filenameExtension(filename);
		
		if (ext == ".txt") {
			// tab separated values (single line)
			ofstream						inFile(filename.c_str());
			
			for (int i = 0; i < vec.size(); ++i) 
				inFile >> vec[i];
		} else {
			throw std::runtime_error("loading vectors from .txt files only for the time being");
		}
	}


    /// weighted add of two vectors of same size (debug version contains an assertion on vector size)
	template<class Vec, class T>
	Vec weightedCombination(const Vec& v1, T weight1, const Vec& v2, T weight2) {
		assert(v1.size() == v2.size());
		
		Vec ret(v1.size());
		for (size_t i = 0; i < ret.size(); ++i) {
			ret[i] = v1[i] * weight1 + v2[i] * weight2;
		}
		return ret;
	}


    /// run a filter on a vector (result is saved inline)
	template<class Vec>
	void filterVector(Vec& v, const Vec& window) {
		Vec									temp(v.size());
		int									halfSize = window.size() / 2;
		for (size_t i = 0; i < v.size(); ++i) {
			temp[i] = 0;
			double winSum = 0;
			for (size_t j = 0; j < window.size(); ++j) {
				if ((i + j - halfSize >= 0) && (i + j - halfSize < v.size())) {
					temp[i] += window[j] * v[i + j - halfSize];
					winSum += window[j];
				}
			}
			if (winSum > 0)
				temp[i] /= winSum;
		}
		v.swap(temp);
	}


	// ***********************************************************************************************
	// struct NeighbourStruct
	//
	/// Used in class Neighbourhood. Holds relative position of single neighbour and a direction 
	/// vector pointing from reference point to neighbour.
	///
	struct NeighbourStruct {
		typedef Pattern::Vector<int, 3> RelativeVec;
		typedef Pattern::Vector<double, 3> DirectionVec;
		
		// difference between neighbour and observed cell (neighbour pos. - cell pos.)
		RelativeVec							dif;
		// direction = normalized difference
		DirectionVec						dir;
		
		NeighbourStruct(const RelativeVec& df, const DirectionVec& dr) : dif(df), dir(dr) {}
		NeighbourStruct(const RelativeVec& d) : dif(d) {
			double							len = sqrt((double)sqrLength(dif));
			dir[0] = (double)dif[0] / len;
			dir[1] = (double)dif[1] / len;
			dir[2] = (double)dif[2] / len;
		}
	};


	// ***********************************************************************************************
	// class Neighbourhood
	//
	/// Defines a neighbourhood of any given point in space (as a vector of relative positions of 
	/// neighbours to the reference point). 
	///
	class Neighbourhood {
	public:
		typedef Pattern::Vector<int, 3> RelativeVec;
		
	protected:
		std::vector<NeighbourStruct>		neighbours;
			
	public:
        /// callback function for 4-neighbourhood
		void callback4N(const RelativeVec& target, const RelativeVec& active) {
			RelativeVec						d = abs(active - target);
			int								minD = min(d[1], d[2]);
			int 								maxD = max(d[1], d[2]);
			if ((minD == 1) && (maxD == 1) && (d[0] == 0)) {
				neighbours.push_back(active - target);
			}			
		}	
		
		/// callback function for 8-neighbourhood
		void callback8N(const RelativeVec& target, const RelativeVec& active) {
			RelativeVec						d = abs(active - target);
			int 								maxD = max(d[1], d[2]);
			if ((maxD == 1) && (d[0] == 0)) {
				neighbours.push_back(active - target);
			}			
		}
		
		/// callback function for 3D equivalent of 4-neighbourhood
		void callback3x2N(const RelativeVec& target, const RelativeVec& active) {
			RelativeVec						d = abs(active - target);
			int								minD = min(min(d[0], d[1]), d[2]);
			int 								maxD = max(max(d[0], d[1]), d[2]);
			if ((minD == 1) && (maxD == 1)) {
				neighbours.push_back(active - target);
			}			
		}
		
		/// callback function for 3D equivalent of 8-neighbourhood
		void callbackCube(const RelativeVec& target, const RelativeVec& active) {
			RelativeVec						d = abs(active - target);
			int 								maxD = max(max(d[0], d[1]), d[2]);
			if (maxD == 1) {
				neighbours.push_back(active - target);
			}			
		}
		
		/// create neighbourhood type by specifying a callback function
		void create(void (Neighbourhood::*callback)(const RelativeVec&, const RelativeVec&)) {		
			RelativeVec						active;
			RelativeVec						target;
			target[0] = target[1] = target[2] = 1;		
			
			for (active[0] = 0; active[0] < 3; ++active[0]) {
				for (active[1] = 0; active[1] < 3; ++active[1]) {
					for (active[2] = 0; active[2] < 3; ++active[2]) {
						(this->*callback)(target, active);
					}
				}
			}
		}
		
		/// returns the vector containing all neighbours under the defined neighbourood
		std::vector<NeighbourStruct>& neighboursVector() {
			return neighbours;
		}
	};


	// ***********************************************************************************************
	// struct PriorityQueueEl
	//
	/// Used in calculation of excitation sequence, to create a priority queue of grid elements to
	/// excite in the next steps
	///
	struct PriorityQueueEl {
		double								time;
		MatrixIndex							index;
		
		PriorityQueueEl(double t, MatrixIndex i) : time(t), index(i) {
		}
		
		bool friend operator< (const PriorityQueueEl& p1, const PriorityQueueEl& p2) {
			return p1.time > p2.time;
		}
		
		bool friend operator== (const PriorityQueueEl& p1, const PriorityQueueEl& p2) {
			return p1.time == p2.time;
		}
	};


	// ***********************************************************************************************
	// class MeasuringPoint
	//
	/// Contains measuring point position and measurement values
	///
	class MeasuringPoint {
		std::vector<double>				measurement;
		Pattern::Vector<double, 3>		position;
		
	public:
		MeasuringPoint(const Pattern::Vector<double, 3>& pos) : position(pos) {}
		
		MeasuringPoint() {}
		
		MeasuringPoint& operator<< (double value) {
			measurement.push_back(value);
			return *this;
		}
		
		void setPosition(double x, double y, double z) {
			position[0] = z;
			position[1] = y;
			position[2] = x;
		}
		
		const Pattern::Vector<double, 3>& getPosition() const {
			return position;
		}
		
		const std::vector<double>& getMeasurement() const {
			return measurement;
		}
		
		std::vector<double>& getMeasurement() {
			return measurement;
		}
		
		Pattern::Vector<double, 3>& getPosition() {
			return position;
		}
	};


	// ***********************************************************************************************
	// class ActionPotential
	//
	/// Action potential represented with exteded Wohlfart model. Includes activation time.
	/// 
	struct ActionPotential {
		WohlfartPlus wohl;
		double at; // activation time
		
	public:
		void init(const double* wohlK, double newAt);
		void init(const ActionPotential& ap, double newAt);
		double operator() (double time) const;
		double& operator[] (size_t kIndex) {return wohl.getK()[kIndex];}
		double operator[] (size_t kIndex) const {return wohl.getK()[kIndex];}
		const double* getK() const { return wohl.getK(); }
		size_t size() const { return wohl.numParams; }
	};


    /// DEPRICATED: read (from stream) a function of type 1.9*2 + 2*3 + 4;
	bool readFunction(istringstream& oneFuncSS, double& factor, size_t& index);


    /// Helper struct used in function saveMeasurements. std::vector of saveVecElements can be supplied to function exportVectors.
	struct saveVecElement {
		string name;
		string comment;
		const std::vector<double>* data;
		
		size_t size() const {
			return data->size();
		}
		
		double operator[] (size_t n) const {
			return (*data)[n];
		}
	};
	
	
	// ***********************************************************************************************
	// class Simulation
	//
	/// holding all the data one needs to run the simulation
	/// 
	class Simulation {
	    /// RelativeVec is a vector of relative position (such as the position of one cell relative to another cell)
		typedef Pattern::Vector<int, 3> RelativeVec;
		typedef ActionPotential AP;
		typedef MeasuringPoint MP;
		
		ShapeMatrix						shape;
		ExcitationTransferMatrix		transferMatrix;
		Neighbourhood					nhood;
		Neighbourhood 					excitaionNhood;
		Region<RelativeVec>				simSpace;
		std::vector<AP>					aps;
		std::vector<MP>					mps;
		size_t 							targetNumOfAps;
		double							timeStep;
		double							timeStepInverse;
		double 							startTime;
		double 							synchronization;
		size_t							layersInExcitationSequence;
		
	public:	
		Simulation() : mps(0), targetNumOfAps(12), synchronization(1.0), 
			layersInExcitationSequence(0) 
		{}
		
	protected:
		Simulation(const Simulation&);
		void operator=(const Simulation&);
		
	public:
		/// setup time step
		void setup(double tStep);
		/// load shape matrix and set targetNumOfAps
		void loadShape(const ShapeMatrix& s);
		/// load transfer matrix throw exception on insufficient size of matrix)
		void loadTransferMatrix(const string& filename);
		/// set APs, destroying source AP vector in the process
		void setApsDestructive(std::vector<AP>& newAps);
		/// calculate exc. sequence and output it to file
		void calculateExcitationSequence(const string& output);
		/// load excitation sequence from file (can also output it to another file)
		void loadExcitationSequence(const string& input, const string& output);
		/// load file containing positions of measuring points
		void loadMeasuringPoints(const string& filename);
		/// save ECGs to file
		void saveMeasurements(const std::string& filename, const std::string& comment = "");
		/// output some settings for visual assertion
		void printSettings(const Settings& s);
		/// run the simulation, main function of Simulation class
		/// also sets ap indices for all the cells in the model
		void run(double startT, double totalTime);
		/// run the approximation of an ECG (ENDO - delayed EPI); may only be run before
		/// the real simulation, because it uses layer APs
		void runApproximation(double startT, double totalTime, double delay, std::vector<double>& result);
		
		inline size_t getTargetNumOfAps() const { return targetNumOfAps; }
		inline const std::vector<AP>& getAps() const { return aps; }
		inline double valueOfAP(double time, const ShapeElement& se) { return aps[se.apIndex](time); }
		inline Neighbourhood& neighbourhood() { return nhood; }
		inline MP& getMeasurment(size_t n) { return mps[n]; }
		inline const MP& getMeasurment(size_t n) const { return mps[n]; }
		inline size_t numMeasurements() const { return mps.size(); }
		inline const ShapeMatrix& getShape() const { return shape; }
		inline ShapeMatrix& getShape() { return shape; }
		
	private:
		/// helper for excitation sequence
		void exciteElement(priority_queue<PriorityQueueEl>& pq);
		/// transform variable \ref aps, which should hold layer aps at this point, to hold individual
        /// aps for each cell
		void setApIndices();
	};


	// ***********************************************************************************************
	// class InputLoader
	//
	/// Loads heart shape from disk or generates a default one (for debugging).
	/// @todo should be refactored into single function of Simulation.
	///
	class InputLoader {
		ShapeMatrix matrix;
		
	public:
        /// generate a heart shape (has no external parameters for tweaking the shape)
		void generateTestShape() {
			MatrixIndex matrixSize;
			matrixSize[2] = 120;
			matrixSize[1] = 160;
			matrixSize[0] = 1;
			const int dLow = 50;
			const int dHigh = 76;
			matrix.resize(matrixSize);
			
			// fill matrix with test data
			size_t mid1 = matrixSize[1] / 2, mid2 = matrixSize[2] / 4;
			for (size_t i = 0; i < matrixSize[1]; ++i) {
				for (size_t j = 0; j < matrixSize[2]; ++j) {
					double d = sqrt((double)sqr(i - mid1) + (double)sqr(j - mid2));
					if (d > dLow) {
						if (d < dHigh) {
							matrix[0][i][j] = ShapeElement((size_t)(d - dLow)/2);
						} else {
							matrix[0][i][j] = ShapeElement(0);
						}
					} else {
						matrix[0][i][j] = ShapeElement(0);
					}
				}
			}
			
			MatrixIndex								start;
			start[0] = 0;
			start[2] = matrix.size()[2] / 4;
			for (start[1] = matrix.size()[1] / 2; start[1] < matrix.size()[1]; ++start[1]) {
				if (matrix[start].layer > 0) {
					break;
				}
			}
			
			if (start[1] == matrix.size()[1]) {
				// oops, no starting point found!
				throw std::runtime_error("Could not create starting point for excitation sequence");
			}
			
			matrix[start].layer = matrix[start].layer + ShapeElement::layerStartingPoint;
			
			// export generated shape as autoGeneratedShape.txt
			exportShape(matrix, "autoGeneratedShape.matrix");
		}
		
		/// load heart shape form *.matrix file
		void loadShape(const string& filename) {
			if (filename == "")
				generateTestShape();
			else
				importMatrix(matrix, filename, &shapeLoad);
		}
		
		/// autumatic conversion operator
		operator ShapeMatrix&() {
			return matrix;
		}
	};
	
} // namespace SimLib


#endif // SIMULATOR_H_INCLUDED
