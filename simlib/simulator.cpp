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

#include "simulator.h"
#include "timer.h"


namespace SimLib {

	void exportShape(const ShapeMatrix& shape, const string& filename, bool excit) {
		string ext = filenameExtension(filename);
		if (ext == ".matrix") {
			ofstream						outFile(filename.c_str());

			outFile << (excit ? "Excitation" : "Shape") << " file (next line specifies matrix type (2D) size (X x Y x Z) and elemet separator type (tab) after that comes data (x are rows, y are lines, z are blocks - separated with empty line)\n";
			if (shape.size()[0] > 1)
				outFile << "3D " << shape.size()[2] << " x " << shape.size()[1] << " x " << shape.size()[0];
			else
				outFile << "2D " << shape.size()[2] << " x " << shape.size()[1];
			if (excit)
				outFile << " tab\n";
			else
				outFile << "\n";

			for (int i = 0; i < (int)shape.size()[0]; ++i) {
				for (int j = 0; j < (int)shape.size()[1]; ++j) {
					for (int k = 0; k < (int)shape.size()[2]; ++k) {
						if (excit) {
							if (k > 0)
								outFile << '\t';
							outFile << std::fixed << std::setprecision(3) << shape[i][j][k].excitationDelay;
						} else {
							outFile << ((shape[i][j][k].layer & ShapeElement::layerStartingPoint) > 0 ?
								'X' : (shape[i][j][k].layer > 9 ? (char)(shape[i][j][k].layer - 10 + 'A') :
								(char)('0' + (unsigned char)shape[i][j][k].layer)));
						}
					}
					outFile << "\n";
				}
				outFile << "\n";
			}
		} else {
			throw std::runtime_error("can only export to .matrix for the time being");
		}
	}


	// parse str for formula of type 1.9*2 + 2*3 + 4;
	bool readFunction(istringstream& oneFuncSS, double& factor, size_t& index) {
		// loop through formula
		if (!oneFuncSS.eof()) {
			oneFuncSS >> factor;
			if (oneFuncSS.fail()) {
				// reading first number failed
				factor = 0.0;
				return false;
			} else {
				// first number read ok
				// read operator (skip any whitespace before)
				char							useOperator(' ');
				while ((useOperator == ' ') && !oneFuncSS.eof())
					oneFuncSS >> useOperator;
				if (useOperator == '*') {
					// expect second number (index) next
					oneFuncSS >> index;
					if (oneFuncSS.fail()) {
						index = 0;
					} else {
						// expect operator or eof
						useOperator = ' ';
						while ((useOperator == ' ') && !oneFuncSS.eof())
							oneFuncSS >> useOperator;
					}
				} else {
					if (floor(factor) != ceil(factor)) {
						index = 0;
					} else {
						index = (size_t)factor;
						factor = 1.0;
					}
				}
			}
			return true;
		} else {
			factor = 0.0;
			return false;
		}
	}


	/// ***********************************************************************************************
	/// class ActionPotential
	void ActionPotential::init(const double* wohlK, double newAt) {
		wohl.setK(wohlK);
		wohl.getK()[8] -= newAt;
		at = newAt;
	}


	void ActionPotential::init(const ActionPotential& ap, double newAt) {
		wohl.setK(ap.getK());
		wohl.getK()[8] -= newAt;
		at = newAt;
	}


	double ActionPotential::operator() (double time) const {
		return wohl[time - at];
	}


	/// ***********************************************************************************************
	/// class Simulation
	///
	void Simulation::setup(double tStep) {
		timeStep = tStep;
		timeStepInverse = 1.0 / timeStep;
	}

	void Simulation::loadShape(const ShapeMatrix& s) {
		shape = s;
		RelativeVec						zero;
		zero[0] = zero[1] = zero[2] = 0;
		simSpace.set(zero, shape.size());

		// count number of required APs
		MatrixIndex						index;
		size_t							maxLayerNum = 0;
		for (index[0] = 0; index[0] < shape.size()[0]; ++index[0]) {
			for (index[1] = 0; index[1] < shape.size()[1]; ++index[1]) {
				for (index[2] = 0; index[2] < shape.size()[2]; ++index[2]) {
					if (maxLayerNum < (shape[index].layer & ~ShapeElement::layerStartingPoint))
						maxLayerNum = (shape[index].layer & ~ShapeElement::layerStartingPoint);
				}
			}
		}
		targetNumOfAps = maxLayerNum;
	}

	void Simulation::loadTransferMatrix(const string& filename) {
		importMatrix(transferMatrix, filename);
		if ((transferMatrix.size()[0] < targetNumOfAps) || (transferMatrix.size()[1] < targetNumOfAps)) {
			throw std::runtime_error("loaded transfer matrix too small");
		}
	}

	void Simulation::setApsDestructive(std::vector<AP>& newAps) {
		aps.swap(newAps);
	}

	void Simulation::exciteElement(priority_queue<PriorityQueueEl>& pq) {
		while (!pq.empty()) {
			const PriorityQueueEl active = pq.top();
			pq.pop();
			//cerr << "work on " << active.index << "\n";

			// if element on the top of the queeu is not excited yet...
			if (shape[active.index].excitationDelay == 0) {
				// set elemet excitement delay
				shape[active.index].excitationDelay = active.time;
				//cerr << " set delay to " << active.time << "\n";

				// process its neighbours
				size_t layer = shape[active.index].layer;

				// repeat for each neighbour
				//FOREACH(const NeighbourStruct& ns, excitaionNhood.neighboursVector()) {
				for (const NeighbourStruct& ns : excitaionNhood.neighboursVector()) {
					RelativeVec			index2 = (RelativeVec)active.index - ns.dif;
					// is neighbour location within constraints?
					if (simSpace.inside(index2) && (shape[index2].layer > 0) && (shape[index2].excitationDelay == 0)) {
						// valid unexcited neighbour -> push it on queue
						if ((layer >= transferMatrix.size()[1]) || (shape[index2].layer >= transferMatrix.size()[1])) {
							std::ostringstream str;
							str << "transfer (conduction) matrix does not define layer " << std::max(layer, (shape[index2].layer));
							throw std::runtime_error(str.str());
						}
						double lag = transferMatrix[layer][shape[index2].layer];
						lag *= sqrt(sqrLength(ns.dif));
						pq.push(PriorityQueueEl(active.time + lag, index2));
					}
				}
			}
		}
	}

	void Simulation::calculateExcitationSequence(const string& output) {
		// seek starting point
		priority_queue<PriorityQueueEl> starts;
		if (shape.size()[0] > 1)
			excitaionNhood.create(&Neighbourhood::callbackCube);
		else
			excitaionNhood.create(&Neighbourhood::callback8N);
		MatrixIndex index;

		// find starting points and check number of layers
		layersInExcitationSequence = 0;
		for (index[0] = 0; index[0] < shape.size()[0]; ++index[0]) {
			for (index[1] = 0; index[1] < shape.size()[1]; ++index[1]) {
				for (index[2] = 0; index[2] < shape.size()[2]; ++index[2]) {
					if ((shape[index].layer & ShapeElement::layerStartingPoint) > 0) {
						starts.push(PriorityQueueEl(1, index));
						shape[index].layer = shape[index].layer - ShapeElement::layerStartingPoint;
					}
					if (layersInExcitationSequence < shape[index].layer)
						layersInExcitationSequence = shape[index].layer;
					//cout << (char)(shape[index].layer + '0');
				}
				//cout << "\n";
			}
			//cout << "\n";
		}

		if (starts.size() == 0) {
			// oops, no starting point found!
			throw std::runtime_error("Could not find starting point for excitation sequence");
		}

		// start excitement (recursion)
		exciteElement(starts);

		// output result
		if (output != "")
			exportShape(shape, output, true);
	}

	void Simulation::loadExcitationSequence(const string& input, const string& output) {
		// load the sequence
		ifstream inFile(input.c_str());
		string ext = filenameExtension(input);
		layersInExcitationSequence = 0;

		if (ext == ".matrix") {
			if (!inFile.is_open()) {
				throw std::runtime_error(string("could not open file ") + input + " to load excitation sequence");
			}

			typedef ShapeMatrix::SizeVector SizeVector;
			typedef ShapeMatrix::DataType DataType;

			// skip comment line
			inFile.ignore(100000000, '\n');
			SizeVector						newSize;
			for (size_t i = 0; i < SizeVector::dimension; ++i)
				newSize[i] = 1;

			// read matrix header (dimensionality)
			string							dimStr;
			size_t							dimension;
			inFile >> dimStr;
			if ( ((dimStr[1] == 'D') || (dimStr[1] == 'd')) && ((dimStr[0] > '0') &&
				(dimStr[0] <= ('0' + SizeVector::dimension))) ) {
				// dimension ok
				dimension = dimStr[0] - '0';
			} else {
				throw std::runtime_error("error while loading matrix, could not read dimensionality");
			}

			// read matrix header (size vector)
			int writeDim = shape.dimensionality() - 1;
			newSize[0] = newSize[1] = 1;
			inFile >> newSize[writeDim];
			for (--dimension; dimension > 0; --dimension) {
				inFile.ignore(100, 'x');
				--writeDim;
				inFile >> newSize[writeDim];
			}

			// read matrix header (separator type)
			string							separatorTypeStr;
			getline(inFile, separatorTypeStr, '\n');
			istringstream					separatorTypeSS(separatorTypeStr);
			string							separatorType;
			separatorTypeSS >> separatorType;
			char								separator = (separatorType == "tab" ? '\t' : (
				(separatorType == "comma") || (separatorType == ",") ? ',' :
				(separatorType == "space" ? ' ' : 0)));

			if (inFile.fail()) {
				throw std::runtime_error(string("reading file ") + input + " failed (while reading size)");
			}

			if ((shape.size()[0] != newSize[0]) || (shape.size()[1] != newSize[1]) || (shape.size()[2] != newSize[2]))
				throw std::runtime_error(string("failed to read excitation sequence due to its incompatible size"));

			for (int i = 0; i < (int)shape.size()[0]; ++i) {
				for (int j = 0; j < (int)shape.size()[1]; ++j) {
					for (int k = 0; k < (int)shape.size()[2]; ++k) {
						if ((shape[i][j][k].layer & ShapeElement::layerStartingPoint) > 0)
							shape[i][j][k].layer = shape[i][j][k].layer - ShapeElement::layerStartingPoint;
						if (layersInExcitationSequence < shape[i][j][k].layer)
							layersInExcitationSequence = shape[i][j][k].layer;
						inFile >> shape[i][j][k].excitationDelay;
						if (separator == ',')
							inFile.ignore(1);
					}
				}
			}
		} else {
			throw std::runtime_error("can only import from .matrix for the time being");
		}

		// output result
		if (output != "")
			exportShape(shape, output, true);
	}

	void Simulation::loadMeasuringPoints(const string& filename) {
		ifstream							inFile(filename.c_str());
		if (inFile.is_open())
			mps.resize(0);

		while (!inFile.eof()) {
			Pattern::Vector<double, 3>		point;
			inFile >> point[2];
			inFile.ignore(1000,',');
			inFile >> point[1];
			inFile.ignore(1000,',');
			inFile >> point[0];
			inFile.ignore(1000,'\n');
			if (inFile.fail()) {
				break;
			} else {
				mps.push_back(point);
			}
		}

		for (size_t i = 0; i < mps.size(); ++i)
			cerr << " -> " << mps[i].getPosition()[0] << ", " << mps[i].getPosition()[1] << ", " << mps[i].getPosition()[2] << "\n";

		if (mps.size() == 0) {
			throw std::runtime_error(string("failed to open file containing measuring points: ") + filename);
        }
	}
	
	void Simulation::setMeasuringPoints(const std::vector<PositionVec>& points) {
		mps.resize(points.size());
		for (size_t i = 0; i < mps.size(); ++i)
			mps[i].setPosition(points[i]);
	}
	
	void Simulation::getMeasuringPoints(std::vector<PositionVec>& points) {
		points.resize(mps.size());
		for (size_t i = 0; i < mps.size(); ++i)
			points[i] = mps[i].getPosition();
	}

	void Simulation::saveMeasurements(const std::string& filename, const std::string& comment) {
		std::vector<saveVecElement>		saveVec(mps.size());

		for (size_t i = 0; i < mps.size(); ++i) {
			std::ostringstream		nameStr;
			nameStr << mps[i].getPosition()[2] << "," << mps[i].getPosition()[1] <<
				"," << mps[i].getPosition()[0];
			saveVec[i].name = nameStr.str();
			saveVec[i].data = &(mps[i].getMeasurement());
			//mps[i].saveToFile(nameStr.str(), 1, timeStep);
		}

		if (mps.size() > 0)
			exportVectors(saveVec, filename, startTime, timeStep, comment);
	}

	void Simulation::printSettings(const Settings& s) {
		// shape
		cout << "model: " << s.inputShapeFilename << " (" << shape.size() << ")\n";

		// neighbourhood
		cout << "neighbourhood: ";
		if (s.neighbourhoodType == "2D4") {
			cout << "2D, 4 neighbours\n";
		} else if (s.neighbourhoodType == "2D8") {
			cout << "2D, 8 neighbours\n";
		} else if (s.neighbourhoodType == "3D4") {
			cout << "3D, 6 neighbours\n";
		} else if ((s.neighbourhoodType == "3D8") || (s.neighbourhoodType == "cube")) {
			cout << "3D, 26 neighbours (full cube)\n";
		}

		// simulation
		cout << "simulation time step = " << s.simulationTimeStep << "\n" << "simulation length = "
			<< s.simulationLength << "\n" << "total steps = "
			<< s.simulationLength / s.simulationTimeStep << "\n";
	}

	// new implementation, the longest loop is now innermost; much improved execution speed
	void Simulation::run(double startT, double totalTime) {
		setApIndices();

		startTime = startT;

		// create a shape with 0 border
		ShapeMatrix borderedShape;
		typedef ShapeMatrix::SizeVector SizeVector;
		borderedShape.resize(shape.size() + SizeVector(SizeVector::ConstInit(2)), ShapeMatrix::ZeroInit());
		SizeVector::ConstInit init1(1);
		SizeVector sizeVector1(init1);
		inject(shape, borderedShape, sizeVector1);

		// clear measuring points data
		for (size_t i = 0; i < mps.size(); ++i)
			mps[i].getMeasurement().resize(0);

		// we start with 0, then go up to (but not including) totalTime
		size_t totalSteps = (size_t)ceil(totalTime / timeStep);
		// run simulation only on the inside of the bordered shape, so that neighbourhood will always
		// be valid
		Region<NeighbourStruct::RelativeVec> rgn(sizeVector1, sizeVector1 + shape.size());
		typedef Pattern::Vector<double, 3> Vec3;
		std::vector<Vec3> dipole(totalSteps);
		std::vector<double> cellApValue(dipole.size());
		std::vector<std::vector<double> > measurementVectors(mps.size());
		for (size_t i = 0; i < measurementVectors.size(); ++i)
			measurementVectors[i].resize(totalSteps, 0.0);

		totalTime += startTime;

		PrecisionTimer t_cellApValue;
		PrecisionTimer t_neigApValue;
		PrecisionTimer t_measuring;
		// loop through the model cells, starting at cell (1, 1 [, 1])
		NeighbourStruct::RelativeVec	index(NeighbourStruct::RelativeVec::ConstInit(1));

		for (bool iOk = true; iOk; iOk = rgn.incIndex(index)) {

			// only use the cell if it is a part of the model (it is not empty)
			if (borderedShape[index].layer > 0) {
				// precalculate AP values for the cell
				{
					::ScopeTimer t(t_cellApValue);
					double simTime = startTime;
					for (size_t i = 0; i<cellApValue.size(); ++i) {
						cellApValue[i] = valueOfAP(simTime, borderedShape[index]);
						simTime += timeStep;
					}
				}

				// calculate vector D (dipole vector of cell)

				// start by defining a 3D vector, set to (0, 0, 0)
				Vec3 D(Pattern::Vector<double, 3>::ConstInit(0));
				std::fill(dipole.begin(), dipole.end(), D);

				// loop through cell neighbours
				//FOREACH(const NeighbourStruct& ns, nhood.neighboursVector()) {
                for (const NeighbourStruct& ns : nhood.neighboursVector()) {
					::ScopeTimer t(t_neigApValue);
					// only use neighbours that are not empty
					NeighbourStruct::RelativeVec		index2 = index - ns.dif;
					if (borderedShape[index2].layer > 0) {

						// calculate dipole time series
						size_t next = 0;
						for (double simTime = startTime; simTime < totalTime; simTime += timeStep) {
							double valueDif = valueOfAP(simTime, borderedShape[index2]) - cellApValue[next];
							dipole[next] += (ns.dif * valueDif);
							++next;
						}
					}
				}

				// have dipole vector, calculate scalar product with measuring points and
				// append it to their values
				::ScopeTimer t(t_measuring);
				for (size_t i = 0; i < mps.size(); ++i) {
					Pattern::Vector<double, 3> mPos = mps[i].getPosition() -
						(Pattern::Vector<double, 3>)index;
					mPos *= 1.0 / pow3Length(mPos);
					for (size_t t = 0; t < dipole.size(); ++t) {
						measurementVectors[i][t] += mPos * dipole[t];
					}
				}
			}
		}
//		std::cout << " - cell AP value " << t_cellApValue.totalSeconds() << "\n";
//		std::cout << " - neighbours' AP values " << t_neigApValue.totalSeconds() << "\n";
//		std::cout << " - measuring points " << t_measuring.totalSeconds() << "\n";

		// push measuring point values
		for (size_t i = 0; i < mps.size(); ++i) {
			for (size_t t = 0; t < dipole.size(); ++t) {
				mps[i] << measurementVectors[i][t];
			}
		}
	}

	void Simulation::runApproximation(double startT, double totalTime, double delay, std::vector<double>& result) {
		result.resize((size_t)(totalTime / timeStep), 0.0);

		for (size_t i = 0; i < result.size(); ++i) {
			double simTime = startT + i*timeStep;
			result[i] = aps[0](simTime + delay) - aps.back()(simTime);
		}
	}

	void Simulation::setApIndices() {
		assert(targetNumOfAps == aps.size());
		assert(layersInExcitationSequence = targetNumOfAps);

		/// build indices
		// map of translations from activation time into ap index for all the layers
		typedef std::map<double, size_t> IndexMap;
		std::vector<IndexMap> indexMap(targetNumOfAps);

		size_t nextIndex = 0;

		MatrixIndex						index;
		for (index[0] = 0; index[0] < shape.size()[0]; ++index[0]) {
			for (index[1] = 0; index[1] < shape.size()[1]; ++index[1]) {
				for (index[2] = 0; index[2] < shape.size()[2]; ++index[2]) {
					ShapeElement& s = shape[index];
					if (s.layer == 0)
						continue;

					IndexMap::iterator it = indexMap[s.layer-1].find(s.excitationDelay);
					if (it != indexMap[s.layer-1].end()) {
						s.apIndex = it->second;
					} else {
						s.apIndex = nextIndex;
						indexMap[s.layer-1][s.excitationDelay] = nextIndex;
						++nextIndex;
					}
				}
			}
		}

		/// map indices to APs
		std::vector<AP>	newAps(nextIndex);

		// loop through layers
		for (size_t apLayer = 0; apLayer < indexMap.size(); ++apLayer) {
			// for each layer, use base AP and all possible delays to create all possible cellular APs
			IndexMap::const_iterator mapI;
			for (mapI = indexMap[apLayer].begin(); mapI != indexMap[apLayer].end(); ++mapI) {
				newAps[mapI->second].init(aps[apLayer], mapI->first);
			}
		}

		/*std::ofstream kFile("C:/Cygwin/home/mdepolli/k.txt");
		for (size_t j = 0; j < aps.size(); ++j) {
			for (size_t i=0; i<9; ++i)
				kFile << aps[j].getK()[i] << "\t";
			kFile << "\n";
		}//*/

		newAps.swap(aps);
		//std::cout << "total of " << aps.size() << " different aps found\n";

		/*
		std::ofstream kFile("C:/Cygwin/home/mdepolli/k.txt");
		for (size_t j = 0; j < aps.size(); j += 70) {
			for (size_t i=0; i<9; ++i)
				kFile << aps[j].getK()[i] << "\t";
			kFile << "\n";
		}//*/
	}

} // namespace SimLib
