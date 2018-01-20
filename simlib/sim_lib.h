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

#ifndef SIM_LIB_H_INCLUDED
#define SIM_LIB_H_INCLUDED

/**
    wrapper of simulator functions into a library
**/

#include "simulator.h"


namespace SimLib {

	/// ***********************************************************************************************
	/// class EkgSim
	/// 
	/// to use EkgSim, use:
	///     EkgSim ekgSime;
	///     ekgSim.loadSettings();
	///     ekgSim.applySettings();
	///     ekgSim.loadShape();             % heart model shape (contains a 3D matrix of integer values, representing layers; negative numbers represent excitation starting points)
	///     ekgSim.loadTransferMatrix();    % transfer matrix containing conduction velocities between layers
	///     ekgSim.loadMeasuringPoints();   % measuring points are triplets defining 3D measuring position relative to the heart model origin (which equals the corner which is read first in the shape file)
	///     ekgSim.simExcitationSequence(); % either loads or calculates the excitation sequence - beware, heart model must contain at least one starting point
	///     ekgSim.run();                   % run the simulation as it has been set up
	///     std::vector = ekgSim.getMeasurement(n); or
    ///         ekgSim.saveMeasurement(comment); 
    ///         ekgSim.saveMeasurementAs(filename, comment);    % get and store measurements - results of the simulation
	/// 
	class EkgSim {
	public:
		typedef Simulation::PositionVec PositionVec;
	private:
//		StreamRedirector redirect;
		Settings settings;
		Simulation sim;
		// measuring positions as loaded from file
		std::vector<PositionVec> originalMeasuringPositions;
		// measuring positions as given to the optimizer
		std::vector<PositionVec> measuringPositions;
		// u, v pairs for each measuring positions (vectors u and v define the plane of movement for measuring point)
		std::vector<PositionVec> vectorU, vectorV;
		
	public:
		EkgSim() /*: redirect(cerr)*/ {
		}
		
		void loadSettings(const char *fname = "settings.ini") {
			ScopeTimerForLogging tm(std::cerr, "loading settings                            ");
			Ini::File ini(fname);
			if (!ini.ok() || ini.notFound()) {
				if (ini.empty())
					std::cout << "warning, " << fname << " is empty. ";
				else if (ini.notFound())
					std::cout << "warning, " << fname << " not found. ";
				else 
					std::cout << "warning, unknown error while opening " << fname << ". ";
			} else {
				settings.loadFromIni(ini);
			}
		}
		
		Settings& getSettings() {
			return settings;
		}
		
		/// applies the selected timestep and neighbourhood to the simulator
		void applySettings() {
			// define neighbourhood
			if (settings.neighbourhoodType == "2D4") {
				sim.neighbourhood().create(&Neighbourhood::callback4N);
			} else if (settings.neighbourhoodType == "2D8") {
				sim.neighbourhood().create(&Neighbourhood::callback8N);
			} else if (settings.neighbourhoodType == "3D4") {
				sim.neighbourhood().create(&Neighbourhood::callback3x2N);
			} else if ((settings.neighbourhoodType == "3D8") || (settings.neighbourhoodType == "cube")) {
				sim.neighbourhood().create(&Neighbourhood::callbackCube);
			} else {
				throw std::runtime_error("\n   unknown neighbourhood (only know of these: 2D4, 3D4, 2D8, 3D8 (cube))");
			}
			
			// set time step
			sim.setup(settings.simulationTimeStep);
		}
		
		/// Set APs (a destructive operation for the parameter, whose value is cleared after the operation).
		/// Destructive operation is used for its speed.
		void setApsDestructive(std::vector<ActionPotential>& destructible) {
			sim.setApsDestructive(destructible);
		}
		
		/// get APs (not destructable, APs may only be read, not modified)
		const std::vector<ActionPotential>& getAps() const {
			return sim.getAps();
		}
		
		/// Load transfer matrix - a matrix of excitation delays (the inverse of transfer/excitation speeds).
		void loadTransferMatrix() {
			ScopeTimerForLogging tm(std::cerr, "loading transfer matrix                     ");
			sim.loadTransferMatrix(settings.inputTransferFilename);
		}
		
		/// load the vectors defining measuring points (locations in which we wish to measure ECGs)
		void loadMeasuringPoints() {
			ScopeTimerForLogging tm(std::cerr, "loading measuring points                    ");
			cerr << "measuring points: \n";
			sim.loadMeasuringPoints(settings.inputPointsFilename);
			// store local copy of measuring points and calculate their uvw space 
			// uvw space is the space defined by a plane perpendicular to the position vector of an individual measuring point
			sim.getMeasuringPoints(originalMeasuringPositions);
			measuringPositions.resize(originalMeasuringPositions.size());
			vectorU.resize(originalMeasuringPositions.size());
			vectorV.resize(originalMeasuringPositions.size());
			// note that since measuring points are provided as (z,y,x), the code below also uses such coordinate order
			PositionVec x;
			x[0] = 0; x[1] = 0; x[2] = 1;
			for (size_t i = 0; i < originalMeasuringPositions.size(); ++i) {
				PositionVec &p = originalMeasuringPositions[i];
				measuringPositions[i] = p;
				vectorU[i] = cross(x, p);
				vectorV[i] = cross(vectorU[i], p);
				// since the coordinates are in reverse order, vector product is negative
				vectorU[i] *= -1;
				// but vector v is already ok, since it was made from negative vector u ...
				normalize(vectorU[i]);
				normalize(vectorV[i]);
				
				std::cout << " u" << i << " = " << vectorU[i][2] << "," << vectorU[i][1] << "," << vectorU[i][0] << "\n";
				std::cout << " v" << i << " = " << vectorV[i][2] << "," << vectorV[i][1] << "," << vectorV[i][0] << "\n";
			}
		}
		
		/// move the measuring points from original position by requested relative amount, in uvw space
		void moveMeasuringPoints(const std::vector<PositionVec>& uvwDisplacement) {
			if (uvwDisplacement.size() == measuringPositions.size()) {
				for (size_t i = 0; i < originalMeasuringPositions.size(); ++i) {
					// measuring position is displaced by u*vectorU +v*vectorV + w*w; w is always zero
					measuringPositions[i] = originalMeasuringPositions[i] + uvwDisplacement[i][0] * vectorU[i] + uvwDisplacement[i][1] * vectorV[i]; 
				}
				sim.setMeasuringPoints(measuringPositions);
			} else {
				std::cerr << uvwDisplacement.size() << " != " << measuringPositions.size() << " " << originalMeasuringPositions.size() << std::endl;
				throw std::runtime_error("invalid vector of measuring point displacements");
			}
		}
		
		/// load heart shape
		void loadShape() {
			ScopeTimerForLogging tm(std::cerr, "loading shape                               ");
			InputLoader	tempInLoader;
			tempInLoader.loadShape(settings.inputShapeFilename);
			sim.loadShape(tempInLoader);
		}
		
		/// determines how many APs are required (as many as there are defined layers in the heart model)
		size_t requiredAps() const {
			return sim.getTargetNumOfAps();
		}
		
		/// simulate the excitation sequence
		void simExcitationSequence() {
			if (settings.inputExcitationSequenceFilename == "") {
				ScopeTimerForLogging		tm(std::cerr, "calculating excitation sequence             ");
				sim.calculateExcitationSequence(settings.outputExcitationSequence);
			} else {
				ScopeTimerForLogging		tm(std::cerr, "loading excitation sequence                 ");
				sim.loadExcitationSequence(settings.inputExcitationSequenceFilename, settings.outputExcitationSequence);
			}
			
			// generate APs and ap indices in shape matrix
		}
		
		void printSettings() {
			sim.printSettings(settings);
		}
		
		void printMessages(std::ostream& out) {
			/*out << "\n" << redirect.messageBuffer.str() << std::endl;*/
		}
		
		/// run full simulation
		void run() {
			sim.run(settings.simulationStart, settings.simulationLength);
		}
		
		/// run approximate simulation (string model)
		void runApproximation(double delay, std::vector<double>& result) {
			sim.runApproximation(settings.simulationStart, settings.simulationLength, delay, result);
		}
		
		/// returns the specified measurement (simulation result)
		inline const std::vector<double>& getMeasurement(size_t n) const {
			return sim.getMeasurment(n).getMeasurement();
		}
		
		/// returns the number of defined measuremnts (simulation results)
		inline size_t numMeasurements() const {
			return sim.numMeasurements();
		}
		
		/// saves a single result of the simulation
		void saveMeasurement(const std::string& comment = "") {
			ScopeTimerForLogging tm(std::cerr, "saving measurement                          ");
			sim.saveMeasurements(settings.outputFilename, comment);
		}
		
		/// saves all simulation results into a single file (*.column)
		void saveMeasurementAs(const char* fname, const std::string& comment = "") {
			ScopeTimerForLogging tm(std::cerr, "saving measurement                          ");
			sim.saveMeasurements(fname, comment);
		}
		
		const ShapeMatrix& getModelShape() const {
		    return sim.getShape();
		}
		
		ShapeMatrix& getModelShape() {
		    return sim.getShape();
		}

	private:
		EkgSim(const EkgSim&);
		void operator=(const EkgSim&);
	};

} // namespace SimLib;


typedef SimLib::EkgSim EkgSim;

#endif // SIM_LIB_H_INCLUDED
