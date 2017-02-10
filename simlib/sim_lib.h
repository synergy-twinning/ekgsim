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
//		StreamRedirector redirect;
		Settings settings;
		Simulation sim;
		
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
