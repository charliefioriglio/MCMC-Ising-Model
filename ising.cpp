#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include "random.h"

// Lattice dimensions
const int Nx = 100;
const int Ny = 100;
const int N = Nx * Ny;

// Coupling constant
const double J = 1.0;

// Class to represent the 2D Ising Model
class IsingModel {
public:
    std::vector<std::vector<int>> spins;  // Spin lattice: +1 or -1
    double beta;                           // Inverse temperature: 1/(k_B * T)
    double energy;                         // Total energy of the system
    double magnetization;                  // Total magnetization (sum of spins)
    random_number_generator& rng;
    
    // Constructor
    IsingModel(double kT, random_number_generator& rng_ref) : rng(rng_ref) {
        beta = 1.0 / kT;
        
        // Initialize lattice with all spins up (+1)
        spins.resize(Nx, std::vector<int>(Ny, 1));
        
        // Initial energy and magnetization
        energy = computeTotalEnergy();
        magnetization = computeTotalMagnetization();
    }
    
    // Periodic BC
    inline int periodicX(int i) const {
        return (i + Nx) % Nx;
    }
    
    inline int periodicY(int j) const {
        return (j + Ny) % Ny;
    }
    
    // Compute the total energy
    double computeTotalEnergy() const {
        double E = 0.0;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                // Only count right and down neighbors to avoid double counting
                int right = periodicX(i + 1);
                int down = periodicY(j + 1);
                E -= J * spins[i][j] * spins[right][j];
                E -= J * spins[i][j] * spins[i][down];
            }
        }
        return E;
    }
    
    double computeTotalMagnetization() const {
        double M = 0.0;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                M += spins[i][j];
            }
        }
        return M;
    }
    
    // Calculate energy change when flipping spin at (i, j)
    // ΔE = 2 * J * σ_{i,j} * (sum of neighbors)
    double computeDeltaE(int i, int j) const {
        int sum_neighbors = spins[periodicX(i+1)][j] + 
                           spins[periodicX(i-1)][j] + 
                           spins[i][periodicY(j+1)] + 
                           spins[i][periodicY(j-1)];
        return 2.0 * J * spins[i][j] * sum_neighbors;
    }
    
    // Perform one step
    bool metropolisStep() {
        // Choose a random spin
        int i = rng.uniform_int(0, Nx);
        int j = rng.uniform_int(0, Ny);
        
        // Calculate energy change
        double dE = computeDeltaE(i, j);
        
        // Accept or reject the flip
        bool accept = false;
        if (dE <= 0) {
            // Always accept if energy decreases
            accept = true;
        } else {
            // Accept with probability exp(-beta * dE)
            double prob = std::exp(-beta * dE);
            if (rng.uniform() < prob) {
                accept = true;
            }
        }
        
        if (accept) {
            // Flip spin
            spins[i][j] *= -1;
            // Update energy and magnetization
            energy += dE;
            magnetization += 2.0 * spins[i][j];
        }
        
        return accept;
    }
    
    // Perform N steps (one sweep)
    void sweep() {
        for (int k = 0; k < N; k++) {
            metropolisStep();
        }
    }
    
    double getMagnetizationPerSpin() const {
        return magnetization / N;
    }
    
    double getEnergyPerSpin() const {
        return energy / N;
    }
    
    // Save lattice config
    void saveConfiguration(const std::string& filename) const {
        std::ofstream file(filename);
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                file << spins[i][j];
                if (i < Nx - 1) file << " ";
            }
            file << "\n";
        }
        file.close();
    }
};

// Function to run simulation and save results
void runSimulation(double kT, int thermalization_steps, int measurement_steps, 
                   const std::string& config_file, const std::string& data_file,
                   random_number_generator& rng) {
    
    IsingModel model(kT, rng);
    
    // Open file for energy and magnetization data
    std::ofstream dataFile(data_file);
    dataFile << "# sweep energy_per_spin magnetization_per_spin\n";
    
    // Initial state
    dataFile << std::fixed << std::setprecision(6);
    dataFile << 0 << " " << model.getEnergyPerSpin() << " " 
             << model.getMagnetizationPerSpin() << "\n";
    
    // Thermalization
    for (int sweep = 1; sweep <= thermalization_steps; sweep++) {
        model.sweep();
        
        // Save every 10 sweeps
        if (sweep % 10 == 0) {
            dataFile << sweep << " " << model.getEnergyPerSpin() << " " 
                     << model.getMagnetizationPerSpin() << "\n";
        }
    
    }
    
    // Measurement stats
    double sum_M = 0.0;
    double sum_M2 = 0.0;
    double sum_E = 0.0;
    double sum_E2 = 0.0;
    int n_measurements = 0;
    
    for (int sweep = 1; sweep <= measurement_steps; sweep++) {
        model.sweep();
        
        // Collect measurements every sweep
        double m = model.getMagnetizationPerSpin();
        double e = model.getEnergyPerSpin();
        sum_M += m;
        sum_M2 += m * m;
        sum_E += e;
        sum_E2 += e * e;
        n_measurements++;
        
        // Save data every 10 sweeps
        if (sweep % 10 == 0) {
            dataFile << (thermalization_steps + sweep) << " " << e << " " << m << "\n";
        }
    }
    
    dataFile.close();
    
    // Averages
    double mean_M = sum_M / n_measurements;
    double mean_M2 = sum_M2 / n_measurements;
    double mean_E = sum_E / n_measurements;
    double mean_E2 = sum_E2 / n_measurements;
    
    double var_M = mean_M2 - mean_M * mean_M;
    double var_E = mean_E2 - mean_E * mean_E;
    
    double chi = model.beta * N * var_M;  // χ = β * N * Var(m)
    double Cv = model.beta * model.beta * N * var_E;  // Cv = β² * N * Var(e)
    
    // Save final configuration
    model.saveConfiguration(config_file);
}

int main(int argc, char* argv[]) {
    // Parse for seed (optional)
    uint64_t seed = 0;
    bool use_fixed_seed = false;
    
    if (argc > 1) {
        seed = std::stoull(argv[1], nullptr, 16);
        use_fixed_seed = true;
    }
    
    // Random number generator
    random_number_generator rng = use_fixed_seed ? 
        random_number_generator(seed) : random_number_generator();
    
    // Simulation parameters
    int thermalization_steps = 5000;  // Number of sweeps for thermalization
    int measurement_steps = 10000;    // Number of sweeps for measurement
    
    // Tests
    // k_B*T = 2
    runSimulation(2.0, thermalization_steps, measurement_steps, 
                  "config_kT2.dat", "data_kT2.dat", rng);
    
    // k_B*T = 3
    runSimulation(3.0, thermalization_steps, measurement_steps, 
                  "config_kT3.dat", "data_kT3.dat", rng);
    
    
    return 0;
}
