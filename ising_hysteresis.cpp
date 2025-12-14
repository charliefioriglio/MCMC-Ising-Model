#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>
#include "random.h"

// Lattice dimensions
const int Nx = 50;
const int Ny = 50;
const int N = Nx * Ny;

// Coupling constant
const double J = 1.0;

// Ising Model with external field
class IsingModelWithField {
public:
    std::vector<std::vector<int>> spins;
    double beta;
    double H;              // External magnetic field
    double energy;
    double magnetization;
    random_number_generator& rng;
    
    IsingModelWithField(double kT, double H_init, random_number_generator& rng_ref, int init_spin = 1) 
        : rng(rng_ref) {
        beta = 1.0 / kT;
        H = H_init;
        
        // Initialize lattice with specified spin direction
        spins.resize(Nx, std::vector<int>(Ny, init_spin));
        
        energy = computeTotalEnergy();
        magnetization = computeTotalMagnetization();
    }
    
    inline int periodicX(int i) const {
        return (i + Nx) % Nx;
    }
    
    inline int periodicY(int j) const {
        return (j + Ny) % Ny;
    }
    
    double computeTotalEnergy() const {
        double E = 0.0;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                // Neighbor interaction (only right and down to avoid double counting)
                int right = periodicX(i + 1);
                int down = periodicY(j + 1);
                E -= J * spins[i][j] * spins[right][j];
                E -= J * spins[i][j] * spins[i][down];
                // External field term
                E -= H * spins[i][j];
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
    
    // Set external field and update energy
    void setField(double H_new) {
        energy -= (H_new - H) * magnetization;
        H = H_new;
    }
    
    // Calculate energy change when flipping spin at (i, j)
    // ΔE = 2Jσ(Σneighbors) + 2Hσ
    double computeDeltaE(int i, int j) const {
        int sum_neighbors = spins[periodicX(i+1)][j] + 
                           spins[periodicX(i-1)][j] + 
                           spins[i][periodicY(j+1)] + 
                           spins[i][periodicY(j-1)];
        return 2.0 * J * spins[i][j] * sum_neighbors + 2.0 * H * spins[i][j];
    }
    
    // Metropolis step
    inline void metropolisStep() {
        int i = rng.uniform_int(0, Nx);
        int j = rng.uniform_int(0, Ny);
        
        double dE = computeDeltaE(i, j);
        
        bool accept = false;
        if (dE <= 0) {
            accept = true;
        } else {
            if (rng.uniform() < std::exp(-beta * dE)) {
                accept = true;
            }
        }
        
        if (accept) {
            spins[i][j] *= -1;
            energy += dE;
            magnetization += 2.0 * spins[i][j];
        }
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
};

// Run hysteresis simulation at a given temperature
void runHysteresis(double kT, double H_max, int n_H_steps, int sweeps_per_H,
                   int measurement_sweeps, const std::string& filename,
                   random_number_generator& rng) {
    
    // Start with all spins up (+1) and strong positive field
    IsingModelWithField model(kT, H_max, rng, 1);
    
    // Thermalize at initial H
    for (int i = 0; i < 1000; i++) {
        model.sweep();
    }
    
    std::ofstream outFile(filename);
    outFile << "# H  M/N  (forward sweep then backward sweep)\n";
    outFile << std::fixed << std::setprecision(6);
    
    // Forward sweep: H from +H_max to -H_max
    for (int step = 0; step <= n_H_steps; step++) {
        double H = H_max - 2.0 * H_max * step / n_H_steps;
        model.setField(H);
        
        // Equilibrate at this H value
        for (int i = 0; i < sweeps_per_H; i++) {
            model.sweep();
        }
        
        // Measure average magnetization
        double sum_M = 0.0;
        for (int i = 0; i < measurement_sweeps; i++) {
            model.sweep();
            sum_M += model.getMagnetizationPerSpin();
        }
        double mean_M = sum_M / measurement_sweeps;
        
        outFile << H << " " << mean_M << "\n";
    }
    
    // Backward sweep: H from -H_max to +H_max
    for (int step = 1; step <= n_H_steps; step++) {
        double H = -H_max + 2.0 * H_max * step / n_H_steps;
        model.setField(H);
        
        // Equilibrate at this H value
        for (int i = 0; i < sweeps_per_H; i++) {
            model.sweep();
        }
        
        // Measure average magnetization
        double sum_M = 0.0;
        for (int i = 0; i < measurement_sweeps; i++) {
            model.sweep();
            sum_M += model.getMagnetizationPerSpin();
        }
        double mean_M = sum_M / measurement_sweeps;
        
        outFile << H << " " << mean_M << "\n";
    }
    
    outFile.close();
}

int main(int argc, char* argv[]) {
    uint64_t seed = 0;
    bool use_fixed_seed = false;
    
    if (argc > 1) {
        seed = std::stoull(argv[1], nullptr, 16);
        use_fixed_seed = true;
    }
    
    random_number_generator rng = use_fixed_seed ? 
        random_number_generator(seed) : random_number_generator();
    
    // Theoretical critical temperature
    double kT_c = 2.0 * J / std::log(1.0 + std::sqrt(2.0));
    
    // Hysteresis parameters
    double H_max = 3.0;
    int n_H_steps = 100;        // Number of H steps in each direction
    int sweeps_per_H = 100;     // Equilibration sweeps at each H
    int measurement_sweeps = 200; // Measurement sweeps at each H

    // Run at different temperatures
    // Below T_c (ferromagnetic, strong hysteresis expected)
    runHysteresis(1.0, H_max, n_H_steps, sweeps_per_H, measurement_sweeps, 
                  "hysteresis_kT1.0.dat", rng);
    
    runHysteresis(1.5, H_max, n_H_steps, sweeps_per_H, measurement_sweeps, 
                  "hysteresis_kT1.5.dat", rng);
    
    runHysteresis(2.0, H_max, n_H_steps, sweeps_per_H, measurement_sweeps, 
                  "hysteresis_kT2.0.dat", rng);
    
    // Near T_c
    runHysteresis(2.27, H_max, n_H_steps, sweeps_per_H, measurement_sweeps, 
                  "hysteresis_kT2.27.dat", rng);
    
    // Above T_c (paramagnetic, minimal hysteresis expected)
    runHysteresis(2.5, H_max, n_H_steps, sweeps_per_H, measurement_sweeps, 
                  "hysteresis_kT2.5.dat", rng);
    
    runHysteresis(3.0, H_max, n_H_steps, sweeps_per_H, measurement_sweeps, 
                  "hysteresis_kT3.0.dat", rng);
                  
    return 0;
}
