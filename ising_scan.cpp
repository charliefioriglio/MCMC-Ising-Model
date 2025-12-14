#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include "random.h"

// Lattice dimensions
const int Nx = 50;
const int Ny = 50;
const int N = Nx * Ny;

// Coupling constant
const double J = 1.0;

// Class for 2D Ising Model
class IsingModel {
public:
    std::vector<std::vector<int>> spins;
    double beta;
    double energy;
    double magnetization;
    random_number_generator& rng;
    
    // Precomputed exponentials for efficiency
    double exp_table[5];  // For dE = -8, -4, 0, 4, 8 (indices 0-4)
    
    IsingModel(double kT, random_number_generator& rng_ref) : rng(rng_ref) {
        beta = 1.0 / kT;
        
        // Initialize lattice with all spins up
        spins.resize(Nx, std::vector<int>(Ny, 1));
        
        // Precompute exponentials for Metropolis algorithm
        // dE can only be -8J, -4J, 0, 4J, or 8J for 2D Ising
        for (int i = 0; i < 5; i++) {
            int dE = -8 + 4 * i;  // -8, -4, 0, 4, 8
            exp_table[i] = std::exp(-beta * dE * J);
        }
        
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
    
    // Fast Metropolis step using precomputed exponentials
    inline void metropolisStep() {
        int i = rng.uniform_int(0, Nx);
        int j = rng.uniform_int(0, Ny);
        
        int sum_neighbors = spins[periodicX(i+1)][j] + 
                           spins[periodicX(i-1)][j] + 
                           spins[i][periodicY(j+1)] + 
                           spins[i][periodicY(j-1)];
        
        // dE = 2 * J * σ * Σneighbors
        // Since σ = ±1 and Σneighbors ∈ {-4,-2,0,2,4}
        // dE ∈ {-8,-4,0,4,8} * J
        int dE_index = (spins[i][j] * sum_neighbors + 4) / 2;  // Maps to 0,1,2,3,4
        
        double dE = 2.0 * J * spins[i][j] * sum_neighbors;
        
        bool accept = false;
        if (dE <= 0) {
            accept = true;
        } else {
            // Use precomputed exponential
            if (rng.uniform() < exp_table[dE_index]) {
                accept = true;
            }
        }
        
        if (accept) {
            spins[i][j] *= -1;
            energy += dE;
            magnetization += 2.0 * spins[i][j];
        }
    }
    
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

// Structure for measurement results
struct Results {
    double kT;
    double mean_absM;
    double mean_E;
    double chi;        // Magnetic susceptibility
    double Cv;         // Heat capacity
    double var_M;
    double var_E;
};

// Run simulation at a single temperature
Results runAtTemperature(double kT, int thermalization_sweeps, int measurement_sweeps,
                         random_number_generator& rng) {
    IsingModel model(kT, rng);
    
    // Thermalization
    for (int sweep = 0; sweep < thermalization_sweeps; sweep++) {
        model.sweep();
    }
    
    // Measurement phase
    double sum_M = 0.0, sum_M2 = 0.0, sum_absM = 0.0;
    double sum_E = 0.0, sum_E2 = 0.0;
    
    for (int sweep = 0; sweep < measurement_sweeps; sweep++) {
        model.sweep();
        
        double m = model.getMagnetizationPerSpin();
        double e = model.getEnergyPerSpin();
        
        sum_M += m;
        sum_M2 += m * m;
        sum_absM += std::abs(m);
        sum_E += e;
        sum_E2 += e * e;
    }
    
    // Calculate averages
    Results res;
    res.kT = kT;
    res.mean_absM = sum_absM / measurement_sweeps;
    res.mean_E = sum_E / measurement_sweeps;
    
    double mean_M = sum_M / measurement_sweeps;
    double mean_M2 = sum_M2 / measurement_sweeps;
    double mean_E = sum_E / measurement_sweeps;
    double mean_E2 = sum_E2 / measurement_sweeps;
    
    res.var_M = mean_M2 - mean_M * mean_M;
    res.var_E = mean_E2 - mean_E * mean_E;
    
    // χ = β * N * Var(m) where m = M/N
    // For extensive quantity: χ = β * <(ΔM)²> = β * N² * Var(m)
    // Per spin: χ/N = β * N * Var(m)
    res.chi = model.beta * N * res.var_M;
    
    // Cv = β² * N * Var(e) where e = E/N
    // Per spin: Cv/N = β² * N * Var(e)
    res.Cv = model.beta * model.beta * N * res.var_E;
    
    return res;
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
    
    // Temperature range
    double kT_min = 0.0;
    double kT_max = 4.0;
    int n_temperatures = 500;  // Number of temperature points
    
    // MCMC parameters
    // Use adaptive sweep counts based on the temperature. For temperatures close to T_c, we need more sweeps.
    // For low temperatures, fewer thermalizationsweeps are needed.
    int thermalization_sweeps_low = 2000;
    int thermalization_sweeps_tc = 25000;
    int thermalization_sweeps_high = 5000;

    int measurement_sweeps = 10000;     // ~12.5 million steps for measurement
    
    // Theoretical T_c
    double kT_c = 2.0 * J / std::log(1.0 + std::sqrt(2.0));

    // Output file
    std::ofstream outFile("temperature_scan.dat");
    outFile << "# k_B*T  |<M>|/N  <E>/N  chi/N  Cv/N  Var(M)  Var(E)\n";
    outFile << std::fixed << std::setprecision(6);
    
    std::vector<Results> all_results;
    
    for (int i = 0; i < n_temperatures; i++) {
        double kT = kT_min + (kT_max - kT_min) * i / (n_temperatures - 1);
        
        // Determine thermalization sweeps based on temperature
        int thermalization_sweeps;
        if (kT < kT_c - 0.5) {
            thermalization_sweeps = thermalization_sweeps_low;
        } else if (kT > kT_c + 0.5) {
            thermalization_sweeps = thermalization_sweeps_high;
        } else {
            thermalization_sweeps = thermalization_sweeps_tc;
        }
        Results res = runAtTemperature(kT, thermalization_sweeps, measurement_sweeps, rng);
        all_results.push_back(res);
        
        outFile << res.kT << " " << res.mean_absM << " " << res.mean_E << " "
                << res.chi << " " << res.Cv << " " << res.var_M << " " << res.var_E << "\n";
        
        std::cout << "T = " << res.kT 
                  << ", <|M|>/N = " << res.mean_absM 
                  << ", <E>/N = " << res.mean_E 
                  << ", χ/N = " << res.chi 
                  << ", Cv/N = " << res.Cv 
                  << ", Var(M) = " << res.var_M 
                  << ", Var(E) = " << res.var_E 
                  << std::endl;
    }
    
    outFile.close();
    
    // Find approximate T_c from peaks
    double max_chi = 0, kT_chi_peak = 0;
    double max_Cv = 0, kT_Cv_peak = 0;
    
    for (const auto& res : all_results) {
        if (res.chi > max_chi) {
            max_chi = res.chi;
            kT_chi_peak = res.kT;
        }
        if (res.Cv > max_Cv) {
            max_Cv = res.Cv;
            kT_Cv_peak = res.kT;
        }
    }
    
    return 0;
}
