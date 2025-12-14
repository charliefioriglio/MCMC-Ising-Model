import numpy as np
import matplotlib.pyplot as plt
import os


# Input to this script: multiple data files named temperature_scan_<size>.dat
# Each file contains data for a different lattice size (e.g. now it's set to 10, 20, 35, 50, 80)

def main():

    # Load data
    # specify the size list to read (10, 20, 35, 50, 80)
    sizes = [10, 20, 35, 50, 80]
    # load all data into a dictionary
    data_dict = {}
    for size in sizes:
        filename = f'temperature_scan_{size}.dat'
        if not os.path.exists(filename):
            print(f"Missing data file: {filename}")
            print("Please run the C++ simulation ising_scan.cpp first")
            return
        data = np.loadtxt(filename, comments='#')
        data_dict[size] = data
    
    # Plot three subplots (magnetization, susceptibility, heat capacity) in one figure, with different markers for different sizes
    # every file is in the format of kT, absM, E, chi, Cv

    # Just plot one plot containing three subplots, each subplot contains data from different sizes
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    markers = ['x', 'o', 's', '^', 'd']
    for i, size in enumerate(sizes):
        data = data_dict[size]
        kT = data[:, 0]
        absM = data[:, 1]
        chi = data[:, 3]
        Cv = data[:, 4]

        # Plot Magnetization
        axes[0].plot(kT, absM, linestyle='none', marker=markers[i], markersize=4, label=f'N={size}')
        
        # Plot Susceptibility
        axes[1].plot(kT, chi, linestyle='none', marker=markers[i], markersize=4, label=f'N={size}')
        
        # Plot Heat Capacity
        axes[2].plot(kT, Cv, linestyle='none', marker=markers[i], markersize=4, label=f'N={size}')


    for ax, title in zip(axes, ['Magnetization', 'Susceptibility', 'Heat Capacity']):
        ax.set_xlabel('kT')
        ax.set_ylabel(title)
        ax.legend()
        ax.grid(True)
    plt.tight_layout()
    plt.savefig('phase_transition_comparison.png', dpi=300)



if __name__ == '__main__':
    main()
