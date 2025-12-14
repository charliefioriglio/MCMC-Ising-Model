import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

def load_configuration(filename):
    """Load spin configuration from file."""
    return np.loadtxt(filename, dtype=int)

def load_data(filename):
    """Load energy and magnetization data from file."""
    data = np.loadtxt(filename, comments='#')
    sweeps = data[:, 0]
    energy = data[:, 1]
    magnetization = data[:, 2]
    return sweeps, energy, magnetization

def plot_configuration(config, kT, filename):
    """Plot spin configuration as a heatmap."""
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Create custom colormap: blue for -1, red for +1
    cmap = mcolors.ListedColormap(["#46bd31", "#e60d0d"])
    bounds = [-1.5, 0, 1.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    
    im = ax.imshow(config.T, cmap=cmap, norm=norm, origin='lower')
    
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(f'2D Ising Model Configuration at $k_BT = {kT}$', fontsize=14)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, ticks=[-1, 1])
    cbar.set_label('Spin', fontsize=12)
    
    # Calculate and display magnetization
    M = np.mean(config)
    ax.text(0.02, 0.98, f'$\\langle M \\rangle / N = {M:.4f}$', 
            transform=ax.transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Configuration plot saved to {filename}")

def plot_evolution(sweeps, energy, magnetization, kT, thermalization_sweeps, filename):
    """Plot energy and magnetization evolution."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    # Energy plot
    ax1 = axes[0]
    ax1.plot(sweeps, energy, 'b-', linewidth=0.5, alpha=0.7)
    ax1.axvline(x=thermalization_sweeps, color='r', linestyle='--', 
                label='End of thermalization')
    ax1.set_ylabel('Energy per spin $E/N$', fontsize=12)
    ax1.set_title(f'Evolution of Energy and Magnetization ($k_BT = {kT}$)', fontsize=14)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Calculate and show equilibrium average
    eq_mask = sweeps >= thermalization_sweeps
    eq_energy = energy[eq_mask]
    mean_E = np.mean(eq_energy)
    ax1.axhline(y=mean_E, color='g', linestyle='-', alpha=0.5, 
                label=f'$\\langle E/N \\rangle = {mean_E:.4f}$')
    ax1.legend(loc='upper right')
    
    # Magnetization plot
    ax2 = axes[1]
    ax2.plot(sweeps, magnetization, 'b-', linewidth=0.5, alpha=0.7)
    ax2.axvline(x=thermalization_sweeps, color='r', linestyle='--', 
                label='End of thermalization')
    ax2.set_xlabel('MCMC Sweeps', fontsize=12)
    ax2.set_ylabel('Magnetization per spin $M/N$', fontsize=12)
    ax2.grid(True, alpha=0.3)
    
    # Calculate and show equilibrium average
    eq_mag = magnetization[eq_mask]
    mean_M = np.mean(eq_mag)
    ax2.axhline(y=mean_M, color='g', linestyle='-', alpha=0.5, 
                label=f'$\\langle M/N \\rangle = {mean_M:.4f}$')
    ax2.legend(loc='upper right')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Evolution plot saved to {filename}")
    
    return mean_E, mean_M

def plot_combined_configurations(config2, config3, filename):
    """Plot both configurations side by side."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Create custom colormap: blue for -1, red for +1
    cmap = mcolors.ListedColormap(["#46bd31", "#e60d0d"])
    bounds = [-1.5, 0, 1.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    
    # kT = 2
    ax1 = axes[0]
    im1 = ax1.imshow(config2.T, cmap=cmap, norm=norm, origin='lower')
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('y', fontsize=12)
    ax1.set_title(f'$k_BT = 2$', fontsize=14)
    M2 = np.mean(config2)
    ax1.text(0.02, 0.98, f'$\\langle M \\rangle / N = {M2:.4f}$', 
             transform=ax1.transAxes, fontsize=12, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # kT = 3
    ax2 = axes[1]
    im2 = ax2.imshow(config3.T, cmap=cmap, norm=norm, origin='lower')
    ax2.set_xlabel('x', fontsize=12)
    ax2.set_ylabel('y', fontsize=12)
    ax2.set_title(f'$k_BT = 3$', fontsize=14)
    M3 = np.mean(config3)
    ax2.text(0.02, 0.98, f'$\\langle M \\rangle / N = {M3:.4f}$', 
             transform=ax2.transAxes, fontsize=12, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Add single colorbar
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im1, cax=cbar_ax, ticks=[-1, 1])
    cbar.set_label('Spin', fontsize=12)
    
    fig.suptitle('2D Ising Model Configurations in Thermal Equilibrium', fontsize=16, y=1.02)
    
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Combined configuration plot saved to {filename}")

def main():
    # Check if data files exist
    data_files = ['config_kT2.dat', 'config_kT3.dat', 'data_kT2.dat', 'data_kT3.dat']
    missing_files = [f for f in data_files if not os.path.exists(f)]
    
    if missing_files:
        print("Missing data files:", missing_files)
        print("Please run the C++ simulation first: ./ising")
        return
    
    thermalization_sweeps = 5000  # Should match the C++ simulation
    
    # Load data for kT = 2
    print("\n=== Loading data for k_B*T = 2 ===")
    config2 = load_configuration('config_kT2.dat')
    sweeps2, energy2, magnetization2 = load_data('data_kT2.dat')
    
    # Load data for kT = 3
    print("\n=== Loading data for k_B*T = 3 ===")
    config3 = load_configuration('config_kT3.dat')
    sweeps3, energy3, magnetization3 = load_data('data_kT3.dat')
    
    # Plot individual configurations
    print("\n=== Generating plots ===")
    plot_configuration(config2, 2, 'config_kT2.png')
    plot_configuration(config3, 3, 'config_kT3.png')
    
    # Plot combined configurations
    plot_combined_configurations(config2, config3, 'configurations_combined.png')
    
    # Plot evolution
    mean_E2, mean_M2 = plot_evolution(sweeps2, energy2, magnetization2, 2, 
                                       thermalization_sweeps, 'evolution_kT2.png')
    mean_E3, mean_M3 = plot_evolution(sweeps3, energy3, magnetization3, 3, 
                                       thermalization_sweeps, 'evolution_kT3.png')
    
    # Print summary
    print("\n" + "="*50)
    print("SUMMARY OF RESULTS")
    print("="*50)
    print(f"\nFor k_B*T = 2 (below critical temperature ~2.269):")
    print(f"  Mean energy per spin: <E>/N = {mean_E2:.4f}")
    print(f"  Mean magnetization per spin: <M>/N = {mean_M2:.4f}")
    
    print(f"\nFor k_B*T = 3 (above critical temperature ~2.269):")
    print(f"  Mean energy per spin: <E>/N = {mean_E3:.4f}")
    print(f"  Mean magnetization per spin: <M>/N = {mean_M3:.4f}")
    
    print("\n" + "="*50)
    print("The critical temperature for the 2D Ising model is:")
    print("  T_c = 2J/(k_B * ln(1+√2)) ≈ 2.269 J/k_B")
    print("\nAt T < T_c: System is ferromagnetic (ordered, |<M>| ≈ 1)")
    print("At T > T_c: System is paramagnetic (disordered, <M> ≈ 0)")
    print("="*50)

if __name__ == '__main__':
    main()
