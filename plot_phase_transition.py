import numpy as np
import matplotlib.pyplot as plt
import os

def main():
    # Check if data file exists
    if not os.path.exists('temperature_scan_10.dat'):
        print("Missing data file: temperature_scan.dat")
        print("Please run the C++ simulation ising_scan.cpp first")
        return
    
    # Load data
    data = np.loadtxt('temperature_scan_10.dat', comments='#')
    kT = data[:, 0]
    absM = data[:, 1]
    E = data[:, 2]
    chi = data[:, 3]
    Cv = data[:, 4]
    
    # Theoretical critical temperature
    J = 1.0
    kT_c = 2.0 * J / np.log(1 + np.sqrt(2))
    
    # Find peaks (need to remove the 0K point for both chi and Cv since they are nans there)
    chi_peak_idx = np.argmax(chi[1:]) + 1
    Cv_peak_idx = np.argmax(Cv[1:]) + 1
    kT_chi_peak = kT[chi_peak_idx]
    kT_Cv_peak = kT[Cv_peak_idx]
    
    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot 1: Magnetization
    ax1 = axes[0]
    ax1.plot(kT, absM, linestyle='none', marker='x', color='b', markersize=3)
    ax1.axvline(x=kT_c, color='r', linestyle='--', linewidth=1.5, label=f'$T_c$ (theory) = {kT_c:.3f}')
    ax1.set_xlabel(r'$k_B T / J$', fontsize=14)
    ax1.set_ylabel(r'$|\langle M \rangle| / N$', fontsize=14)
    ax1.set_title('Magnetization', fontsize=14)
    ax1.set_xlim([kT.min(), kT.max()])
    ax1.set_ylim([0, 1.05])
    ax1.legend(loc='upper right', fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Susceptibility
    ax2 = axes[1]
    ax2.plot(kT, chi, linestyle='none', marker='x', color='b', markersize=3)
    ax2.axvline(x=kT_c, color='r', linestyle='--', linewidth=1.5, label=f'$T_c$ (theory) = {kT_c:.3f}')
    ax2.axvline(x=kT_chi_peak, color='g', linestyle=':', linewidth=1.5, label=f'Peak at {kT_chi_peak:.3f}')
    ax2.set_xlabel(r'$k_B T / J$', fontsize=14)
    ax2.set_ylabel(r'$\chi / N$', fontsize=14)
    ax2.set_title('Magnetic Susceptibility', fontsize=14)
    ax2.set_xlim([kT.min(), kT.max()])
    ax2.legend(loc='upper right', fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Heat Capacity
    ax3 = axes[2]
    ax3.plot(kT, Cv, linestyle='none', marker='x', color='b', markersize=3)
    ax3.axvline(x=kT_c, color='r', linestyle='--', linewidth=1.5, label=f'$T_c$ (theory) = {kT_c:.3f}')
    ax3.axvline(x=kT_Cv_peak, color='g', linestyle=':', linewidth=1.5, label=f'Peak at {kT_Cv_peak:.3f}')
    ax3.set_xlabel(r'$k_B T / J$', fontsize=14)
    ax3.set_ylabel(r'$C_v / N$', fontsize=14)
    ax3.set_title('Heat Capacity', fontsize=14)
    ax3.set_xlim([kT.min(), kT.max()])
    ax3.legend(loc='upper right', fontsize=11)
    ax3.grid(True, alpha=0.3)
    
    plt.suptitle('Phase Transition of the 2D Ising Model', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig('phase_transition.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\nPlot saved to phase_transition.png")
    
    # Also create individual plots for the report
    # Magnetization plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(kT, absM, linestyle='none', marker='x', color='b', markersize=3)
    ax.axvline(x=kT_c, color='r', linestyle='--', linewidth=2, label=f'$k_BT_c$ (theory) = {kT_c:.3f}')
    ax.set_xlabel(r'$k_B T / J$', fontsize=14)
    ax.set_ylabel(r'$|\langle M \rangle| / N$', fontsize=14)
    ax.set_title('Mean Magnetization vs Temperature', fontsize=14)
    ax.set_xlim([kT.min(), kT.max()])
    ax.set_ylim([0, 1.05])
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('magnetization_vs_T.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Plot saved to magnetization_vs_T.png")
    
    # Susceptibility plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(kT, chi, linestyle='none', marker='x', color='b', markersize=3)
    ax.axvline(x=kT_c, color='r', linestyle='--', linewidth=2, label=f'$k_BT_c$ (theory) = {kT_c:.3f}')
    ax.axvline(x=kT_chi_peak, color='g', linestyle=':', linewidth=2, label=f'Measured peak: {kT_chi_peak:.3f}')
    ax.set_xlabel(r'$k_B T / J$', fontsize=14)
    ax.set_ylabel(r'$\chi / N$', fontsize=14)
    ax.set_title('Magnetic Susceptibility vs Temperature', fontsize=14)
    ax.set_xlim([kT.min(), kT.max()])
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('susceptibility_vs_T.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Plot saved to susceptibility_vs_T.png")
    
    # Heat capacity plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(kT, Cv, linestyle='none', marker='x', color='b', markersize=3)
    ax.axvline(x=kT_c, color='r', linestyle='--', linewidth=2, label=f'$k_BT_c$ (theory) = {kT_c:.3f}')
    ax.axvline(x=kT_Cv_peak, color='g', linestyle=':', linewidth=2, label=f'Measured peak: {kT_Cv_peak:.3f}')
    ax.set_xlabel(r'$k_B T / J$', fontsize=14)
    ax.set_ylabel(r'$C_v / N$', fontsize=14)
    ax.set_title('Heat Capacity vs Temperature', fontsize=14)
    ax.set_xlim([kT.min(), kT.max()])
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('heat_capacity_vs_T.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Plot saved to heat_capacity_vs_T.png")
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY OF PHASE TRANSITION ANALYSIS")
    print("="*60)
    print(f"\nTheoretical critical temperature: k_B*T_c = {kT_c:.4f}")
    print(f"\nFrom susceptibility peak: k_B*T = {kT_chi_peak:.4f}")
    print(f"  Deviation from theory: {abs(kT_chi_peak - kT_c):.4f} ({100*abs(kT_chi_peak - kT_c)/kT_c:.2f}%)")
    print(f"\nFrom heat capacity peak: k_B*T = {kT_Cv_peak:.4f}")
    print(f"  Deviation from theory: {abs(kT_Cv_peak - kT_c):.4f} ({100*abs(kT_Cv_peak - kT_c)/kT_c:.2f}%)")
    print("="*60)

if __name__ == '__main__':
    main()
