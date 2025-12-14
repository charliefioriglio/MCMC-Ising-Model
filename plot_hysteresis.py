import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def load_hysteresis_data(filename):
    """Load hysteresis data from file."""
    data = np.loadtxt(filename, comments='#')
    H = data[:, 0]
    M = data[:, 1]
    
    # Split into forward and backward sweeps
    n = len(H) // 2
    H_forward = H[:n+1]
    M_forward = M[:n+1]
    H_backward = H[n+1:]
    M_backward = M[n+1:]
    
    return H_forward, M_forward, H_backward, M_backward

def plot_single_hysteresis(filename, kT, ax=None, show_arrows=True):
    """Plot a single hysteresis curve."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    
    H_f, M_f, H_b, M_b = load_hysteresis_data(filename)
    
    # Plot forward sweep (H decreasing, blue)
    ax.plot(H_f, M_f, 'b-', linewidth=1.5, label='Forward ($H$ decreasing)')
    
    # Plot backward sweep (H increasing, red)
    ax.plot(H_b, M_b, 'r-', linewidth=1.5, label='Backward ($H$ increasing)')
    
    # Add arrows to show direction
    if show_arrows:
        # Forward sweep arrow (at middle of curve)
        idx_f = len(H_f) // 2
        ax.annotate('', xy=(H_f[idx_f], M_f[idx_f]), 
                   xytext=(H_f[idx_f-5], M_f[idx_f-5]),
                   arrowprops=dict(arrowstyle='->', color='blue', lw=2))
        
        # Backward sweep arrow
        idx_b = len(H_b) // 2
        ax.annotate('', xy=(H_b[idx_b], M_b[idx_b]), 
                   xytext=(H_b[idx_b-5], M_b[idx_b-5]),
                   arrowprops=dict(arrowstyle='->', color='red', lw=2))
    
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    ax.set_xlabel(r'External Field $H$', fontsize=12)
    ax.set_ylabel(r'$\langle M \rangle / N$', fontsize=12)
    ax.set_title(f'Hysteresis Curve at $k_BT = {kT}$', fontsize=14)
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlim([-3.5, 3.5])
    ax.set_ylim([-1.1, 1.1])
    ax.grid(True, alpha=0.3)
    
    return ax

def main():
    # Find all hysteresis data files
    files = sorted(glob.glob('hysteresis_kT*.dat'))
    
    if not files:
        print("No hysteresis data files found.")
        print("Please run ising_hysteresis.cpp first")
        return
    
    # Extract temperatures from filenames
    temperatures = []
    for f in files:
        # Extract kT value from filename like "hysteresis_kT1.0.dat"
        kT_str = f.replace('hysteresis_kT', '').replace('.dat', '')
        temperatures.append(float(kT_str))
    
    # Theoretical critical temperature
    kT_c = 2.0 / np.log(1 + np.sqrt(2))

    # Plot individual hysteresis curves
    for f, kT in zip(files, temperatures):
        fig, ax = plt.subplots(figsize=(8, 6))
        plot_single_hysteresis(f, kT, ax)
        
        # Add note about T relative to T_c
        if kT < kT_c:
            note = f'$T < T_c$ (Ferromagnetic)'
        elif kT > kT_c:
            note = f'$T > T_c$ (Paramagnetic)'
        else:
            note = f'$T \\approx T_c$'
        ax.text(0.02, 0.02, note, transform=ax.transAxes, fontsize=11,
               verticalalignment='bottom', bbox=dict(boxstyle='round', 
               facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        outfile = f.replace('.dat', '.png')
        plt.savefig(outfile, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Plot saved to {outfile}")
    
    # Create comparison plot with multiple temperatures
    # Select 3 representative temperatures: below, near, and above T_c
    selected_temps = []
    selected_files = []
    
    for f, kT in zip(files, temperatures):
        if kT in [1.0, 1.5, 2.5]:  # Try to get these specific temperatures
            selected_temps.append(kT)
            selected_files.append(f)
    
    # If we don't have exactly these, just take first 3
    if len(selected_files) < 3 and len(files) >= 3:
        selected_files = files[:3]
        selected_temps = temperatures[:3]
    
    if len(selected_files) >= 3:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        for ax, f, kT in zip(axes, selected_files[:3], selected_temps[:3]):
            H_f, M_f, H_b, M_b = load_hysteresis_data(f)
            
            ax.plot(H_f, M_f, 'b-', linewidth=1.5, label='Forward')
            ax.plot(H_b, M_b, 'r-', linewidth=1.5, label='Backward')
            
            ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
            ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
            
            ax.set_xlabel(r'$H$', fontsize=12)
            ax.set_ylabel(r'$\langle M \rangle / N$', fontsize=12)
            ax.set_title(f'$k_BT = {kT}$', fontsize=14)
            ax.legend(loc='lower right', fontsize=9)
            ax.set_xlim([-3.5, 3.5])
            ax.set_ylim([-1.1, 1.1])
            ax.grid(True, alpha=0.3)
            
            # Add T < T_c or T > T_c label
            if kT < kT_c:
                ax.text(0.05, 0.95, '$T < T_c$', transform=ax.transAxes, 
                       fontsize=11, verticalalignment='top')
            else:
                ax.text(0.05, 0.95, '$T > T_c$', transform=ax.transAxes, 
                       fontsize=11, verticalalignment='top')
        
        fig.suptitle('Hysteresis Curves of 2D Ising Model at Different Temperatures', 
                    fontsize=14, y=1.02)
        plt.tight_layout()
        plt.savefig('hysteresis_comparison.png', dpi=150, bbox_inches='tight')
        plt.close()
        print("Comparison plot saved to hysteresis_comparison.png")
    
    # Create a single combined plot with all temperatures
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = plt.cm.coolwarm(np.linspace(0, 1, len(files)))
    
    for i, (f, kT) in enumerate(zip(files, temperatures)):
        H_f, M_f, H_b, M_b = load_hysteresis_data(f)
        
        # Combine for continuous loop
        H_loop = np.concatenate([H_f, H_b])
        M_loop = np.concatenate([M_f, M_b])
        
        ax.plot(H_loop, M_loop, color=colors[i], linewidth=1.5, 
               label=f'$k_BT = {kT}$')
    
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    ax.set_xlabel(r'External Field $H$', fontsize=14)
    ax.set_ylabel(r'$\langle M \rangle / N$', fontsize=14)
    ax.set_title('Hysteresis Curves at Different Temperatures', fontsize=14)
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlim([-3.5, 3.5])
    ax.set_ylim([-1.1, 1.1])
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('hysteresis_all_temps.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Combined plot saved to hysteresis_all_temps.png")

if __name__ == '__main__':
    main()
